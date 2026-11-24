#!/usr/bin/env python3
"""Simulate and re-infer a one-population demographic model and gamma DFE.

One invocation = one replicate, written to its own TSV file (one header
row + one data row). Failed replicates are written too (status = failed)
so the failure rate can be reported rather than silently dropped.

Changes from the previous version:
  1. The nonsynonymous SFS is simulated under the TRUE demography; the DFE
     is then inferred under the FITTED demography, so demographic error
     carries into DFE error as it would in a real analysis.
  2. theta_NS used for DFE inference is estimated from the synonymous data
     (fitted theta_syn x 2.31), not fixed at its true value.
  3. Both fits use several random starting points (not perturbations of
     the true values) and keep the best log-likelihood.
  4. Fits that end at a parameter bound are flagged.

The dadi sample size (ns) is the number of sampled haplotypes.
"""

import argparse
import csv
import os
import sys
import time

import dadi
import dadi.DFE as DFE
import nlopt
import numpy as np


TRUE_ALPHA = 0.2
TRUE_BETA = 1000.0
N_ANCESTRAL = 1000
MU = 1e-8
L_NEUTRAL = 1000e6
NS_SYN_RATIO = 2.31
THETA_0 = 4 * N_ANCESTRAL * MU * L_NEUTRAL
THETA_NS = THETA_0 * NS_SYN_RATIO

DEMOG_LOWER = np.array([1e-5, 1e-5])
DEMOG_UPPER = np.array([2.0, 2.0])
SEL_LOWER = np.array([1e-3, 1e-2])
SEL_UPPER = np.array([1.0, 50000.0])

# Random starting points are drawn log-uniformly from these ranges
DEMOG_START_RANGE = (np.array([1e-3, 1e-3]), np.array([1.0, 1.0]))
SEL_START_RANGE = (np.array([0.02, 10.0]), np.array([1.0, 10000.0]))

GAMMA_BOUNDS = (1e-4, 2000)
GAMMA_PTS = 100

HEADER = [
    "sample_haplotypes", "replicate", "seed", "status", "error_message",
    "nu_true", "T_true", "nu_fit", "T_fit",
    "theta_syn_true", "theta_syn_fit", "ll_demog", "demog_starts_ok", "demog_at_bound",
    "alpha_true", "beta_true", "alpha_fit", "beta_fit",
    "theta_ns_true", "theta_ns_used", "ll_dfe", "dfe_starts_ok", "dfe_at_bound",
    "sec_demog_fit", "sec_caches", "sec_dfe_fit", "sec_total",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("ns", type=int, help="Number of sampled haplotypes")
    p.add_argument("nu", type=float, help="Current size / ancestral size")
    p.add_argument("T", type=float, help="Time since size change (units of 2*N_ancestral generations)")
    p.add_argument("replicate", type=int, help="Replicate number")
    p.add_argument("seed", type=int, help="Random seed")
    p.add_argument("output", help="Output TSV for this replicate (overwritten)")
    p.add_argument("--n-starts", type=int, default=10,
                   help="Random starting points per optimization (default 10)")
    p.add_argument("--cpus", type=int, default=1,
                   help="CPUs for building the selection caches (default 1)")
    return p.parse_args()


def write_row(path, values):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    row = [values.get(col, "NA") for col in HEADER]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(HEADER)
        writer.writerow(row)


def clip_negative(fs):
    """Extrapolated spectra can have tiny negative entries; set them to 0
    so Poisson sampling does not fail."""
    fs = fs.copy()
    fs[fs < 0] = 0
    return fs


def log_uniform(rng, lo, hi):
    return np.exp(rng.uniform(np.log(lo), np.log(hi)))


def at_bound(x, lo, hi, tol=0.01):
    """True if any parameter is within ~1% (log scale) of a bound."""
    lx = np.log(x)
    return int(np.any(np.abs(lx - np.log(lo)) < tol) or np.any(np.abs(lx - np.log(hi)) < tol))


def multistart(fit_one, loglik, rng, start_range, n_starts):
    """Run fit_one from n_starts random starts; return the best fit."""
    best_p, best_ll, n_ok, last_err = None, -np.inf, 0, None
    for _ in range(n_starts):
        p0 = log_uniform(rng, *start_range)
        try:
            res = fit_one(p0)
            popt = np.asarray(res[0] if isinstance(res, tuple) else res, dtype=float)
            if popt.shape != p0.shape or not np.all(np.isfinite(popt)):
                continue
            ll = float(loglik(popt))
        except Exception as exc:  # one bad start should not end the replicate
            last_err = exc
            continue
        if not np.isfinite(ll):
            continue
        n_ok += 1
        if ll > best_ll:
            best_p, best_ll = popt, ll
    if best_p is None:
        raise RuntimeError(f"All {n_starts} starts failed. Last error: {last_err}")
    return best_p, best_ll, n_ok


def main(args, out):
    if args.ns < 4:
        raise ValueError("ns must be at least 4 haplotypes")
    if args.nu <= 0 or args.T <= 0:
        raise ValueError("nu and T must be positive")

    np.random.seed(args.seed)                  # used by dadi's .sample()
    rng = np.random.default_rng(args.seed + 1)  # used for starting points

    demog_params = np.array([args.nu, args.T], dtype=float)
    selection_params = np.array([TRUE_ALPHA, TRUE_BETA], dtype=float)
    ns = [args.ns]
    pts_l = [40, 50, 60] if args.ns < 10 else [args.ns * 10, args.ns * 12, args.ns * 15]

    # --- 1. Simulate the synonymous (neutral) SFS under the true demography
    demo_model_ex = dadi.Numerics.make_extrap_func(dadi.Demographics1D.two_epoch)
    model_true = clip_negative(demo_model_ex(demog_params, ns, pts_l))
    data_sfs_dem = (model_true * THETA_0).fold().sample()

    # --- 2. Infer demography (multinomial: theta is a free scaling factor)
    def fit_demog(p0):
        return dadi.Inference.opt(
            p0, data_sfs_dem, demo_model_ex, pts_l,
            lower_bound=list(DEMOG_LOWER), upper_bound=list(DEMOG_UPPER),
            algorithm=nlopt.LN_NELDERMEAD, maxeval=10000,
        )

    def ll_demog(p):
        return dadi.Inference.ll_multinom(demo_model_ex(p, ns, pts_l).fold(), data_sfs_dem)

    t0 = time.time()
    demog_fit, demog_ll, demog_ok = multistart(
        fit_demog, ll_demog, rng, DEMOG_START_RANGE, args.n_starts)
    out["sec_demog_fit"] = round(time.time() - t0, 1)

    # Estimated synonymous theta -> theta_NS used for DFE inference
    model_fit = demo_model_ex(demog_fit, ns, pts_l).fold()
    theta_syn_fit = float(dadi.Inference.optimal_sfs_scaling(model_fit, data_sfs_dem))
    theta_ns_fit = theta_syn_fit * NS_SYN_RATIO

    out.update(nu_fit=demog_fit[0], T_fit=demog_fit[1], theta_syn_fit=theta_syn_fit,
               ll_demog=demog_ll, demog_starts_ok=demog_ok,
               demog_at_bound=at_bound(demog_fit, DEMOG_LOWER, DEMOG_UPPER))

    # --- 3. Simulate the nonsynonymous SFS under the TRUE demography
    t0 = time.time()
    cache_kwargs = dict(pts=pts_l, gamma_bounds=GAMMA_BOUNDS, gamma_pts=GAMMA_PTS,
                        verbose=False, cpus=args.cpus)
    spectra_true = DFE.Cache1D(demog_params, ns, DFE.DemogSelModels.two_epoch_sel,
                               **cache_kwargs)
    data_sfs_sel = spectra_true.integrate(
        selection_params, None, DFE.PDFs.gamma, THETA_NS, None)
    data_sfs_sel = clip_negative(data_sfs_sel).fold().sample()

    # --- 4. Infer the DFE under the FITTED demography and estimated theta_NS
    spectra_fit = DFE.Cache1D(demog_fit, ns, DFE.DemogSelModels.two_epoch_sel,
                              **cache_kwargs)
    out["sec_caches"] = round(time.time() - t0, 1)
    for name, cache in (("true", spectra_true), ("fitted", spectra_fit)):
        min_val = float(np.min(cache.spectra))
        if min_val < -1e-3:
            raise RuntimeError(f"{name}-demography cache has a negative value: {min_val}")

    def fit_sel(p0):
        return dadi.Inference.opt(
            p0, data_sfs_sel, spectra_fit.integrate, pts=None,
            func_args=[DFE.PDFs.gamma, theta_ns_fit],
            lower_bound=list(SEL_LOWER), upper_bound=list(SEL_UPPER),
            verbose=0, multinom=False,
            algorithm=nlopt.LN_NELDERMEAD, maxeval=10000,
        )

    def ll_sel(p):
        model = spectra_fit.integrate(p, None, DFE.PDFs.gamma, theta_ns_fit, None).fold()
        return dadi.Inference.ll(model, data_sfs_sel)

    t0 = time.time()
    sel_fit, sel_ll, sel_ok = multistart(
        fit_sel, ll_sel, rng, SEL_START_RANGE, args.n_starts)
    out["sec_dfe_fit"] = round(time.time() - t0, 1)

    out.update(alpha_fit=sel_fit[0], beta_fit=sel_fit[1], theta_ns_used=theta_ns_fit,
               ll_dfe=sel_ll, dfe_starts_ok=sel_ok,
               dfe_at_bound=at_bound(sel_fit, SEL_LOWER, SEL_UPPER),
               status="ok", error_message="")


if __name__ == "__main__":
    args = parse_args()
    out = dict(
        sample_haplotypes=args.ns, replicate=args.replicate, seed=args.seed,
        status="failed", error_message="",
        nu_true=args.nu, T_true=args.T, theta_syn_true=THETA_0,
        alpha_true=TRUE_ALPHA, beta_true=TRUE_BETA, theta_ns_true=THETA_NS,
    )
    start = time.time()
    try:
        main(args, out)
    except Exception as exc:
        out["error_message"] = str(exc).replace("\t", " ").replace("\n", " ")
        print(f"Replicate failed: {exc}", file=sys.stderr)
    out["sec_total"] = round(time.time() - start, 1)
    print(f"Replicate {args.replicate} (ns={args.ns}): {out['status']} in "
          f"{out['sec_total'] / 60:.1f} min", file=sys.stderr)
    # Always write a row, so failures are counted rather than dropped
    write_row(args.output, out)
