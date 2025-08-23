import dadi
import dadi.DFE as DFE
from sys import argv
import nlopt
import sys
# input: params
# output: bias for two DFE parameters

# argv[1] => ns
# argv[2] => nu
# argv[3] => T

# generate folded SFS data
demog_params = [float(argv[2]), float(argv[3])]
selection_params = [0.2, 1000.]
ns = [int(argv[1])]
N=1000
mu=1e-8
L=1000e6
theta0 = 4*N*mu*L
theta_ns = theta0 * 2.31 # human params
if ns[0] < 10:
    pts_l = [40, 50, 60]
else:
    pts_l = [ns[0]*10, ns[0]*12, ns[0]*15]

# demographic model definition
demo_model = dadi.Demographics1D.two_epoch
demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)
data_sfs_dem = demo_model_ex(demog_params, ns, pts_l)*theta0
try:
    data_sfs_dem = data_sfs_dem.fold().sample()
except ValueError:
    sys.exit(1)

# demographic model fit
lower_bound, upper_bound = [1e-5, 1e-5], [2, 2.]
p0 = dadi.Misc.perturb_params(demog_params, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
demog_params_fit = dadi.Inference.opt(p0, data_sfs_dem, demo_model_ex, pts_l, 
                          lower_bound=lower_bound,
                          upper_bound=upper_bound, algorithm=nlopt.LN_NELDERMEAD, maxeval=10000)

# selection spectra cache
spectra = DFE.Cache1D(demog_params_fit[0], ns, DFE.DemogSelModels.two_epoch_sel, pts=pts_l, 
                      gamma_bounds=(1e-5, 500), gamma_pts=50, verbose=True, cpus=1)

data_sfs_sel = spectra.integrate(selection_params, None, DFE.PDFs.gamma, theta_ns, None)
try:
    data_sfs_sel = data_sfs_sel.fold().sample()
except ValueError:
    print("Value error with inferred selection params", selection_params)
    sys.exit(1)

# perform inference
sel_params = [0.2, 1000.]
lower_bound, upper_bound = [1e-3, 1e-2], [1, 50000.]
p0 = dadi.Misc.perturb_params(sel_params, lower_bound=lower_bound,
                              upper_bound=upper_bound)
popt = dadi.Inference.opt(p0, data_sfs_sel, spectra.integrate, pts=None,
                          func_args=[DFE.PDFs.gamma, theta_ns],
                          lower_bound=lower_bound, upper_bound=upper_bound, 
                          verbose=len(sel_params), multinom=False, algorithm=nlopt.LN_NELDERMEAD)

popt[0] - selection_params
with open("results/"+argv[1]+"_"+argv[2]+"_"+argv[3]+".txt", "a") as f:
    print(*list(popt[0] - selection_params), *list(popt[0]), *selection_params, 
          *list(demog_params_fit[0] - demog_params), *list(demog_params_fit[0]),*demog_params, file=f)
