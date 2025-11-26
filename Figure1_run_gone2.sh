#!/bin/bash
#SBATCH --job-name=EW_gone
#SBATCH --output=/scratch1/marjanak/wolf_msmc2/logs/gone.out
#SBATCH --error=/scratch1/marjanak/wolf_msmc2/logs/gone.err
#SBATCH --partition=qcb
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100G


module purge
module load gcc/13.3.0
module load bcftools

# Filter for biallelic SNPs

bcftools view -m2 -M2 -v snps EW.chr1.gvcf.gz > EW.chr1.vcf


# Run Gone2 with all 6 individuals

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -r 1 -t24 -o EW_chr1_r1

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -r 1.5 -t24 -o EW_chr1_r1.5

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -r 0.75 -t24 -o EW_chr1_r0.75


# Run Gone2 with 3 individuals

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -i 3 -r 1 -t24 -o EW_chr1_i3_r1

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -i 3 -r 1.5 -t24 -o EW_chr1_i3_r1.5

/scratch1/marjanak/wolf_msmc2/GONE2/gone2 EW.chr1.vcf -i 3 -r 0.75 -t24 -o EW_chr1_i3_r0.75