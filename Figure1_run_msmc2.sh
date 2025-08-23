#!/bin/bash
#SBATCH --job-name=ew_msmc2_pipeline
#SBATCH --output=/scratch1/marjanak/wolf_msmc2/pipeline.out
#SBATCH --error=/scratch1/marjanak/wolf_msmc2/pipeline.err
#SBATCH --time=48:00:00
#SBATCH --partition=qcb
#SBATCH --cpus-per-task=12
#SBATCH --mem=36G

set -e  # Exit on any error

module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load bcftools

# File paths
sample_file=/project/jazlynmo_738/DataRepository/Canids/metaData/DogEWolf_CanFam3.1_Mooney2023MBE/EW_indivList_downsample_N6.txt
gvcf=/project/jazlynmo_738/DataRepository/Canids/Invariant/DogEWolf_CanFam3.1_Mooney2023MBE/final_master_allsites_6files_gvcf/MasterEVERYTHING_6_bespokefixed_v4_minRGQ1minGQ20minDP10maxHet0.99missing2Ten_3All15_5_v4_mergeGaTKfiltered_varnonvar_TenBCTMIREWAll15PGAW_jointcalled_chr1.vcf.gz
an=12

# Output directories
outdir=/scratch1/marjanak/wolf_msmc2/split_EW_chr1
mkdir -p $outdir
mkdir -p /scratch1/marjanak/wolf_msmc2/multihet
mkdir -p /scratch1/marjanak/wolf_msmc2/msmc2_chr1_rtest
mkdir -p /scratch1/marjanak/wolf_msmc2/logs

echo "=== Ethiopian Wolf MSMC2 Analysis Pipeline ==="

# Step 0: Create negative mask file
echo "Step 0: Creating negative mask file..."
zgrep -w '^chr1' /project/jazlynmo_738/DataRepository/Reference_Files/CanFam3.1/CanFam3.1_masks/canfam3.1.rep.map.genic.bed.gz > /scratch1/marjanak/wolf_msmc2/chr1_neg_mask.bed

# Step 1: Filter samples and apply AN filter
echo "Step 1: Filtering VCF data..."
bcftools view -S "$sample_file" "$gvcf" | bcftools filter -i "AN >= $an" -Oz > /scratch1/marjanak/wolf_msmc2/EW.chr1.gvcf.gz

# Step 2: Create BED file (proper 3-column BED format)
echo "Step 2: Creating mask BED file..."
bcftools query -f '%CHROM\t%POS0\t%POS\n' /scratch1/marjanak/wolf_msmc2/EW.chr1.gvcf.gz > /scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed

# Step 3: Index the resulting filtered GVCF
echo "Step 3: Indexing filtered GVCF..."
bcftools index --threads 8 -t /scratch1/marjanak/wolf_msmc2/EW.chr1.gvcf.gz

# Step 4: Split VCF by sample
echo "Step 4: Splitting VCF by sample..."
while IFS= read -r sample; do
    echo "Processing sample: $sample"
    bcftools view -s $sample -Oz -o $outdir/${sample}.EW.chr1.gvcf.gz /scratch1/marjanak/wolf_msmc2/EW.chr1.gvcf.gz
    bcftools index --threads 8 -t $outdir/${sample}.EW.chr1.gvcf.gz
done < "$sample_file"

# Step 5: Activate conda environment for MSMC2
echo "Step 5: Activating MSMC2 environment..."
eval "$(conda shell.bash hook)"
conda activate /home1/marjanak/.conda/envs/mymsmc2

# Step 6: Create input files for MSMC2 using 6 individuals (12 haplotypes)
echo "Step 6: Creating MSMC2 input for 6 individuals..."
generate_multihetsep.py \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--negative_mask=/scratch1/marjanak/wolf_msmc2/chr1_neg_mask.bed \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW2.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW3.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW4.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW7.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW9.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW10.EW.chr1.gvcf.gz > /scratch1/marjanak/wolf_msmc2/multihet/EW.chr1.multihetsep.txt

# Step 7: Run MSMC2 using 6 individuals (12 haplotypes) with different r values
echo "Step 7: Running MSMC2 for 6 individuals with different r values..."
r_values=(0.1 1 10)
for RHO in "${r_values[@]}"; do
    echo "Running MSMC2 for EW chr1 with -r = $RHO"
    msmc2_Linux -t 12 -r $RHO \
      -o /scratch1/marjanak/wolf_msmc2/msmc2_chr1_rtest/EW.chr1.r${RHO}.msmc2 \
      -I 0-1,2-3,4-5,6-7,8-9,10-11 \
      /scratch1/marjanak/wolf_msmc2/multihet/EW.chr1.multihetsep.txt
    echo "Completed MSMC2 run for -r = $RHO"
done

# Step 8: Create input files for MSMC2 using 2 individuals (4 haplotypes)
echo "Step 8: Creating MSMC2 input for 2 individuals..."
generate_multihetsep.py \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--mask=/scratch1/marjanak/wolf_msmc2/EW.chr1.mask.bed \
--negative_mask=/scratch1/marjanak/wolf_msmc2/chr1_neg_mask.bed \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW2.EW.chr1.gvcf.gz \
/scratch1/marjanak/wolf_msmc2/split_EW_chr1/EW3.EW.chr1.gvcf.gz > /scratch1/marjanak/wolf_msmc2/multihet/EW.n2.chr1.multihetsep.txt

# Step 9: Run MSMC2 using 2 individuals (4 haplotypes) with different r values
echo "Step 9: Running MSMC2 for 2 individuals with different r values..."
for RHO in "${r_values[@]}"; do
    echo "Running MSMC2 for EW chr1 with -r = $RHO"
    msmc2_Linux -t 12 -r $RHO \
      -o /scratch1/marjanak/wolf_msmc2/msmc2_chr1_rtest/EW.n2.chr1.r${RHO}.msmc2 \
      -I 0-1,2-3 \
      /scratch1/marjanak/wolf_msmc2/multihet/EW.n2.chr1.multihetsep.txt
    echo "Completed MSMC2 run for -r = $RHO"
done

echo "=== Pipeline Complete ==="
echo "All results are in /scratch1/marjanak/wolf_msmc2/"