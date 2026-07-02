Run GWAS analyses with GEMMA
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2026-07-01 <br><br>


## Start with all-ramet dataset (366 ramets)
### Location-only covariate model
Start with a model that only has the one 'location' covariate.

<br>

`run_gemma_allramet_location.slurm`:
```bash
#!/bin/bash

#SBATCH --job-name=run_gemma_allramet_location_2026-07-01
#SBATCH --output=%A_%a_%x.out
#SBATCH --error=%A_%a_%x.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jtoy@odu.edu
#SBATCH --partition=main
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=5-00:00:00
#SBATCH --cpus-per-task=1

set -euo pipefail


# Load GEMMA module
module load container_env gemma/0.98.5

# Define paths
BASEDIR=/archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas

# Output directories
OUTPREFIX='allramet_location'
OUTBASE=${BASEDIR}/gemma_results/lmm/${OUTPREFIX}
mkdir -p "${OUTBASE}"


# Dense, non-LD-pruned, MAF>0.05 PLINK prefixes for association testing
BFILE=${BASEDIR}/genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet

# Kinship matrices generated from subsetted, LD-pruned, MAF>0.05 PLINK files
KIN=${BASEDIR}/kinship_matrix/output/pacu_allramet_gk1_kinship.cXX.txt

# Covariate files: no header, same sample order as .fam, first column = intercept
COV=${BASEDIR}/covariates/pacu_allramet_location.cov

# Match this to your upstream max missingness threshold.
# GEMMA has its own SNP filters, so setting this prevents GEMMA's default missingness behavior from silently imposing a stricter filter than intended.
MISS_GUARD=0.2
MAF_GUARD=0.05

# Move to working directory
cd "${OUTBASE}"

# Run GEMMA LMM
# -----------------------------
crun.gemma gemma \
  -bfile "${BFILE}" \
  -n 1 \
  -k "${KIN}" \
  -c "${COV}" \
  -maf "${MAF_GUARD}" \
  -miss "${MISS_GUARD}" \
  -lmm 4 \
  -outdir "${OUTBASE}" \
  -o "${OUTPREFIX}_lmm"

echo "Finished ${OUTPREFIX} GEMMA run"
```
