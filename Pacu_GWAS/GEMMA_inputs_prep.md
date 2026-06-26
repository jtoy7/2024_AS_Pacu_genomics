Preparation of input files for GEMMA GWAS analyses
================
**AUTHOR:** Jason A. Toy  
**DATE:** 2026-06-25 <br><br>


Required GEMMA input files:
1. Genotype file
    - bed/bim/fam format (PLINK1)
    - Use variants with MAF > 0.05 (not LD-pruned)
2. Phenotype file
    - This can be combined with PLINK1 files by encoding phenotype in column 6 of .fam file
3. Kinship matrix
    - Can be generated in GEMMA as a preliminary step
    - Generate this using the genotype files with only the individuals used for each analysis (i.e., full roster/clone-pruned)
    - Will need to use a different set of variants for this, since this relatedness measure should be calculated using the MAF > 0.05, **LD-pruned** variant set to avoid influence of large linkage blocks
4. Covariate file
    - Can create this from the LMM dataset files
    - Must contain a column of 1s if one wants to include an intercept term
    - Location may need to be coded as multiple binary variables (with one omitted as reference)

Note: Sample order must match exactly across the .fam file, covariate file, and kinship matrix
    
<br>

To create these input files, we first need a list of individuals to keep for each analysis:
- Full roster of ramets with complete data (P. acuta only; start with `Pacu_complete_multivariate_dataset`)
- Clone-pruned list (P. acuta only; start with `Pacu_complete_multivariate_dataset_clonepruned`)

<br>

### Start by uploading complete multivariate datasets to the cluster workspace:
```bash
cd ~/winhome/OneDrive - Old Dominion University/BarshisLab/FieldWork/2024-012_AmSamoa/jt/StressExperiments/AllPAMResults

cp ./Pacu_complete_multivariate_dataset* jtoy@wahab.hpc.odu.edu:/archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/starting files
```
This uploads the two files: `Pacu_complete_multivariate_dataset` and `Pacu_complete_multivariate_dataset_clonepruned`.

<br>

### Find the relevant PLINK files for the kinship matrix calculation
We want the dataset that's filtered for MAF > 0.05 and LD-pruned. So it looks like we want the `pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes` PLINK files (bim/bed/fam). These are in the directory, `/archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf`.

There is a version of the .bim file with the suffix, `_original`. This is the one we want. The current file named `pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.bim` was altered to change the chromosome names to a format accepted by ADMIXTURE. So first we want to change it back to the original version:
```bash
# First copy the current bim file to note that it's for use with ADMIXTURE
cp pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.bim pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_admixture.bim

# Then create a copy of the `_original` file that has the simplified name to match the .bed and .fam files (this will overwrite the previous one used for ADMIXTURE analysis)
cp pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_original.bim pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.bim
```

Now let's check a few things about the file:
```bash
# check sample name format used in plink files
head pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.fam
```
```
0       2024_ALOF_Pspp_X1_1     0       0       0       -9
0       2024_ALOF_Pspp_X2_1     0       0       0       -9
0       2024_ALOF_Pspp_X4_1     0       0       0       -9
0       2024_ALOF_Pver_01_1     0       0       0       -9
0       2024_ALOF_Pver_02_1     0       0       0       -9
0       2024_ALOF_Pver_03_1     0       0       0       -9
0       2024_ALOF_Pver_04_1     0       0       0       -9
0       2024_ALOF_Pver_05_1     0       0       0       -9
0       2024_ALOF_Pver_06_1     0       0       0       -9
0       2024_ALOF_Pver_07_1     0       0       0       -9
```
So this is the naming convention that I need to use for my sample lists.


```bash
# Check number of samples in this dataset
wc -l pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.fam
```
```
396 pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.fam
```
So we've confirmed this includes all samples, including extras/non-P.acuta samples and samples with missing data.

<br>

Let's copy the PLINK files for this dataset to our new `gemma` directory to start in a clean workspace:
```bash
mkdir ../../gemma_gwas/starting_files/

cp pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.bim ../../gemma_gwas/starting_files/
cp pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.bed ../../gemma_gwas/starting_files/
cp pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes.fam ../../gemma_gwas/starting_files/
```

<br>

### Start creating gemma input files
Now let's create our sample lists from the multivariate data files:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/starting_files

# First double-check lengths of files (note that this will include the header line)
wc -l Pacu_complete_multivariate_dataset*
```
```
   133 Pacu_complete_multivariate_dataset_clonepruned.tsv
   367 Pacu_complete_multivariate_dataset.tsv
```

```bash
# Now isolate the sample names that match the PLINK .fam file from each multivariate data file and save as sample lists (while removing header)
mkdir ../sample_lists

cut -f14 Pacu_complete_multivariate_dataset.tsv | tail -n +2 > ../sample_lists/sample_list_allramet.txt
cut -f14 Pacu_complete_multivariate_dataset_clonepruned.tsv | tail -n +2 > ../sample_lists/sample_list_cpruned.txt

wc -l ../sample_lists/*.txt
```
```
 366 sample_list_allramet.txt
 132 sample_list_cpruned.txt
```
