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
- Full roster of ramets with complete data (P. acuta only; start with `combined_all_complete.tsv`)
- Clone-pruned list (P. acuta only; start with `combined_cpruned_complete.tsv`)