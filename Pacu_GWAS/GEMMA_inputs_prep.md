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

### Create GEMMA input files

#### Sample lists

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

<br>

Now append "_1" to the sample names to match the PLINK file names exactly, add a column of zeros to these sample lists to create the FID column required by PLINK:
```bash
cd ../sample_lists

awk 'BEGIN{OFS="\t"} {print 0, $1"_1"}' sample_list_allramet.txt | sort -k2,2 > sample_list_allramet.keep
awk 'BEGIN{OFS="\t"} {print 0, $1"_1"}' sample_list_cpruned.txt | sort -k2,2 > sample_list_cpruned.keep
```

<br>

#### Kinship matrix

Now use these keep files to filter the PLINK dataset to create two new datasets
```bash
mkdir /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/kinship_matrix
cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/sample_lists

module load plink/1.9-20240319

# create a subset bed/bim/fam with PLINK for all-ramet dataset
crun.plink plink --bfile ../starting_files/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes \
    --keep sample_list_allramet.keep \
    --make-bed \
    --out ../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet \
    --allow-extra-chr

# create a subset bed/bim/fam with PLINK for clone-pruned dataset
crun.plink plink --bfile ../starting_files/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes \
    --keep sample_list_cpruned.keep \
    --make-bed \
    --out ../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned \
    --allow-extra-chr
```
```
191515 MB RAM detected; reserving 95757 MB for main workspace.
76737 variants loaded from .bim file.
396 people (0 males, 0 females, 396 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.nosex
.
--keep: 366 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 366 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.999809.
76737 variants and 366 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.bed
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.bim
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.fam
```

```
191515 MB RAM detected; reserving 95757 MB for main workspace.
76737 variants loaded from .bim file.
396 people (0 males, 0 females, 396 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.nosex
.
--keep: 132 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 132 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.999806.
76737 variants and 132 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.bed
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.bim
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.fam
... done.
```

<br>

Even though it doesn't use phenotype info for kinship matrix calculation, GEMMA still uses the phenotype column (6) in the .fam file to determine which samples have phenotype data. If it is left as all `-9` values, it will think there are no acceptable samples and not run the calculation, so we need to add dummy values of `1` to column 6 for all rows:
```bash
# change column 6 to all 1s, saving to temp file
awk '{$6=1; print}' pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.fam > pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet_temp.fam

# replace old fam file
mv pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet_temp.fam pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.fam


# change column 6 to all 1s, saving to temp file
awk '{$6=1; print}' pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.fam > pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned_temp.fam

# replace old fam file
mv pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned_temp.fam pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned.fam
```


#### Run kinship matrix calculation

The kinship matrix can be calculated in GEMMA using two different formulas, specified by the `-gk 1` or `-gk 2` flags.
- `-gk 1` calculates the "centered" relatedness matrix
- `-gk 2` calculates the standardized relatedness matrix

Here is what the manual says about choosing between these two methods:
```
Which of the two relatedness matrix to choose will largely depend on the underlying genetic architecture of the given trait. Specifically, if SNPs with lower minor allele frequency tend to have larger effects (which is inversely proportional to its genotype variance), then the standardized genotype matrix is preferred. If the SNP effect size does not depend on its minor allele frequency, then the centered genotype matrix is preferred. In our previous experience based on a limited examples, we typically find the centered genotype matrix provides better control for population structure in lower organisms, and the two matrices seem to perform similarly in humans.
```

Given this information, we'll use the `gk -1` option.

Run separate commands for each sample set:
```bash
# load GEMMA module
module load gemma/0.98.5

cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/kinship_matrix

# All-ramet centered kinship matrix
crun.gemma gemma \
  -bfile pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet \
  -gk 1 \
  -o pacu_allramet_gk1_kinship

# Clone-pruned centered kinship matrix
crun.gemma gemma \
  -bfile pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned \
  -gk 1 \
  -o pacu_cpruned_gk1_kinship
```

Sanity-check number of columns and rows:
```bash
cd output/

awk -F'\t' '{print NF; exit}' pacu_allramet_gk1_kinship.cXX.txt
awk -F'\t' '{print NF; exit}' pacu_cpruned_gk1_kinship.cXX.txt
```
```
366
132
```
```bash
wc -l pacu_allramet_gk1_kinship.cXX.txt
wc -l pacu_cpruned_gk1_kinship.cXX.txt
```
```
366 pacu_allramet_gk1_kinship.cXX.txt
132 pacu_cpruned_gk1_kinship.cXX.txt
```

Check matrix values:
```bash
head -n 1 pacu_cpruned_gk1_kinship.cXX.txt
```
```
0.3331046075    0.1467563404    -0.0007884663961        0.00573508739   -0.009913133077 0.01050127062   0.1514684088    0.06014035357   -8.357245469e-05        -0.014067081    -0.003500508583 -0.009276607383 0.0074283106    -0.005072324836 -0.005139677284 0.002722573563-0.008593489156 -0.009024034227 -0.008805736329 0.00237475902   -0.01194846075  -0.008941188182 -0.01635688812  -0.0003801671075        -0.002938598056 -0.01050495034  -0.00109994283  -0.005349028314 -0.00741654491  -0.01209094568  -0.00216093229  -0.003488041747 0.01139302152 -0.004076413842 0.006748477049  0.0004874612054 0.002598724382  0.0008935294198 -0.004789230831 0.04035801429   -0.004877672518 -0.0006662554793        -0.003541483038 -0.00109466281  0.01096429047   -0.0089917366   -0.005117233662 -0.0003827003327        -0.005912744856       -0.00672011565  -0.01020037398  -0.01087350333  -0.01143875218  -0.006653930484 -0.008897244144 -0.006847056823 -0.005179293769 -0.01800952836  -0.01414116651  -0.007469640459 -0.00705662535  -0.007787608072 0.004663355887  -0.01609071802  -0.008015500514       0.01193139783   -0.01544327211  -0.003207651329 -0.01180935646  -0.01304580413  -0.01489821639  -0.01864614718  -0.01343906865  -0.007813679174 -0.0003259342743        0.0009369991041 0.001470236711  -0.01393127883  -0.01051821634  0.001871514497  -0.003976918868       0.0006422498605 0.007323840652  -0.002435003379 0.001740831438  -0.004414869636 -0.006302660455 -0.00476185707  -0.01274924143  -0.005314228729 -0.01064049245  0.003620819988  -0.01592394587  0.009369607079  -0.01994103472  -0.00204706666  -0.004332370291 -0.01743142238        -0.01273455955  -0.01159115415  -0.0117262198   -0.007516934037 -0.0152939517   -0.01800007294  0.005155766301  -0.02077770417  -0.01807996255  -0.005220310993 -0.008795981377 -0.01084583361  -0.01752108683  -0.006839932919 -0.01082742307  -7.827984778e-05      -0.0008150402194        -0.01644121992  0.005455194634  -0.00193357774  -0.004057212881 -0.01084906964  -0.01687463964  -0.0003590878558        -0.01125119925  -0.008216145399 -0.003381647796 -0.007919231556 -0.008676136546 -0.01350314017  0.001519034547  0.006608206188        -0.001631594911 0.000895314042
```

<br>

#### Extract phenotype (ED50) values

Phenotypes will ultimately be incorporated into the genotype plink files as column 6 of the .fam file, but first let's just create a separate file with sample names and ED50 values.

