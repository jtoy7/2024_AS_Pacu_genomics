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
    --maf 0.05 \
    --geno 0.2 \
    --make-bed \
    --out ../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet \
    --allow-extra-chr

# create a subset bed/bim/fam with PLINK for clone-pruned dataset
crun.plink plink --bfile ../starting_files/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes \
    --keep sample_list_cpruned.keep \
    --maf 0.05 \
    --geno 0.2 \
    --make-bed \
    --out ../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_cpruned \
    --allow-extra-chr
```
```
385573 MB RAM detected; reserving 192786 MB for main workspace.
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
0 variants removed due to missing genotype data (--geno).
1775 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
74962 variants and 366 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.bed
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.bim
+
../kinship_matrix/pver_all_QDPSB_MISSMAF05filtered_ld_pruned_0.2_genotypes_gemmakinship_allramet.fam
... done.
```

```
385573 MB RAM detected; reserving 192786 MB for main workspace.
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
0 variants removed due to missing genotype data (--geno).
5470 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
71267 variants and 132 people pass filters and QC.
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
cd ../kinship

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

```
GEMMA 0.98.5 (2021-08-25) by Xiang Zhou, Pjotr Prins and team (C) 2012-2021
Reading Files ...
## number of total individuals = 366
## number of analyzed individuals = 366
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs/var        =    74962
## number of analyzed SNPs         =    74888
Calculating Relatedness Matrix ...
================================================== 100%
**** INFO: Done.
```

```
GEMMA 0.98.5 (2021-08-25) by Xiang Zhou, Pjotr Prins and team (C) 2012-2021
Reading Files ...
## number of total individuals = 132
## number of analyzed individuals = 132
## number of covariates = 1
## number of phenotypes = 1
## number of total SNPs/var        =    71267
## number of analyzed SNPs         =    71182
Calculating Relatedness Matrix ...
================================================== 100%
**** INFO: Done.
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
head -n 1 pacu_allramet_gk1_kinship.cXX.txt
```
```
0.330159512     -0.008517511097 0.001705728144  -0.005536317944 -0.005506641966 0.3198330552    -0.01094451341  -0.00482286193  -0.007258048335 0.3198924163    0.3195511362    0.3192312378    -0.008995532777 0.3198236157    0.3196094524    0.0006076378323 0.3195005186     -0.0004122151771        0.3192805375    -0.009100769669 -0.007639882727 -0.007086411794 -0.001303310415 -0.0087949061   -0.008817983454 -0.01063231606  -0.008181619938 -0.004782401125 0.001361587478  0.3190330872    -0.008197416228 -0.008699421346 0.0003230732171  0.3193601319    -0.001587759632 -0.0002316496607        0.3201631974    0.000415382126  -0.007438399093 -0.01135611229  -0.01135900825  -0.01581205214  -0.0156855935   -0.01567233682  -0.01568383293  -0.01590066806  -0.009105172871 -0.01582349913  -0.01002825064   -0.01561863988  -0.01607612731  -0.009267837753 -0.009216700448 -0.01155654772  -0.005467793122 -0.003097052198 -0.006926585654 -0.009301100492 -0.01547799253  -0.002503804009 -0.003293025602 -0.009083630093 -0.01590861035  -0.007230183532 -0.01553127001  -0.01536307621   -0.01552923253  -0.007365573052 -0.006515542564 -0.01558747923  -0.002442580803 -0.00709922321  -0.01110434643  -0.006795193557 -0.005328379393 -0.01617416256  -0.01040728306  -0.006979800424 -0.01048298885  -0.01111130796  -0.01006348922  -0.004917265756  -0.02214251417  -0.01188555487  -0.01065043844  -0.02224885268  -0.02231562112  -0.01126695624  -0.004767847188 -0.02183839544  -0.02234097285  -0.01098294336  -0.0111903052   -0.02235766483  -0.02294728899  -0.005572805811 -0.02244488967  -0.005719843691 -0.02239007162   -0.02278201078  -0.02286069516  -0.02212629888  -0.005577470768 -0.006531918887 -0.003327391552 -0.007454771663 -0.004635059226 -0.005588265366 -0.007357249982 0.005762120559  0.005100073739  -0.0226356376   -0.01095551368  -0.01460569187  -0.0161090637    -0.00322973702  -0.01464512999  -0.01492045424  -0.01588637212  -0.01595747305  -0.01595624699  -0.01690125422  -0.01620904922  -0.01679515589  -0.0163910543   -0.01688169315  -0.0166893198   -0.01647743723  -0.007723146942 -0.01678798711  -0.01676502526  -0.01628675673   -0.009721939108 -0.01555804519  -0.003201882741 -0.01616431881  -0.006100949987 -0.01487370391  -0.01644654117  -0.01647391546  -0.01631030786  -0.01365842075  -0.01382553095  -0.01492727775  -0.01500446211  -0.01461529794  -0.006040670243 -0.01480978251   -0.01373323937  -0.01472956604  -0.01484122797  -0.01479911161  -0.01599580504  -0.01246593948  -0.009008892206 -0.01250993453  -0.01847593942  -0.01798555815  -0.008857631123 -0.008507085317 -0.004969132307 -0.0157258086   -0.004370796516 -0.007468488216 -0.003889143499  -0.001418097212 -0.01196291139  -0.001476086809 -0.005448028348 -0.007241514292 -0.01284668133  -0.00179698672  -0.01258729188  -0.005564689204 -0.01677440224  -0.004677855048 -0.005862485541 -0.01611638309  -0.005349148995 -0.001663067442 -0.008111249606  -0.004732266823 -0.009469418077 -0.003482917282 -0.003335560235 -0.005526043952 -0.009966627237 -0.005452354159 -0.002452266371 -0.01233085406  -0.01323518327  -0.01886711095  -0.007308799226 -0.006580130042 -0.006340771927 -0.009875930014 -0.006350085017 -0.004792123558  -0.008124146311 -0.0066292935   -0.006285612828 -0.008835941531 -0.008395421348 0.00121403089   -0.00786738171  -0.01099736077  -0.003632722157 -0.01071414261  -0.01068721618  -0.01024038779  -0.004902069498 -0.03196193616  -0.03177375878  -0.03208533989   -0.03132938657  -0.03153986608  -0.03150874125  -0.03160036605  -0.03147846654  -0.0109572295   -0.0316240171   -0.01046090727  -0.01043167955  -0.03160621869  -0.0316665058   -0.03133338369  -0.03186423984  -0.01120372465  -0.01132564075  -0.03158515774  -0.03204846543   -0.01488079643  -0.03137791646  -0.0114572864   -0.01460420204  -0.03191956411  -0.01444386746  -0.01439014954  -0.006889417853 -0.006723342923 -0.009010002308 -0.0146668521   -0.01452442859  -0.01302927816  -0.01459298851  -0.01024020935  -0.01437263363   -0.01442021809  -0.01152542603  -0.01108447148  0.004939392544  0.005577911253  0.0003381294458 -0.01149619125  0.002851987734  -0.0001695698154        0.0006030837182 0.005611324259  0.00351679602   -0.01091083501  -0.006096143708 -0.00568273806  -0.003325999586  -0.01141726742  -0.006113192957 -0.005719786165 -0.006310087322 -0.00578670239  9.00150452e-05  0.000234748484  -0.0001691253378        0.00492076083   0.005253080113  0.0009644669136 0.0003784934692 0.0005069682931 0.000212831603  -0.007222146638 -0.0003354180809 0.0003959132557 0.001101396026  -0.01139072986  0.0005790703991 -0.01043772464  -3.453802833e-05        -0.007716632611 0.009855058778  -0.01296012854  -0.006120507453 -0.01365745554  -0.01307852583  -0.01313116585  -0.00040194422  -1.685430201e-05        -0.0003064645377 0.0003834779377 -0.009799794388 -0.01344381605  -0.01333858716  -0.01347416369  0.00016900295   -0.0001811763414        -0.01313861543  -0.01390955527  0.0001073658628 -0.01319320323  -0.01328151542  -0.01294727282  -0.01304090657  -0.009998975151 -0.006008239557  -0.005855006611 -0.00647416383  0.001948880867  -0.01275758119  0.001653156501  -0.01265663804  -0.01314438382  -0.0136016169   -0.01304712309  -0.01311008104  -0.01272501533  -0.01292355255  -0.01258474777  -0.01340864281  -0.01311388825  -0.01307921813   0.001886116076  -0.01479388157  -0.01449778385  -0.001647796595 -0.01491643582  -0.01163371593  -0.002406760305 -0.007294480848 -0.02232609756  -0.007412488538 -0.0192043379   -0.01908388399  -0.00603004776  -0.02293218626  -0.0215975281   -0.002572866626 -0.001774185259  -0.02284959972  -0.01631907252  -0.01262716891  -0.007818455927 -0.006543034121 -0.02201154548  -0.01134633719  -0.01233951613  -0.01922548601  -0.01942113178  -0.02218045127  -0.02242278061  -0.01705996432  -0.02277860993  -0.02234908527  -0.02214882611   -0.02198576551  -0.02231590919  -0.01912578283  -0.01875959063  -0.002383899949 -0.02186212796
```
```bash
head -n 1 pacu_cpruned_gk1_kinship.cXX.txt
```
```
0.3505563974    0.1540979463    -0.0003646042459        0.005884475869  -0.01041461316  0.01127301119   0.1613493435    0.06361349237   -9.024814757e-05        -0.01490564618  -0.003706341264 -0.009959531896 0.007991108736  -0.00545981687  -0.005215901486 0.00293605876    -0.008663304937 -0.009408758014 -0.009198206329 0.002683406429  -0.01289384598  -0.009372225921 -0.01731984188  -0.0005005163824        -0.002846764577 -0.0110269177   -0.001093345009 -0.005810517565 -0.00768050475  -0.01292029089  -0.002007041544 -0.003666754292  0.01184389173   -0.004410974299 0.007148056074  0.0007089711705 0.003057848614  0.0009017207712 -0.005133401352 0.04151732073   -0.004778848324 -0.0007017643392        -0.003594631522 -0.001248366555 0.01159537654   -0.009313520719 -0.005618706319 -0.0001367098657 -0.006127532327 -0.007229438552 -0.01108306885  -0.01160266434  -0.0120922397   -0.007589426624 -0.009421568128 -0.007037237332 -0.005495578974 -0.01914331254  -0.01494364425  -0.007845800842 -0.007685291455 -0.008150756295 0.004868727321  -0.01697069055  -0.00843726986   0.01253244158   -0.01644527014  -0.003449558066 -0.01248462514  -0.01384410347  -0.01584897188  -0.01990813232  -0.01398031316  -0.008272488926 -0.0001611381795        0.001302913743  0.001540284287  -0.01487316109  -0.01106918052  0.002192942499  -0.004101310138  0.0009115213662 0.007708122879  -0.002745720619 0.001889767185  -0.004186835197 -0.006997995934 -0.004880214535 -0.01341396285  -0.00525920027  -0.01133461948  0.003148707193  -0.01726314372  0.01037381301   -0.0211052247   -0.00214047197  -0.00442150414   -0.01836780177  -0.01370037308  -0.01228237133  -0.01248730892  -0.007796664091 -0.01625431369  -0.01937686223  0.005543325613  -0.0219926367   -0.01925285982  -0.005718741438 -0.009228146312 -0.01129729325  -0.01868245203  -0.00716208103  -0.01130404297  0.0003675666901  -0.0007738850549        -0.01748626799  0.005455524936  -0.002072704552 -0.004310470395 -0.01154640007  -0.01771006289  -8.364057549e-05        -0.01211928784  -0.008966715406 -0.00349116149  -0.008413773679 -0.009075606828 -0.0144475281   0.001951255648   0.007367550234  -0.001126756192 0.00072044873
```

<br>

#### Extract phenotype (ED50) values

Phenotypes will ultimately be incorporated into the genotype plink files as column 6 of the .fam file, but first let's just create a separate file with sample names and ED50 values.

Extract sample name and ED50 values from multivariate dataset:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/starting_files
mkdir ../phenotypes

cut -f9,14 Pacu_complete_multivariate_dataset.tsv | awk 'BEGIN{OFS="\t"} {print $2"_1", $1}' | tail -n +2 | sort -k1,1 > ../phenotypes/Pacu_complete_ed50s.txt
cut -f9,14 Pacu_complete_multivariate_dataset_clonepruned.tsv | awk 'BEGIN{OFS="\t"} {print $2"_1", $1}' | tail -n +2 | sort -k1,1 > ../phenotypes/Pacu_complete_cpruned_ed50s.txt
```
```bash
cd ../phenotypes/
head *
```
```
==> Pacu_complete_cpruned_ed50s.txt <==
2024_ALOF_Pver_03_1     37.7280845350833
2024_ALOF_Pver_06_1     37.7970495591192
2024_ALOF_Pver_07_1     37.6121588898274
2024_ALOF_Pver_08_1     37.6306961285469
2024_ALOF_Pver_09_1     36.9392839584712
2024_ALOF_Pver_10_1     37.5853252077299
2024_ALOF_Pver_19_1     37.5161130960978
2024_ALOF_Pver_23_1     37.9794565750908
2024_ALOF_Pver_24_1     37.8598758999653
2024_ALOF_Pver_27_1     38.046780127904

==> Pacu_complete_ed50s.txt <==
2024_ALOF_Pver_02_1     37.274612836513
2024_ALOF_Pver_03_1     37.7280845350833
2024_ALOF_Pver_04_1     37.6325184444698
2024_ALOF_Pver_05_1     37.871523866322804
2024_ALOF_Pver_06_1     37.7970495591192
2024_ALOF_Pver_07_1     37.6121588898274
2024_ALOF_Pver_08_1     37.6306961285469
2024_ALOF_Pver_09_1     36.9392839584712
2024_ALOF_Pver_10_1     37.5853252077299
2024_ALOF_Pver_11_1     37.4656074378858
```

#### Prepare PLINK files for genotype + phenotype data

For the genotype files, we want to use the **non-LD-pruned**, MAF > 0.05 dataset.
We still want to use the MAF > 0.05 cutoff because very low-frequency alleles introduce a lot of uncertainty. The genotype calls themselves become uncertain, and even when the genotype calls are correct, any association at a very rare variant is supported by only a few minor allele copies, so the estimated effect is fragile and highly sensitive to individual samples, sample-specific covariate patterns, missingness, or genotyping error.

So the dataset we want is `pver_all_QDPSB_MISSMAF05filtered_genotypes`
This dataset has 1,122,273 SNPs.

It is currently in PLINK2 format (pgen, psam, pvar), so we need to convert it to PLINK1 format:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf

module load plink/2.0-2026.01.10

crun.plink plink2 \
  --pfile pver_all_QDPSB_MISSMAF05filtered_genotypes \
  --make-bed \
  --out pver_all_QDPSB_MISSMAF05filtered_genotypes
```

<br>

Now use the keep files created previously to filter the PLINK dataset to create two new datasets. This includes commands to refilter based on MAF and variant missingness:
```bash
mkdir /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/genotypes
cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/sample_lists

module load plink/1.9-20240319

# create a subset bed/bim/fam with PLINK for all-ramet dataset
crun.plink plink --bfile /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf/pver_all_QDPSB_MISSMAF05filtered_genotypes \
    --keep sample_list_allramet.keep \
    --maf 0.05 \
    --geno 0.2 \
    --make-bed \
    --out ../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet \
    --allow-extra-chr

# create a subset bed/bim/fam with PLINK for clone-pruned dataset
crun.plink plink --bfile /archive/barshis/barshislab/jtoy/pver_gwas/hologenome_mapped_all/vcf/pver_all_QDPSB_MISSMAF05filtered_genotypes \
    --keep sample_list_cpruned.keep \
    --maf 0.05 \
    --geno 0.2 \
    --make-bed \
    --out ../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned \
    --allow-extra-chr
```
```
385573 MB RAM detected; reserving 192786 MB for main workspace.
1122273 variants loaded from .bim file.
396 people (0 males, 0 females, 396 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.nosex .
--keep: 366 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 366 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.999976.
1 variant removed due to missing genotype data (--geno).
14404 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
1107868 variants and 366 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.bed +
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.bim +
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.fam ... done.
```
```
385573 MB RAM detected; reserving 192786 MB for main workspace.
1122273 variants loaded from .bim file.
396 people (0 males, 0 females, 396 ambiguous) loaded from .fam.
Ambiguous sex IDs written to
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.nosex .
--keep: 132 people remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 132 founders and 0 nonfounders present.
Calculating allele frequencies... done.
Total genotyping rate in remaining samples is 0.999976.
2 variants removed due to missing genotype data (--geno).
66744 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
1055527 variants and 132 people pass filters and QC.
Note: No phenotypes present.
--make-bed to
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.bed +
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.bim +
../genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.fam ... done.
```
So the all-ramet and clone-pruned genotype files now have `1,107,868` and `1,055,527` SNPs, respectfully.

<br>

Now we need to combine the .fam file with the phenotype file we created earlier, so that the phenotypes are matched to their correct samples in the 6th column of the .fam file. Let's do this in R.

`add_phenotypes_to_fam_file_for_gemma.R`:
```r
# Add phenotype values to .fam files for use as genotype/phenotype files in GEMMA analysis
# Created: 2026-06-29
# Last updated: 2026-06-29
# Jason A. Toy


rm(list = ls())

setwd("/archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/")

library(tidyverse)



# Load in files

# PLINK .fam files
fam_allramet <- read_delim("genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.fam", col_names = FALSE, delim = " ")
fam_cpruned <- read_delim("genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.fam", col_names = FALSE, delim = " ")

str(fam_allramet)

# phenotype value files
ed50_allramet <- read_delim("phenotypes/Pacu_complete_ed50s.txt", col_names = FALSE, delim = "\t")
ed50_cpruned <- read_delim("phenotypes/Pacu_complete_cpruned_ed50s.txt", col_names = FALSE, delim = "\t")

str(ed50_allramet)

# Create new .fam files by replacing column 6 with matched phenotype values
pheno_allramet <- left_join(fam_allramet %>% select(-X6), ed50_allramet, by = c("X2" = "X1"))
pheno_cpruned <- left_join(fam_cpruned %>% select(-X6), ed50_cpruned, by = c("X2" = "X1"))

str(pheno_allramet)

# Write new .fam files
write_delim(pheno_allramet, "genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet_withpheno.fam", delim = " ", col_names = FALSE)
write_delim(pheno_cpruned, "genotypes/pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned_withpheno.fam", delim = " ", col_names = FALSE)
```

<br>

Now we need to rename the old and new fam files so that the new fam file names match the .bed and .bim file names:
```bash
cd /archive/barshis/barshislab/jtoy/pver_gwas/gemma_gwas/genotypes

## all-ramet ##

# rename old fam file
mv pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.fam pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet_nopheno.fam

# rename new fam file to match bed and bim
mv pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet_withpheno.fam pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.fam



## clone-pruned ##

# rename old fam file
mv pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.fam pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned_nopheno.fam

# rename new fam file to match bed and bim
mv pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned_withpheno.fam pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.fam
```

<br>

Double-check new .fam files:
```bash
head *[t,d].fam
```
```
==> pver_all_QDPSB_MISSMAF05filtered_genotypes_allramet.fam <==
0 2024_ALOF_Pver_02_1 0 0 0 37.274612836513
0 2024_ALOF_Pver_03_1 0 0 0 37.7280845350833
0 2024_ALOF_Pver_04_1 0 0 0 37.6325184444698
0 2024_ALOF_Pver_05_1 0 0 0 37.871523866322804
0 2024_ALOF_Pver_06_1 0 0 0 37.7970495591192
0 2024_ALOF_Pver_07_1 0 0 0 37.6121588898274
0 2024_ALOF_Pver_08_1 0 0 0 37.6306961285469
0 2024_ALOF_Pver_09_1 0 0 0 36.9392839584712
0 2024_ALOF_Pver_10_1 0 0 0 37.5853252077299
0 2024_ALOF_Pver_11_1 0 0 0 37.4656074378858

==> pver_all_QDPSB_MISSMAF05filtered_genotypes_cpruned.fam <==
0 2024_ALOF_Pver_03_1 0 0 0 37.7280845350833
0 2024_ALOF_Pver_06_1 0 0 0 37.7970495591192
0 2024_ALOF_Pver_07_1 0 0 0 37.6121588898274
0 2024_ALOF_Pver_08_1 0 0 0 37.6306961285469
0 2024_ALOF_Pver_09_1 0 0 0 36.9392839584712
0 2024_ALOF_Pver_10_1 0 0 0 37.5853252077299
0 2024_ALOF_Pver_19_1 0 0 0 37.5161130960978
0 2024_ALOF_Pver_23_1 0 0 0 37.9794565750908
0 2024_ALOF_Pver_24_1 0 0 0 37.8598758999653
0 2024_ALOF_Pver_27_1 0 0 0 38.046780127904
```
