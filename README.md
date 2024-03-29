# RAINBOWR
###   Reliable Association INference By Optimizing Weights with R
#### Author : Kosuke Hamazaki (hamazaki@ut-biomet.org)
#### Date : 2019/03/25 (Last update: 2022/01/31)

## NOTE!!!!
### The paper for `RAINBOWR` has been published in PLOS Computational Biology (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007663). If you use this `RAINBOWR` in your paper, please cite `RAINBOWR` as follows:
#### Hamazaki, K. and Iwata, H. (2020) RAINBOW: Haplotype-based genome-wide association study using a novel SNP-set method. PLOS Computational Biology, 16(2): e1007663.


### The stable version for `RAINBOWR` package is now available at the [CRAN (Comprehensive R Archive Network)](https://cran.r-project.org/package=RAINBOWR).

### Please check the change in RAINBOWR with the version update from [NEWS.md](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/NEWS.md).

### The older version of `RAINBOWR` is `RAINBOW`, which is available at https://github.com/KosukeHamazaki/RAINBOW.
##### We changed the package name from `RAINBOW` to `RAINBOWR` because the original package name `RAINBOW` conflicted with the package `rainbow` (https://cran.r-project.org/package=rainbow) when we submitted our package to `CRAN` (https://cran.r-project.org/).

----------

In this repository, the `R` package `RAINBOWR` is available.
Here, we describe how to install and how to use `RAINBOWR`.

----------
## What is `RAINBOWR`
`RAINBOWR`(Reliable Association INference By Optimizing Weights with R) is a package to perform several types of `GWAS` as follows.

- Single-SNP GWAS by `RGWAS.normal` function
- SNP-set (, haplotype-block, or gene-set) GWAS by `RGWAS.multisnp` function (which tests multiple SNPs at the same time)
- Check epistatic (SNP-set x SNP-set interaction) effects by `RGWAS.epistasis` (very slow and less reliable)
- Single-SNP GWAS including tests for interaction with genetic background by `RGWAS.normal` function
- SNP-set (, haplotype-block, or gene-set) GWAS including tests for interaction with genetic background (or epistatic effects with polygenes) by `RGWAS.multisnp` function (which tests multiple SNPs at the same time)


`RAINBOWR` also offers some functions to solve the linear mixed effects model.

- Solve multi-kernel linear mixed effects model with `EM3.general` function (using `gaston`, `MM4LMM`, or `RAINBOWR` packages; fast for `gaston` and `MM4LMM`)
- Solve one-kernel linear mixed effects model with `EMM.cpp` function
- Solve multi-kernel linear mixed effects model with `EM3.cpp` function (for the general kernel, not so fast)
- Solve multi-kernel linear mixed effects model with `EM3.linker.cpp` function (for the linear kernel, fast)

By utilizing these functions, you can estimate the genomic heritability and perform genomic prediction (`GP`).

Finally, `RAINBOWR` offers other useful functions.

- `qq` and `manhattan` function to draw Q-Q plot and Manhattan plot
- `modify.data` function to match phenotype and marker genotype data
- `CalcThreshold` function to calculate thresholds for GWAS results
- `See` function to see a brief view of data (like `head` function, but more useful)
- `genetrait` function to generate pseudo phenotypic values from marker genotype
- `SS_GWAS` function to summarize GWAS results (only for simulation study)
- `estPhylo` and `estNetwork` functions to estimate phylogenetic tree or haplotype network and haplotype effects with non-linear kernels for haplotype blocks of interest.
- `convertBlockList` function to convert haplotype block list estimated by PLINK to the format which can be inputted as a `gene.set` argument in `RGWAS.multisnp`, `RGWAS.multisnp.interaction`, and `RGWAS.epistasis` functions.


## Installation
The stable version of `RAINBOWR` is now available at the [CRAN (Comprehensive R Archive Network)](https://cran.r-project.org/package=RAINBOWR). The latest version of `RAINBOWR` is also available at the `KosukeHamazaki/RAINBOWR` repository in the [`GitHub`](https://github.com/KosukeHamazaki/RAINBOWR), so please run the following code in the R console.

``` r
#### Stable version of RAINBOWR ####
install.packages("RAINBOWR")  


#### Latest version of RAINBOWR ####
### If you have not installed yet, ...
install.packages("devtools")  

### Install RAINBOWR from GitHub
devtools::install_github("KosukeHamazaki/RAINBOWR")
```

If you get some errors via installation, please check if the following packages are correctly installed.
(We removed a dependency on `rgl` package!)


``` r
Rcpp,      # also install `Rtools` for Windows user
plotly,    # Suggests
Matrix,
cluster,
MASS,
pbmcapply,
optimx,
methods,
ape,
stringr,
pegas,
ggplot2,     # Suggests
ggtree,      # Suggests, install from Bioconducter with `BiocManager::install("ggtree")`
scatterpie,  # Suggests
phylobase,   # Suggests
haplotypes,  # Suggests
rrBLUP,
expm,
here,        # Suggests
htmlwidgets,
Rfast,
adegenet,    # Suggests
gaston,
MM4LMM,
furrr,       # Suggests
future,      # Suggests
progressr,   # Suggests
foreach,     # Suggests, but stongly recommend the installation for Windows user to parallel computation
doParallel,  # Suggests
data.table   # Suggests
```

In `RAINBOWR`,  since part of the code is written in `Rcpp` (`C++` in `R`),  please check if you can use `C++` in `R`.
For `Windows` users,  you should install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

If you have some questions about installation, please contact us by e-mail (hamazaki@ut-biomet.org).


##  Usage
First, import `RAINBOWR` package and load example datasets. These example datasets consist of marker genotype (scored with {-1, 0, 1}, 1,536 SNP chip (Zhao et al., 2010; PLoS One 5(5): e10780)), map with physical position, and phenotypic data (Zhao et al., 2011; Nature Communications 2:467). Both datasets can be downloaded from `Rice Diversity` homepage (http://www.ricediversity.org/data/). Also, the dataset includes a list of haplotype blocks from the version 0.1.30. This list was estimated by the PLINK 1.9 (Taliun et al., 2014; BMC Bioinformatics, 15).

``` r
### Import RAINBOWR
require(RAINBOWR)

### Load example datasets
data("Rice_Zhao_etal")
Rice_geno_score <- Rice_Zhao_etal$genoScore
Rice_geno_map <- Rice_Zhao_etal$genoMap
Rice_pheno <- Rice_Zhao_etal$pheno
Rice_haplo_block <- Rice_Zhao_etal$haploBlock

### View each dataset
See(Rice_geno_score)
See(Rice_geno_map)
See(Rice_pheno)
See(Rice_haplo_block)
```
You can check the original data format by `See` function.
Then, select one trait (here, `Flowering.time.at.Arkansas`) for example.

``` r
### Select one trait for example
trait.name <- "Flowering.time.at.Arkansas"
y <- Rice_pheno[, trait.name, drop = FALSE]
```

For  GWAS, first you can remove  SNPs whose MAF <= 0.05 by `MAF.cut` function.

``` r
### Remove SNPs whose MAF <= 0.05
x.0 <- t(Rice_geno_score)
MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
x <- MAF.cut.res$x
map <- MAF.cut.res$map
```

Next, we estimate additive genomic relationship matrix (GRM) by using `calcGRM` function.

``` r
### Estimate genomic relationship matrix (GRM) 
K.A <- calcGRM(genoMat = x)
```

Next, we modify these data into the GWAS format of `RAINBOWR` by `modify.data` function.

``` r
### Modify data
modify.data.res <- modify.data(pheno.mat = y, geno.mat = x, map = map,
                               return.ZETA = TRUE, return.GWAS.format = TRUE)
pheno.GWAS <- modify.data.res$pheno.GWAS
geno.GWAS <- modify.data.res$geno.GWAS
ZETA <- modify.data.res$ZETA

### View each data for RAINBOWR
See(pheno.GWAS)
See(geno.GWAS)
str(ZETA)
```
`ZETA` is a list of genomic relationship matrix (GRM) and its design matrix.

Finally, we can perform `GWAS` using these data.
First, we perform single-SNP GWAS by `RGWAS.normal` function as follows.

``` r
### Perform single-SNP GWAS
normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS,
                           ZETA = ZETA, n.PC = 4, skip.check = TRUE, P3D = TRUE)
See(normal.res$D)  ### Column 4 contains -log10(p) values for markers
### Automatically draw Q-Q plot and Manhattan by default.
```

Next, we perform SNP-set GWAS by `RGWAS.multisnp` function.

``` r
### Perform SNP-set GWAS (by regarding 11 SNPs as one SNP-set)
SNP_set.res <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, 
                              n.PC = 4, test.method = "LR", kernel.method = "linear", 
                              gene.set = NULL, skip.check = TRUE, 
                              test.effect = "additive", window.size.half = 5, window.slide = 11)

See(SNP_set.res$D)  ### Column 4 contains -log10(p) values for markers
```

You can perform SNP-set GWAS with sliding window by setting `window.slide = 1`.
And you can also perform gene-set (or haplotype-block based) GWAS by assigning the following data set to `gene.set` argument. (You can check the example also by `See(Rice_haplo_block)` in R.)

ex.)

|  gene (or haplotype block)   |  marker | 
| :-----: | :------:| 
| haploblock_1    | id1005261 | 
| haploblock_1    | id1005263 | 
| haploblock_2    | id1009557 | 
| haploblock_2    | id1009616 | 
| haploblock_3    | id1020154 | 
| ...    | ... | 

For example, when using the list of haplotype blocks estimated by PLINK,
``` r
### Perform haplotype-block based GWAS (by using hapltype blocks estimated by PLINK)
haplo_block.res <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS, ZETA = ZETA, 
                              n.PC = 4, test.method = "LR", kernel.method = "linear", 
                              gene.set = Rice_haplo_block, skip.check = TRUE, 
                              test.effect = "additive")

See(haplo_block.res$D)  ### Column 4 contains -log10(p) values for markers
```

There is no significant block for this dataset because the number of markers and blocks is too small for this dataset. However, when whole-genome sequencing data is available, the impact of using SNP-set/gene-set/haplotype-block methods becomes larger and we strongly recommend you use these methods. Please see Hamazaki and Iwata (2020, PLOS Comp Biol) for more details of the features of these methods.


### Help
If you want some help before performing `GWAS` with `RAINBOWR`, please see the help for each function by `?function_name`.



## References
Kennedy, B.W., Quinton, M. and van Arendonk, J.A. (1992) Estimation of effects of single genes on quantitative traits. J Anim Sci. 70(7): 2000-2012.

Storey, J.D. and Tibshirani, R. (2003) Statistical significance for genomewide studies. Proc Natl Acad Sci. 100(16): 9440-9445.

Yu, J. et al. (2006) A unified mixed-model method for association mapping that accounts for multiple levels of relatedness. Nat Genet. 38(2): 203-208.

Kang, H.M. et al. (2008) Efficient Control of Population Structure in Model Organism Association Mapping. Genetics. 178(3): 1709-1723.

Kang, H.M. et al. (2010) Variance component model to account for sample structure in genome-wide association studies. Nat Genet. 42(4): 348-354.

Zhang, Z. et al. (2010) Mixed linear model approach adapted for genome-wide association studies. Nat Genet. 42(4): 355-360.

Endelman, J.B. (2011) Ridge Regression and Other Kernels for Genomic Selection with R Package rrBLUP. Plant Genome J. 4(3): 250.

Endelman, J.B. and Jannink, J.L. (2012) Shrinkage Estimation of the Realized Relationship Matrix. G3 Genes, Genomes, Genet. 2(11): 1405-1413.

Su, G. et al. (2012) Estimating Additive and Non-Additive Genetic Variances and Predicting Genetic Merits Using Genome-Wide Dense Single Nucleotide Polymorphism Markers. PLoS One. 7(9): 1-7.

Zhou, X. and Stephens, M. (2012) Genome-wide efficient mixed-model analysis for association studies. Nat Genet. 44(7): 821-824.

Listgarten, J. et al. (2013) A powerful and efficient set test for genetic markers that handles confounders. Bioinformatics. 29(12): 1526-1533.

Lippert, C. et al. (2014) Greater power and computational efficiency for kernel-based association testing of sets of genetic variants. Bioinformatics. 30(22): 3206-3214.

Jiang, Y. and Reif, J.C. (2015) Modeling epistasis in genomic selection. Genetics. 201(2): 759-768.

Hamazaki, K. and Iwata, H. (2020) RAINBOW: Haplotype-based genome-wide association study using a novel SNP-set method. PLOS Computational Biology, 16(2): e1007663.
