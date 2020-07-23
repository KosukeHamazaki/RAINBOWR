## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dpi=96
)

## ---- include=TRUE------------------------------------------------------------
### Import RAINBOWR
require(RAINBOWR)

### Load example datasets
data("Rice_Zhao_etal")
Rice_geno_score <- Rice_Zhao_etal$genoScore
Rice_geno_map <- Rice_Zhao_etal$genoMap
Rice_pheno <- Rice_Zhao_etal$pheno

### View each dataset
See(Rice_geno_score)
See(Rice_geno_map)
See(Rice_pheno)

## ---- include=TRUE------------------------------------------------------------
### Select one trait for example
trait.name <- "Flowering.time.at.Arkansas"
y <- Rice_pheno[, trait.name, drop = FALSE]

## ---- include=TRUE------------------------------------------------------------
### Remove SNPs whose MAF <= 0.05
x.0 <- t(Rice_geno_score)
MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
x <- MAF.cut.res$x
map <- MAF.cut.res$map

## ---- include=TRUE------------------------------------------------------------
### Estimate genomic relationship matrix (GRM) 
K.A <- calcGRM(genoMat = x)

## ---- include=TRUE------------------------------------------------------------
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

## ---- include=TRUE------------------------------------------------------------
### Perform single-SNP GWAS
normal.res <- RGWAS.normal(pheno = pheno.GWAS, geno = geno.GWAS,
                           plot.qq = FALSE, plot.Manhattan = FALSE,
                           ZETA = ZETA, n.PC = 4, P3D = TRUE, count = FALSE)
See(normal.res$D)  ### Column 4 contains -log10(p) values for markers

## ---- echo=FALSE--------------------------------------------------------------
qq(normal.res$D[, 4])
manhattan(normal.res$D)
### Automatically draw Q-Q plot and Manhattan if you set plot.qq = TRUE and plot.Manhattan = TRUE.

## ---- include=TRUE, message=FALSE---------------------------------------------
### Perform SNP-set GWAS (by regarding 11 SNPs as one SNP-set, first 300 SNPs)
SNP_set.res <- RGWAS.multisnp(pheno = pheno.GWAS, geno = geno.GWAS[1:300, ], ZETA = ZETA,
                              plot.qq = FALSE, plot.Manhattan = FALSE, count = FALSE,
                              n.PC = 4, test.method = "LR", kernel.method = "linear", gene.set = NULL,
                              test.effect = "additive", window.size.half = 5, window.slide = 11)

See(SNP_set.res$D)  ### Column 4 contains -log10(p) values for markers

## ---- echo=FALSE--------------------------------------------------------------
qq(SNP_set.res$D[, 4])
manhattan(SNP_set.res$D)
### Automatically draw Q-Q plot and Manhattan if you set plot.qq = TRUE and plot.Manhattan = TRUE.

## ---- include=TRUE, eval=FALSE------------------------------------------------
#  RGWAS.menu()

