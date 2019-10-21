### Import RAINBOWR
require(RAINBOWR)

### Load example datasets
data("Rice_Zhao_etal")

### View each dataset
See(Rice_geno_score)
See(Rice_geno_map)
See(Rice_pheno)

### Select one trait for example
trait.name <- "Flowering.time.at.Arkansas"
y <- as.matrix(Rice_pheno[, trait.name, drop = FALSE])

### Remove SNPs whose MAF <= 0.05
x.0 <- t(Rice_geno_score)
MAF.cut.res <- MAF.cut(x.0 = x.0, map.0 = Rice_geno_map)
x <- MAF.cut.res$x
map <- MAF.cut.res$map


### Estimate genetic relationship matrix
K.A <- rrBLUP::A.mat(x) ### rrBLUP package can be installed by install.packages("rrBLUP")

### Modify data
modify.res <- modify.data(pheno.mat = y, geno.mat = x, return.ZETA = TRUE)
pheno.mat <- modify.res$pheno.modi
ZETA <- modify.res$ZETA


### Solve linear mixed effects model
EMM.res <- EMM.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
(Vu <- EMM.res$Vu)   ### estimated genetic variance
(Ve <- EMM.res$Ve)   ### estimated residual variance
(herit <- Vu / (Vu + Ve))   ### genomic heritability

(beta <- EMM.res$beta)   ### Here, this is an intercept.
u <- EMM.res$u   ### estimated genotypic values
See(u)
