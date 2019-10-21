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


### Estimate additive genetic relationship matrix & epistatic relationship matrix
K.A <- rrBLUP::A.mat(x) ### rrBLUP package can be installed by install.packages("rrBLUP")
K.AA <- K.A * K.A   ### additive x additive epistatic effects


### Modify data
Z <- design.Z(pheno.labels = rownames(y),
              geno.names = rownames(K.A))  ### design matrix for random effects
pheno.mat <- y[rownames(Z), , drop = FALSE]
ZETA <- list(A = list(Z = Z, K = K.A),
             AA = list(Z = Z, K = K.AA))


### Solve multi-kernel linear mixed effects model (2 random efects)
EM3.res <- EM3.cpp(y = pheno.mat, X = NULL, ZETA = ZETA)
(Vu <- EM3.res$Vu)   ### estimated genetic variance
(Ve <- EM3.res$Ve)   ### estimated residual variance
(weights <- EM3.res$weights)   ### estimated proportion of two genetic variances
(herit <- Vu * weights / (Vu + Ve))   ### genomic heritability (additive, additive x additive)

(beta <- EM3.res$beta)   ### Here, this is an intercept.
u <- EM3.res$u   ### estimated genotypic values (additive, additive x additive)
See(u)
