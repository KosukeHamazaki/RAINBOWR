# Oct 21, 2019, RAINBOWR version 0.1.8
## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel), R version 3.6.1


## R CMD check results
#### There were no ERRORs.

#### There was 1 WARNING:
##### * checking data for ASCII and uncompressed saves ... WARNING

Note: 
significantly better compression could be obtained

by using R CMD build --resave-data
old_size new_size compress
Rice_Zhao_etal.RData    207Kb    105Kb  xz



#### There were 4 NOTEs:
##### * checking CRAN incoming feasibility ... NOTE


##### * checking package dependencies ... NOTE
Package in Depends/Imports which should probably only be in LinkingTo:‘RcppEigen’


##### * checking examples ...
** running examples for arch 'i386' ... [13m] NOTE

** running examples for arch 'x64' ... [534s] NOTE
Examples with CPU or elapsed time > 10s




# Oct 22, 2019, RAINBOWR version 0.1.9 (resubmission)
## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel and release), R version 3.6.1


## R CMD check results
#### There were no ERRORs.

#### There were no WARNINGs:
- We resolved the [WARNING above](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#there-was-1-warning) by compressing the dataset using the `usethis::use_data` function and selecting `xz` compression.


#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE
- Still we have this NOTE message, but we think it is not a big problem. This is because this message is a remainder to check that the submission comes actually from his maintainer and not anybody else as the CRAN maintainer Dr. Uwe Ligges mentioned in https://mailman.stat.ethz.ch/pipermail/r-devel/2014-March/068497.html.

- We resolved [another NOTE](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#-checking-package-dependencies--note) by changing options of `Imports` in the `DECRIPTION` file and `import` in the `NAMESPACE` file for the package `RcppEigen`.

- We also resolved [the other NOTEs](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#-checking-examples-) by reducing the amounts of examples using the `\dontrun` function.





# Oct 29, 2019, RAINBOWR version 0.1.10 (resubmission)
## Major changes
We fixed some parts commented by the CRAN manager, Dr. Martina Schmirl.
We will describe the comments and their solution as follows.

> Please always explain all acronyms in the description text.

- We added the formal name of "SNPs", "GWAS", "RAINBOWR", and "LR" in the description text.


> If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form
> authors (year) <doi:...>
> authors (year) <arXiv:...>
> authors (year, ISBN:...)
> or if those are not available: <https:...>
> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking.

- We added the information about our preprint, "Kosuke Hamazaki and Hiroyoshi Iwata (2019) <doi:10.1101/612028>", in the DESCRIPTION file.


> You write information messages to the console that cannot be easily suppressed.
> Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console.
> (except for print() and summary() functions)

- We modified some print()/cat() functions as follows.
- print("ERROR: error messages.") --> stop("Error messages").
- cat("Warning!!: warning messages.") --> warning("Warning messages.")
- print("Print messages.") --> if (verbose) print("Print messages.")
- Now, the messages to the console can be easily suppressed. 


> \dontrun{} should be only used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user.
> Does not seem necessary.
> Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing.
> (You could also replace \dontrun{} with \donttest if it takes longer than 5 sec to be executed, but it would be preferable to have automatic checks for functions.)

- We replaced \dontrun in the examples to \donttest. We wanted to prepare small toy examples, but these small toy examples with their execution time < 5 sec deviated from the real situations in quantitative genetics, so we used \donttest instead.



## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel and release), R version 3.6.1


## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE





# Nov 2, 2019, RAINBOWR version 0.1.11 (resubmission)
## Major changes
We fixed some parts commented by the CRAN manager, Dr. Martina Schmirl.
We will describe the comments and their solution as follows.

> The package cannot be copyright holder of the package. (license file.)

- We changed the copyright holder to "Kosuke Hamazaki".


> Please ensure that your functions do not write by default or in your examples/vignettes/tests in the user's home filespace (including the package directory and getwd()). That is not allowed by CRAN policies. 
> Please only write/save files if the user has specified a directory in the function themselves. Therefore please omit any default path = getwd() in writing functions.
> In your examples/vignettes/tests you can write to tempdir().
> e.g. genetrait.R, ...

- We removed "getwd()" from the `genetrait` function.


> You still write information messages to the console that cannot be easily suppressed. 
> It is more R like to generate objects that can be used to extract the information a user is interested in, and then print() that object.
> Instead of print()/cat() rather use message()/warning()  or if(verbose)cat(..) if you really have to write text to the console.
> (except for print() and summary() functions)

- We modified some print()/cat() functions as follows.
- cat("\n") --> if (count) { cat("\n") }  ### in the functions related to `RGWAS`
- print()  --> if (verbose) { print() }  ### in the `See` function, we added `verbose` argument
- Now, the messages to the console can be easily suppressed.
You can check this from the examples in \dontshow{}.


> However, there is still a dontrun in EMM.cpp.Rd
> Also we would appreciate it to have examples for automatic test so you can wrap unrealistic examles in \dontshow. Then we can have the checks and the users do not see those examples.

- We replaced \dontrun in the example of `EMM.cpp` function to \donttest.
- We also prepared \dontshow examples (unrealistic toy examples) for the functions that had not offered the examples for the CRAN tests (`EM3.cpp`, `EM3.linker.cpp`, `RGWAS.epistasis`, `RGWAS.multisnp`, `RGWAS.twostep.epi`, `RGWAS.twostep`) before.



## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (release), R version 3.6.1


## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE




# Nov 6, 2019, RAINBOWR version 0.1.12 (resubmission)
## Major changes
We fixed some parts commented by the CRAN manager, Dr. Swetlana Herbrandt.
We will describe the comments and their solution as follows.

> please write package names, software names and API names in single quotes (e.g. 'RAINBOWR') in your Description text.

- We changed "By using RAINBOWR ..." in the Description text to "By using 'RAINBOWR' ...".


> You are changing the user's par() settings in your functions. Please ensure with an immediate call of on.exit() that the settings are reset.

- We modified `par(op)` in the `manhattan3` function (in the `functions_for_RGWAS.R` file) to `on.exit(par(op))` to reset the graphical settings.


## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel and release), R version 3.6.1


## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE




# Nov 8, 2019, RAINBOWR version 0.1.13 (resubmission)
## Major changes
We fixed some parts commented by the CRAN manager, Dr. Jelena Saf.
We will describe the comments and their solution as follows.

> Please shorten the title to a maximum of 65 characters. 

- We changed the title to "Genome-Wide Association Study with SNP-Set Methods" (51 characters).


> Please make sure that you do not change the user's options, par or working directory. If you really have to do so, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited, similar to this:

- We have to change the user's options `par`, so we modified the `manhattan3` function (in the `functions_for_RGWAS.R` file) as follows to reset the graphical settings. (We put the `on.exit(par(oldpar))` before the plotting section.)

```
oldpar <- par(no.readonly = TRUE)
on.exit(par(oldpar))
par(mar = c(3, 3, 3, 6), xpd = T)
plot(x, y, cex = pl.size, xlim = c(0, max(cum.pos)), ylim = c(0, max(cum.pos)),
col = col.id[col.num], pch = 1)
```


> The checktime of 664 seconds exceeds the limit of 10 minutes.
> Please add small executable examples in your Rd-files to illustrate the use of the exported function but also enable automatic testing.

- We had already prepared toy examples before, but we changed the toy examples to a smaller one. We also prepared the smaller version of vignettes because we thought that the building vignettes is the most time-consuming part.



## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel and release), R version 3.6.1


## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE





# Nov 14, 2019, RAINBOWR version 0.1.14 (resubmission)
## Major changes
We fixed some parts commented by the CRAN manager, Dr. Brian D. Ripley.
We will describe the comments and their solution as follows.

> See the CRAN check results page for your package, https://cran.r-project.org/web/checks/check_results_RAINBOWR.html .
> You were warned in 'Writing R Extensions' that calling math functions such as log pow sqrt with integer arguments is not portable.

- In the above URL, CRAN Package Check Results for `r-patched-solaris-x86` said:

> EMM_functions.cpp: In function ‘Rcpp::List spectralG_eigen(Rcpp::NumericMatrix, Rcpp::NumericMatrix, bool, bool)’: 
EMM_functions.cpp:1216:25: error: call of overloaded ‘sqrt(const int&)’ is ambiguous
double offset = sqrt(n);

- Then, we changed the C++ code in the `spectralG_eigen` function in the `EMM_functions.cpp` file as follows.

- Before
```
const int n(X.rows()), p(X.cols());   # L1207

double offset = sqrt(n);            # L1215
```

- After
```
const double n(X.rows()), p(X.cols());   # L1207

double offset = std::sqrt(n);            # L1215
```

- We also modified other parts related to the calling math functions (such as log pow sqrt) with integer arguments in the `EMM_functions.cpp` file.

- Actually, we are not sure that this change is correct because we cannot check the error disappears for the `Solaris` OS. So, if we have additional parts to be corrected, please tell us details with examples of correct codes.





## Test environments 
* platform x86_64-apple-darwin15.6.0, R version 3.6.0
* win-builder (devel and release), R version 3.6.1


## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE




# Nov 16, 2019, RAINBOWR version 0.1.15
## Major changes
We fixed some parts related to the treatment of the missing in marker genotypes.

- As the previous version, if the marker genotype has missing values, errors will be occurred in `RGWAS.normal`, `RGWAS.multisnp`, `RGWAS.epistasis`, `RGWAS.twostep`, and `RGWAS.twostep.epi` functions.



- Then, we changed the R code in the `RGWAS.normal`, `RGWAS.multisnp`, and `RGWAS.epistasis` functions as follows.

- Before
```
M.now <- Z.A[not.NA, ] %*% M
```

- After
```
if (sum(is.na(M)) == 0) {
M.now <- Z.A[not.NA, ] %*% M
} else {
M.now <- M[apply(Z.A[not.NA, ], 1, function(x) which(x == 1)), ]
}
```





# Apr 28, 2020, RAINBOWR version 0.1.16
## Major changes
- We added the citation file for RAINBOWR becuase our paper about RAINBOWR had been published in PLOS Computational Biology (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007663). We also changed the `Description` in the `DESCRIPTION` file to update the reference information.



- We added the `estPhylo` function to estimate and plot the phylogenetic tree for the block of interest. This function will also estimate the genotypic values for the block of interest. We also added `ape` and `cluster` packages in the `NAMESPACE` and `DESCRIPTION` files because these packages are used in the `estPhylo` function.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.0.0
* win-builder release, R version 4.0.0
* win-builder oldrelease, R version 3.6.3

## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE







# Apr 29, 2020, RAINBOWR version 0.1.17
## Major changes
- In version 0.1.16, the following NOTE was shown:

```
Found the following (possibly) invalid URLs:
URL: https://cran.r-project.org/web/packages/RAINBOWR/index.html
From: inst/doc/RAINBOWR.html
README.md
Status: 200
Message: OK
CRAN URL not in canonical form
The canonical URL of the CRAN page for a package is
https://CRAN.R-project.org/package=pkgname
```

Then, we fixed this by using the canonical URL of the CRAN package "https://cran.r-project.org/package=RAINBOWR" in `inst/doc/RAINBOWR.html` and `README.md` files.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.0.0
* win-builder release, R version 4.0.0
* win-builder oldrelease, R version 3.6.3

## R CMD check results
The same results as [those of the previous version 0.1.9](https://github.com/KosukeHamazaki/RAINBOWR/blob/master/cran-comments.md#r-cmd-check-results-1) were obtained as follows.

#### There were no ERRORs.

#### There were no WARNINGs:

#### There were 1 NOTE:
##### * checking CRAN incoming feasibility ... NOTE




# Mar 19, 2020, RAINBOWR version 0.1.18
## Major changes
- Important error about `spectralG.cpp` function was fixed thanks to Dr. Ishimori.
- `class(obj) == "try-error"` was modified to `try-error %in% class(obj)` in order to deal with the cases where the obj has more than one class (to avoid warnings).


# Jul 22, 2020, RAINBOWR version 0.1.19
## Major changes
- We added the `calcGRM` function to calculate genomic relationship matrix (GRM). The parts that calculate GRMs in other functions and examples were replaced to use `calcGRM` function. We also added `stringr` package in the `NAMESPACE` and `DESCRIPTION` files because these packages are used in the `calcGRM` function.

- We fixed some parts related to estimate the siginificance of dominance effects in `score.calc.LR`, `score.calc.LR.MC`, `score.calc.score`, `score.calc.score.MC`, `score.epistasis.LR`, `score.epistasis.score` functions.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.0.0
* win-builder release, R version 4.0.2
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 3.6.3

## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.




# Aug 01, 2020, RAINBOWR version 0.1.20
## Major changes
- We added the `estNetwork` function to estimate and plot haplotype network. We also added `pegas` package in the `NAMESPACE` and `DESCRIPTION` files because this package is used in the `estNetwork` function.

- We fixed some parts related to the estimation of haplotype effects in the `estPhylo` function. We also added the optimization part both in the `estPhylo` & `estNetwork` functions.




# Aug 18, 2020, RAINBOWR version 0.1.21
## Major changes
- We modified the `estPhylo` and `estNetwork` functions. Now, these functions can output the ggplot version of phylogenetic trees or haplotype networks. We also added `ggtree`, `ggplot2`, `scatterpie`, `phylobase`, `haplotypes`, and `ggimage` packages in the `NAMESPACE` and `DESCRIPTION` files because these packages are used in the `estPhylo` and `estNetwork` functions.

- We also added the `plotPhyloTree` and `plotHaploNetwork` functions to plot phylogenetic tree / haplotype network from the estimated results.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.0.2
* win-builder release, R version 4.0.3
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 3.6.3

## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.




# Oct 29, 2020, RAINBOWR version 0.1.22
## Major changes
- To deal with the error for the Windows user, we modified the package so that the user can select the number of cores even for solving the mixed-effects model with GEMMA (`EMM1.cpp` function). We also added the argument `n.core` for the functions that includes the mixed-effects model by `EMM1.cpp`.



# Nov 03, 2020, RAINBOWR version 0.1.23
## Major changes
- To deal with the error for the case where the chromosome numbers do not start from chromosome 1, we fix the scripts related to the treatment for genetic map. 


# Nov 11, 2020, RAINBOWR version 0.1.24
## Major changes
- We fixed an error on computing a Gaussian kernel for phylogenetic relationship in `estPhylo` function.
- We also fixed an error in `genesetmap` function related to the case where the markers in haplotype block list are not included in the marker genotype.


# Dec 15, 2020, RAINBOWR version 0.1.25
## Major changes
- We removed a dependency on `rgl` package since the future of the `rgl` package is unstable. Instead, we used `plot_ly` function in `plotly` package to draw 3d plots for results of epistasis GWAS. We added dependency on `plotly`, `here`, and `htmlwidgets` packages.


# Jan 19, 2021, RAINBOWR version 0.1.26
## Major changes
- We added a dependency on `Rfast` package to perform the faster computation of distance matrix. We also added an option of `pamonce = 5` in `cluster::pam` function to save the cost of performing k-medoids analysis.

- We removed a dependency on `pblapply` package since the package is actually not needed in the RAINBOWR package. We also removed `parallel` package from Imports list in DESCRIPTION file because `parallel` package is initially installed in R. We also removed `plotly` package from Imports list, and moved it to Suggests list in the DESCRIPTION file.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.0.2
* win-builder release, R version 4.0.3
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 3.6.3

## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.





# Apr 06, 2021, RAINBOWR version 0.1.27
## Major changes
- We removed `haplotypes` package from Imports list, and moved it to Suggests list in the DESCRIPTION file. This is because the package `phangorn`, which is used in the `haplotypes` package, is now scheduled for archival on 2021-04-19.



# Jan 05, 2022, RAINBOWR version 0.1.28
## Major changes for correcting errors pointed by CRAN maintainers
- We removed `ggtree` package from Imports list, and moved it to Suggests list in the DESCRIPTION file. By doing this, errors caused by the installation of `ggtree` package are now eliminated.


## Major changes (new functions, etc...)
- We implemented new functions for testing the interaction between each SNP and the genetic background. The functions to compute p-values for those effects are `score.calc.int` and `score.calc.int.MC`, and the function to perform SNP-based GWAS including such interaction effects is `RGWAS.normal.interaction`. We also added the to `NAMESPACE` file.

- We implemented a new function `is.diag`, which judges a matrix is diagonal or not.

- We implemented a new function `parallel.compute`, which enables us to perform parallel computation easily with the three different methods: `mclapply`, `furrr`, and `foreach`. This function is utilized in `score.calc.MC`, `score.calc.LR.MC`, `score.calc.score.MC`, `score.calc.epistasis.LR.MC`, `score.calc.epistasis.score.MC`, and `score.calc.int.MC` functions.

- We implemented a new function `EM3.general`, which enables us to solve mixed-effects model with the three different packages: `RAINBOWR`, `gaston`, and `MM4LMM`.

- We also added the four functions above to the NAMESPACE file to be exported correctly.

- We fixed some mistakes in `EM3.cpp` and `EM3.linker.cpp` functions when computing `Vinv`. We also added new arguments of `return.u.always`, `return.u.each`, and `return.Hinv`, and a new return of `u.each` as a `u` in the older version. We also modified a return of `u` to the summation of genotypic values.

- We rewrote the codes in `score.calc`, `score.calc.MC`, and `GWAS_F_test` (C++) functions so that it can test multiple fixed effects simultaneously.

- We fixed some parts in `RGWAS.normal`, `RGWAS.multisnp`, and `RGWAS.epstasis` functions when using `covariate`, `covariate.factor`, and `structure.matrix` arguments.

- We added a new argument `skip.check` in `RGWAS.normal`, `RGWAS.multisnp`, `RGWAS.epstasis`, `RGWAS.twostep`, and `RGWAS.twostep.epi` functions. As default, RAINBOWR checks the type of input data and modifies it into the correct format. However, it will take some time, so if you prepare the correct format of input data, you can skip this procedure by setting `skip.check = TRUE`.

- We introduced a new argument `package.MM` in `RGWAS.normal`, `RGWAS.multisnp`, `RGWAS.epstasis`, `RGWAS.twostep`, `RGWAS.twostep.epi`, `score.calc`, `score.calc.MC`, `score.calc.LR`, `score.calc.LR.MC`, `score.calc.epistasis.LR`, `score.calc.epistasis.LR.MC`, `score.calc.int`, and `score.calc.int.MC` functions. By changing this argument, you can choose which package is used to solve the mixed-effects model in GWAS from the following three packages: `RAINBOWR`, `gaston`, and `MM4LMM`.


- We introduced a new argument `parallel.method` in `RGWAS.normal`, `RGWAS.multisnp`, `RGWAS.epstasis`, `RGWAS.twostep`, `RGWAS.twostep.epi`, `score.calc.MC`, `score.calc.LR.MC`, `score.calc.score.MC`, `score.calc.epistasis.LR.MC`, `score.calc.epistasis.score.MC`, and `score.calc.int.MC` functions. In the older version, you can just use the `pbmcapply::pbmclapply` function for parallel computation, but now you can choose the parallel computation method from the following three methods: `mclapply`, `furrr`, and `foreach`.


- We largely corrected the functions related to epstatic tests; `RGWAS.epistasis`, `RGWAS.twostep.epi`, `score.calc.epistasis.LR`, `score.calc.epistasis.score`, `score.calc.epistasis.LR.MC`, and `score.calc.epistasis.score.MC` functions. This is because how to compute the interaction term between two haplotype blocks was wrong in the older version. We also added a new argument `skip.self.int` to choose whether skipping the computation for the epistatsis of self-interaction.


- We added `gaston` and `MM4LMM` packages to Imports list in the DESCRIPTION file because these packages are required when you use these packages with `package.MM` argument. We also added `lmm.aireml`, `lmm.diago`, and `MMEst` functions to be imported in the NAMESPACE file.


- We added `adegenet`, `furrr`, `future`, `progressr`, `foreach`, and `doParallel` packages to Suggests list in the DESCRIPTION file according to the implementation of new functions and arguments.We also removed `ggplot2`, `ggtree`, `scatterpie`, and `phylobase` packages from Imports list, and moved it to Suggests list in the DESCRIPTION file.




## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.1.2
* win-builder release, R version 4.1.2
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.0.5


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.





# Jan 06, 2022, RAINBOWR version 0.1.29
## Major changes for correcting errors pointed by CRAN maintainers
- We tried to fix the errors reported in the [RAINBOWR check results](https://cran.r-project.org/web/checks/check_results_RAINBOWR.html), and with our PC, now no error is found also in the test whose settings are same as those in the "gcc-ASAN" test.


## Major changes
- We modified the example files corresponding to `RGWAS.normal`, `RGWAS.normal.interaction`, `RGWAS.multisnp`, `RGWAS.epstasis`, `RGWAS.twostep`, and `RGWAS.twostep.epi` functions because there were some modifications with the function arguments.




## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.1.2
* win-builder release, R version 4.1.2
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.0.5


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.








# Feb 01, 2022, RAINBOWR version 0.1.30
## Major changes (new functions, etc...)
- We implemented new functions for testing the interaction between each SNP-set (haplotype block) and the genetic background (or epistasis with polygenes). The functions to compute p-values for those effects are `score.calc.LR.int` and `score.calc.LR.int.MC`, and the function to perform haplotype-block based GWAS including such interaction effects is `RGWAS.multisnp.interaction`. We also added the function to `NAMESPACE` file.

- We implemented a new function `adjustGRM`, which adjusts genomic relationship matrices when there is population structure. The function utilizes the true/estimated sub-population information (population membership) to estimated each variance component corresponding to each sub-population. We added the function to `NAMESPACE` file.

- We implemented a new function `convertBlockList`, which converts a list of haplotype blocks estimaed by PLINK to the format which can be inputted as a `gene.set` argument in `RGWAS.multisnp`, `RGWAS.multisnp.interaction`, and `RGWAS.epistasis` functions. We added the function to `NAMESPACE` file. We also added the `data.table` package to the Suggests list in the DESCRIPTION file.


- We removed the `RGWAS.normal` function from the imported functions in `NAMESPACE` file because the function did not catch up with the latest version of the `RAINBOWR` package.

- We added a new argument `map.gene.set` to `RGWAS.multisnp`, `RGWAS.multisnp.interaction`, `RGWAS.epstasis`, `RGWAS.twostep`, and `RGWAS.twostep.epi` functions. If this argument is NULL, the map will be constructed by `genesetmap` function after the SNP-set GWAS. It will take some time, so you can reduce the computational time by assigning this argument beforehand.

- We also updated the `README.md` and `RAINBOWR.md` for vignettes.




# Jun 25, 2022, RAINBOWR version 0.1.31
## Major changes (new functions, etc...)
- We added `subpop` argument to `calcGRM` function. By utilizing `subpop` argument, you can consider the difference of allele frequencies between sub-populations when computing the genomic relationship matrix. This argument is only valid when NOIA methods are selected.

- We added how to compute the marker effects from the GBLUP results in the example of `EMM.cpp` function.


# Jul 25, 2023, RAINBOWR version 0.1.32
## Major changes (new functions, etc...)
- We modified the error regarding `See` function when showing  `array` objects.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.3.1
* win-builder release, R version 4.3.1
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.2.3


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.






# Sep 11, 2023, RAINBOWR version 0.1.33
## Major changes (new functions, etc...)
- We added the new special sentinel "_PACKAGE" to the `RAINBOWR.R` file in order to automatically add a -package alias.
- We added `max.HE` argument to the `MAF.cut` function so that it can also remove markers with a large heterozygous rate.
- We added some codes to cut MAF for marker genotype beforehand in the `score.calc` and `score.calc.MC` functions.
- We added some codes to cut MAF for `gene.set` argument beforehand inside the `RGWAS.multisnp`, `RGWAS.epistasis`, and `RGWAS.multisnp.interaction` functions.
- We added some code to try performing other function to solve mixed effects model if the log-likelihood is returned as infinity when using the `score.calc.LR`, `score.calc.LR.MC`, and `EM3.general` functions.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.3.1
* win-builder release, R version 4.3.1
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.2.3


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.



# Oct 2, 2023, RAINBOWR version 0.1.34
## Major changes (new functions, etc...)
- We changed the argument name `nCores` to `n.core` and added the argument `parallel.method` in the `estPhylo` and `estNetwork` functions.
- We modified the `estPhylo` and `estNetwork` functions so that they can perform parallel computing even in Windows OS by using the `parallel.compute` function.
- We modified the `estPhylo` and `estNetwork` functions so that they can return a n x h matrix where n is the number of genotypes and h is the number of haplotypes for each block of interest. Also, they were extended so that they could perform k-medoids clustering to define the haplotype when the number of haplotypes is predefined.
- We also modified the `parallel.compute` function so that `parallel.method = 'foreach'` option can return the same values as `parallel.method = 'mclapply'`.


## Test environments 
* platform x86_64-apple-darwin17.0, R version 4.3.1
* win-builder release, R version 4.3.1
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.2.3


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.



# Mar 03, 2024, RAINBOWR version 0.1.35
## Major changes (new functions, etc...)
- We modified the `calcGRM` function so that it can compute the genomic relationship matrix with the input of marker genotype scored with {0, 1, 2}.
- We modified the `See` function so that it can show the shorter results for `list` objects.
- We added the `method = "Sidak"` option in the `CalcThreshold` function.
- We fixed the error "error: no match for 'operator/'" in `EMM_functions.cpp` by explicitly extracting each element from the matrix by using `.coeff`.



## Test environments 
* platform x86_64-apple-darwin20, R version 4.3.3
* win-builder release, R version 4.3.3
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.2.3


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.



# Mar 21, 2024, RAINBOWR version 0.1.36
## Major changes (new functions, etc...)
- We modified the `calcGRM` function so that it can compute the genomic relationship matrix with the input of marker genotype scored with {0, 1} or {-1, 0}.


# Dec 05, 2024, RAINBOWR version 0.1.37
## Major changes (new functions, etc...)
- We modified the `estPhylo`, `plotPhyloTree`, `estNetwork`, and `plotHaploNetwork` functions to fix some errors.

# May 21, 2025, RAINBOWR version 0.1.38
## Major changes (new functions, etc...)
- We enabled the `calcGRM` function to compute GRM by subsets of marker genotype. To achieve this, we added new arguments `batchSize` and `n.core` for the `calcGRM` function.
- We implemented a new function `EM3.cov`, which solves a multi-kernel linear mixed-effects model considering covariance structure between random effects. We added the function to `NAMESPACE` file.

## Test environments 
* platform aarch64-apple-darwin20, R version 4.5.0
* win-builder release, R version 4.5.0
* win-builder devel, R Under development (unstable)
* win-builder oldrelease, R version 4.4.3


## R CMD check results
#### Status: OK

#### There were no ERRORs.

#### There were no WARNINGs.

#### There were no NOTEs.


# June 26, 2025, RAINBOWR version 0.1.39
## Major changes (new functions, etc...)
- We fixed the error in `calcGRM` function with `methodGRM == "linear` option.
- We fixed the error in `estPhylo` and `plotPhyloTree` functions when drawing phylogenetic trees.
- We fixed the error in `estNetwork` and `plotHaploNetwork` functions when drawing haplotype network.
- We added new options for `EM3.cov` so that the user can determine initial weights for the genetic variance from outside.


# July 15, 2025, RAINBOWR version 0.1.40
## Major changes (new functions, etc...)
- We fixed the error in `EM3.cov` function when outputting `y.pred` for the case where Z is not square matrix.
- We fixed the error in `EM3.cov` function when solving the final uni-kernel mixed-effects model.
- We added `forceApproxK` argument in `EM3.cov` function to enable approximate a weighted kernel with a semi-positive definite matrix when the original kernel is not semi-positive definite.
