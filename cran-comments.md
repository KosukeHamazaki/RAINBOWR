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
