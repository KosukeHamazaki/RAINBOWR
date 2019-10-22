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
