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
