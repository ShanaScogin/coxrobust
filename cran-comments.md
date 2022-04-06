## Test environments
* win-builder (R-release, R-devel, R-oldrelease)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
No functionality is changing with this release. Function `coxr` used in `modeLLtest` has been tested (correctly this time). Test was also added to automatically catch if the function downstream package depends on fails.

## Winbuilder check results
Note exists in winbuilder develop about change of maintainer.

## r-hub notes
There are two r-hub notes 1.) a file called 'coxrobust-manual.tex' that I cannot find and 2.) detritus in the temp directory called 'lastMiKTeXException' that I also cannot resolve. I am not sure if these are artifacts from former releases under previous maintainer, but please let me know if you have suggestions.
