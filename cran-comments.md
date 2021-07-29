## Test environments
* local OS X install, R 4.1.0
* ubuntu 16.04.6 LTS (on travis-ci), R 4.1.0
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, Notes, or WARNINGs. 

* [Fixed] The Description field contains <https://doi.org/10.1214/18-EJS1434>). Please rather write <doi:10.1214/18-EJS1434>.

* [Fixed] Please write references in the description of the DESCRIPTION file in the form authors (year) <doi:...> authors (year) <arXiv:...> authors (year, ISBN:...) or if those are not available: authors (year) <https:...> with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add a title as well please put it in quotes: "Title")

* [Fixed] Please add \value to .Rd files regarding exported methods and explain the functions results in the documentation. 

* [Fixed] Some code lines in examples are commented out in HighDim.Rd

* [Fixed] checking CRAN incoming feasibility ... NOTE
Maintainer: 'Xiaorui (Jeremy) Zhu <zhuxiaorui1989@gmail.com>'

New submission

Possibly misspelled words in DESCRIPTION:
  EBIC (8:66)
  Liu (8:72)
  SPSP (5:21)

Found the following (possibly) invalid URLs:
  URL: https://cran.r-project.org/web/checks/check_results_SPSP.html
    From: README.md
    Status: 404
    Message: Not Found
    
* [Fixed] checking top-level files ... NOTE

Non-standard file/directory found at top level:
  'cran-comments.md'
  
## Downstream dependencies

All packages that I could install passed.
