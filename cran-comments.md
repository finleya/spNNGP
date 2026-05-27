## Test environments

* Local: R 4.5.2, x86_64-pc-linux-gnu, Ubuntu 25.10
* GitHub Actions: R release on ubuntu-latest, macos-latest, and windows-latest

## R CMD check results

0 errors | 0 warnings | 1 note

* Local `R CMD check --as-cran` reported:
  `checking for future file timestamps ... NOTE`
  `unable to verify current time`

  This appears to be an environment/time verification issue on the local
  checking machine.

## Reverse dependencies

There is one reverse dependency on CRAN, `cobin`.
