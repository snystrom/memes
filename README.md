
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dremeR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

An R interface to the [MEME Suite](http://meme-suite.org/) family of
tools.

## Installation

### Development Version

You can install the development version of dremeR from
[GitHub](https://github.com/snystrom/dremeR) with:

``` r
# install.packages("remotes")
remotes::install_github("snystrom/dremeR")
```

## Detecting the MEME Suite

dremeR relies on a local install of the [MEME
Suite](http://meme-suite.org/). For installation instructions for the
MEME suite, see the [MEME Suite Installation
Guide](http://meme-suite.org/doc/install.html?man_type=web).

dremeR needs to know the location of the `meme/bin/` directory on your
local machine. You can tell dremeR the location of your MEME suite
install in 4 ways. dremeR will always prefer the more specific
definition if it is a valid path. Here they are ranked from most- to
least-specific:

1.  Manually passing the install path to the `meme_path` argument of all
    dremeR functions
2.  Setting the path using `options(meme_bin = "/path/to/meme/bin/")`
    inside your R script
3.  Setting `MEME_BIN=/path/to/meme/bin/` in your `.Renviron` file
4.  dremeR will try the default MEME install location `~/meme/bin/`

If dremeR fails to detect your install at the specified location, it
will fall back to the next option.

To verify dremeR can detect your MEME install, use
`check_meme_install()` which uses the search herirarchy above to find a
valid MEME install. It will report whether any tools are missing, and
print the path to MEME that it sees. This can be useful for
troubleshooting issues with your install.

``` r
library(dremeR)

# Verify that dremeR detects your meme install
# (returns all green checks if so)
# (I have MEME installed to the default location)
check_meme_install()
#> checking main install
#> ✔ /nas/longleaf/home/snystrom/meme/bin
#> checking util installs
#> ✔ /nas/longleaf/home/snystrom/meme/bin/dreme
#> ✔ /nas/longleaf/home/snystrom/meme/bin/ame
#> ✔ /nas/longleaf/home/snystrom/meme/bin/fimo
#> ✔ /nas/longleaf/home/snystrom/meme/bin/tomtom
```

``` r
# You can manually input a path to meme_path
# If no meme/bin is detected, will return a red X
check_meme_install(meme_path = 'bad/path')
#> checking main install
#> ✖ bad/path
```

## The Core Tools

  - `runDreme()`
  - `runTomTom()`
  - `runAme()`
  - `runFimo()`

## Visualizing Results

# FAQs

### How do I use dremeR/MEME on Windows?

The MEME Suite does not currently support Windows, although it can be
installed under [Cygwin](https://www.cygwin.com/) or the [Windows Linux
Subsytem](https://docs.microsoft.com/en-us/windows/wsl/install-win10)
(WSL). Please note that if MEME is installed on Cygwin or WSL, you must
also run R inside Cygwin or WSL to use dremeR.
