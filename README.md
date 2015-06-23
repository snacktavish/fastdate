# FASTDATE

## Introduction

The aim of this project is to implement a fast divergence times estimation
algorithm, based on the dynamic programming algorithm by &Ouml;rjan
&Aring;kerborg (2008). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.
* be as accurate or more accurate than current bayesian methods.
* implement a variety of fossil calibration methods, such as tip dating and node dating based on (Stadler 2010).

Currently, FASTDATE implements a faster version of the dynamic programming
approach running in quadratic time instead of the original cubic-time in
(&Aring;kerborg et al. 2008), and instead computes the optimal probability
density of a sampled tree with _n_ extant sampled leaves as shown in (Stadler
2009).


## Compilation instructions

Currently, FASTDATE can be compiled using the included makefile:

`make`

## FASTDATE command line options

General options:

* `--help`
* `--version`
* `--show_tree`
* `--method_sampled`
* `--quiet`
* `--threads`
* `--grid`

Input and output options:

* `--tree_file`
* `--out_file`
* `--out_form`

Model parameters:

* `--bd_lambda`
* `--bd_mu`
* `--bd_rho`
* `--rate_mean`
* `--rate_variance`

## Usage example

In the example below, we demonstrate the usage of FASTDATE for estimating
divergence times of a tree file in newick format (`--tree_file`) using the
sampling through time birth death process (`--method_sampled`) defined in
(Stadler 2009).
First, we discretize time into 1000 equally long intervals (`--grid`). We use a
birth-death process with birth rate (`--bd_lambda`) 2, death rate (`--bd_mu`)
1, and rate of sampled extant individuals 0.5 (`--bd_rho`). The prior for the
edge rates is a gamma distribution with mean (`--rate_mean`) 1 and variance
(`--rate_variance`) 1. We output a diagram of the tree with time estimates on
screen (`--show_tree`) and write newick format of the annotated tree in a file
(`--out_file`). Instead of producing a dated tree, we choose as output an
ultrametric tree (`--out_form`).

`./fastdate --method_sampled --show_tree --grid 1000 --bd_lambda 2 --bd_mu 1 --bd_rho 0.5 --rate_mean 2 --rate_variance 1 --tree_file tree.newick --out_file output.newick --out_form ultrametric`

## FASTDATE license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is written in C.

    File       | Description
---------------|----------------
**bd.c**       | Birth-death process related functions.
**fastdate.c** | Main file handling command-line parameters and executing corresponding parts.
**gamma.c**    | Gamma density distribution functions.
**maps.c**     | Character mapping arrays for converting sequences to internal representation.
**tree.c**     | Functions on the tree structure.
**util.c**     | Various common utility functions.
**arch.c**     | Architecture specific code (Mac/Linux).
**dp.c**       | Dynamic programming algorithm.
**Makefile**   | Makefile
**newick.y**   | Bison grammar file for newick binary tree parsing.
**lex.l**      | Flex lexical analyzer.

## Bugs

The source code has not been tested comprehensively yet. All bug reports are highly appreciated.

## The FASTDATE team

The following people have contributed to FASTDATE. In alphabetic order:

* Bastien Boussau
* Diego Darriba
* Tom&aacute;&scaron; Flouri
* Mark Holder
* Paschalia Kapli
* Emily Jane McTavish
* Alexandros Stamatakis

## References

* &Aring;kerborg &Ouml;., Sennblad B., & Lagergren J. (2008) 
**Birth-death prior on phylogeny and speed dating.**
*BMC evolutionary biology*, 8(1): 77.
doi:[10.1186/1471-2148-8-77](http://dx.doi.org/10.1186/1471-2148-8-77)

* Nee S., May RM, Harvey PH. (1994)
**The reconstructed evolutionary process.**
*Philosophical Transactions of the Royal Society of London B: Biological Sciences*, 344(1309): 305-311.
doi:[10.1098/rstb.1994.0068](http://dx.doi.org/10.1098/rstb.1994.0068)

* Gernhard T. (2008)
**The conditioned reconstructed process.**
*Journal of Theoretical Biology*, 253(4): 769-778.
doi:[10.1016/j.jtbi.2008.04.005](http://dx.doi.org/10.1016/j.jtbi.2008.04.005)

* Stadler T. (2009)
**On incomplete sampling under birth-death models and connections to the sampling-based coalescent.**
*Jounral of Theoretical Biology*, 261(1): 58-66.
doi:[10.1016/j.jtbi.2009.07.018](http://dx.doi.org/10.1016/j.jtbi.2009.07.018)

* Stadler T. (2010)
**Sampling-through-time in birth-death trees.**
*Journal of Theoretical Biology*, 267(3): 396-404.
doi:[doi:10.1016/j.jtbi.2010.09.010](http://dx.doi.org/10.1016/j.jtbi.2010.09.010)
