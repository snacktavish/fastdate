# FASTDATE

## Introduction

The aim of this project is to implement a fast divergence times estimation
algorithm, based on the dynamic programming algorithm by &Ouml;rjan
&Aring;kerborg (2008). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets.
* be as accurate or more accurate than current bayesian methods.
* implement a variety of fossil calibration methods, such as tip dating and node dating based on (Stadler 2010).
* auto-correlated rates (Thorne et al. 1998) and (Kishino et al. 2001).

Currently, FASTDATE implements a faster version of the dynamic programming
approach running in quadratic time instead of the original cubic-time in
(&Aring;kerborg et al. 2008), and instead computes the optimal probability
density of a sampled tree with _n_ extant sampled leaves as shown in (Stadler
2009). FASTDATE also implements a C version of the L-BFGS-B optimization method
(Byrd et al. 1995) converted from the Fortran implementations (Zhu et al.
1997). The method is used to optimize birth and death rate in case they are not
specified.


## Compilation instructions

Currently, FASTDATE can be compiled using the included makefile:

`make`

## FASTDATE command line options

General options:

* `--help`
* `--version`
* `--show_tree`
* `--method_relative`
* `--method_nodeprior`
* `--method_tipdates`
* `--quiet`
* `--threads`
* `--grid`

Input and output options:

* `--tree_file`
* `--prior_file`
* `--age_file`
* `--out_file`
* `--out_form`

Model parameters:

* `--bd_lambda`
* `--bd_mu`
* `--bd_rho`
* `--bd_psi`
* `--rate_mean`
* `--rate_variance`
* `--max_age`

## Usage example

In the example below, we demonstrate the usage of FASTDATE for estimating
(relative) divergence times of a tree file in newick format (`--tree_file`)
using the sampled birth death process (`--method_relative`) defined in (Stadler
2009).
First, we discretize time into 1000 equally long intervals (`--grid`). We use a
birth-death process with birth rate (`--bd_lambda`) 2, death rate (`--bd_mu`)
1, and rate of sampled extant individuals 0.5 (`--bd_rho`). The prior for the
edge rates is a gamma distribution with mean (`--rate_mean`) 1 and variance
(`--rate_variance`) 1. We output a diagram of the tree with time estimates on
screen (`--show_tree`) and write newick format of the annotated tree in a file
(`--out_file`). Instead of producing a dated tree, we choose as output an
ultrametric tree (`--out_form`).

`./fastdate --method_relative --show_tree --grid 1000 --bd_lambda 2 --bd_mu 1 --bd_rho 0.5 --rate_mean 2 --rate_variance 1 --tree_file example/small.tre --out_file output.newick --out_form ultrametric`


## Setting node priors

Node priors are set by specifying a file that contains calibration info with the `--prior_file` switch.
Note that priors can be set only when using the switches `--method_nodeprior` or `--method_tipdates`

Currently, the exponential, log-normal and uniform distributions are supported. 
Node priors can have the following format:

```
nodes exp (mean,offset)
nodes ln (mean,variance,offset)
nodes uni (min_age,max_age)
```

Note that `nodes` can be a single (tip) taxon when calibrating tip nodes (using
`--method_tipdates)`, or a comma-separated list of at least two taxa that
define the MRCA (most recent common ancestor) node on which the prior will be
placed. Two taxa are always sufficient to define an MRCA on a rooted tree,
however, more than two, i.e. an arbitrary subset of tips that define the
subtree rooted at the desired MRCA may be specified.

`mean` sets the mean of the distribution, `variance` sets the variance of the
log-normal distribution, and the x-axis is shifted by `offset` (or minimum age)
to the right, i.e.  x-axis starts at `offset`. `exp`, `ln` and `uni`
respectively indicate whether an exponential, log-normal or uniform
distribution is to be set.

Blank lines and any parts of the prior file prefixed with the symbol `#` are
considered as comments and are therefore skipped.




``./fastdate --method_nodeprior --tree_file example/small.tre --prior_file example/node_prior.txt --out_file node_prior.tre 
--bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 1 --show_tree --grid 100 --rate_mean 5 --rate_variance 1``

## Setting tip priors
Tip priors can be passed in the same way as node priors,

```
tip exp (mean,offset)
tip ln (mean,variance,offset)
```

and should be analyzed using the --method_tipdates and 
a value for the rate of sampling of tips in the past, --bd_psi, 
must be provided.

`./fastdate --method_tipdates --tree_file example/small.tre --prior_file example/tip_prior.txt --out_file tip_dates.tre 
--bd_mu 1 --bd_lambda 4 --max_age 100 --bd_rho 0.01 --show_tree --grid 100 --rate_mean 5 --rate_variance 1 --bd_psi 0.1`



## FASTDATE license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is written in C.

    File           | Description
-------------------|----------------
**bd.c**           | Birth-death process related functions.
**fastdate.c**     | Main file handling command-line parameters and executing corresponding parts.
**gamma.c**        | Gamma density distribution functions.
**exp.c**          | Exponensial distribution functions.
**ln.c**           | Log-normal distributions functions.
**lca.c**          | Lowest Common Ancestor computation.
**parse_prior.y**  | Bison grammar for parsing node priors.
**lex_prior.y**    | Flex lexical analyzer for node priors.
**nodeprior.c**    | Setting node priors to nodes.
**maps.c**         | Character mapping arrays for converting sequences to internal representation.
**tree.c**         | Functions on the tree structure.
**util.c**         | Various common utility functions.
**arch.c**         | Architecture specific code (Mac/Linux).
**dp.c**           | Dynamic programming algorithm.
**Makefile**       | Makefile.
**parse_newick.y** | Bison grammar file for newick binary tree parsing.
**lex_newick.l**   | Flex lexical analyzer for newick trees.

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
*Journal of Theoretical Biology*, 261(1): 58-66.
doi:[10.1016/j.jtbi.2009.07.018](http://dx.doi.org/10.1016/j.jtbi.2009.07.018)

* Stadler T. (2010)
**Sampling-through-time in birth-death trees.**
*Journal of Theoretical Biology*, 267(3): 396-404.
doi:[10.1016/j.jtbi.2010.09.010](http://dx.doi.org/10.1016/j.jtbi.2010.09.010)

* Thorne JL, Kishino H., Painter IS. (1998)
**Estimating the rate of evolution of the rate of molecular evolution.**
*Molecular Biology and Evolution*, 15(12): 1647-1657.
[PDF](http://mbe.oxfordjournals.org/content/15/12/1647.full.pdf)

* Kishino H., Thorne JL, Bruno WJ. (2001)
**Performance of a Divergence Time Estimation Method under a Probabilistic Model of Rate Evolution.**
*Molecular Biology and Evolution*, 18(3): 352-361.
[HTML](http://mbe.oxfordjournals.org/content/18/3/352.long)

* Byrd RH, Lu P., Nocedal J., Zhu C. (1995)
**A Limited Memory Algorithm for Bound Constrained Optimization.**
*SIAM Journal on Scientific Computing*, 16(5): 1190-1208.
doi:[10.1137/0916069](http://dx.doi.org/10.1137/0916069)

* Zhu C., Byrd RH, Lu P., Nocedal J. (1997)
**L-BFGS-B: Algorithm 778: L-BFGS-B, FORTRAN routines for large scale bound constrained optimization.**
*ACM Transactions on Mathematical Software*, 23(4): 550-560.
doi:[10.1145/279232.279236](http://dx.doi.org/10.1145/279232.279236)
