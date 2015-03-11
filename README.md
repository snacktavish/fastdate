# FASTDATE

## Introduction

The aim of this project is to implement a fast divergence times estimation algorithm, based on the dynamic
programming algorithm by &Ouml;rjan  &Aring;kerborg (2008). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets
* be as accurate or more accurate than current bayesian methods.
* implement a variety of fossil calibration methods, such as tip dating, node dating and possibly, the fossilized birth-death process (Heath et al. 2014).

## Compilation instructions

Currently, FASTDATE can be compiled using the included makefile:

`make`

## FASTDATE command line options

General options:

* `--help`
* `--version`
* `--divtimes`
* `--show-tree`

Input and output options:

* `--tree-file`
* `--out-file`

Model parameters:

* `--birth-rate`
* `--death-rate`
* `--edge-rate-mean`
* `--edge-rate-variance`

## Usage example

In the example below, we estimate the divergence times (`--divtimes`) of a
dataset using FASTDATE. We use a birth-death process with birth rate
(`--birth-rate`) 2 and death rate (`--death-rate`) 1. The prior for the edge
rates is a gamma distribution with mean (`--edge-rate-mean) 1 and variance
(`edge-rate-variance`) 1. We output a diagram of the tree with time estimates
on screen (`--show-tree`) and write the annotated file in newick format in a
file (`--out-file`).

`./fastdate --divtimes --show-tree --grid-intervals 100 --birth-rate 2 --death-rate 1 --edge-rate-mean 1 --edge-rate-variance 1 --tree-file tree.newick --out-file output.newick`

## FASTDATE license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is currently written in C.

    File     | Description
-------------|------
**bd.c**  | Birth-death process related functions.
**fastdate.c** | Main file handling command-line parameters and executing corresponding parts.
**gamma.c** | Gamma density distribution functions.
**maps.c** | Character mapping arrays for converting sequences to internal representation.
**tree.c** | Functions on the tree structure.
**util.c** | Various common utility functions.
**Makefile** | Makefile
**newick.y** | Bison grammar file for newick binary tree parsing.
**lex.l** | Flex lexical analyzer.

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

* Heath T., Huelsenbeck JP, Stadler T. (2014)
**The fossilized birth–death process for coherent calibration of divergence-time estimates.**
*Proceedings of the National Academy of Sciences*, 111(29): 2957-2966.
doi:[10.1073/pnas.1319091111](http://dx.doi.org/10.1073/pnas.1319091111)
