# FASTDATE

## Introduction

The aim of this project is to implement a fast divergence times estimation algorithm, based on the dynamic
programming algorithm by Oerjan Akerborg (2008). The new tool should:

* have an open source code with an appropriate open source license.
* 64-bit multi-threaded design that handles very large datasets
* be as accurate or more accurate than curent baysian methods.
* implement a variety of fossil calibration methods, such as tip dating, node dating and possibly, the fossilized birth-death process.

## Compilation instructions

Currently, FASTDATE can be compiled using the makefile:

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

In the example below, we will estimate the divergence times (`--divtimes`) of a dataset using FASTDATE. Prior to that
we print the input tree file in ascii art (`--show-treee`).

`./fastdate --divtimes --show-tree --tree-file tree.newick --out-file output.log`

## FASTDATE license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is currently written in C.

    File     | Description
-------------|------
**maps.c** | Character mapping arrays for converting sequences to internal representation.
**fastdate.c** | Main file handling command-line parametrs and executing corresponding parts.
**util.c** | Various common utility functions.
**tree.c** | Functions on the tree structure.
**Makefile** | Makefile
**newick.y** | Bison grammar file for newick binary tree parsing.
**newick.l** | Flex lexical analyzer.

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
