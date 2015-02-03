# FASTDATE

## Introduction

The aim of this project is to:

  * implement a fast divergence times estimation method.
  * extend previous work to provide interval estimates rather than just apporoximations of the MAP time estimate.
  * 64-bit design that handles very large datasets
  * multithreaded

## Compilation instructions

Currently, FASTDATE can be compiled using the makefile:

`make`

## FASTDATE command line options

General options:

* `--help`
* `--version`
* `--divtimes`

Input and output options:

* `--alignment-file`
* `--tree-file`
* `--out-file`

Model parameters:

* `--birth-rate`
* `--death-rate`
* `--edge-rate-mean`
* `--edge-rate-variance`

## Usage example

In the example below, we will estimate the divergence times of a dataset using FASTDATE.

`./fastdate --divtimes --tree-file tree.newick --alignment-file alignment.fasta --out-file output.log`

## FASTDATE license and third party licenses

The code is currently licensed under the GNU Affero General Public License version 3.

## Code

The code is currently written in C.

    File     | Description
-------------|------
**maps.c** | Character mapping arrays for converting sequences to internal representation.
**fastdate.c** | Main file handling command-line parametrs and executing corresponding parts.
**util.c** | Various common utility functions.
**Makefile** | Makefile

## Bugs

The source code has not been tested comprehensively yet. All bug reports are highly appreciated.

## The FASTDATE team

The following people have contributed to FASTDATE. In alphabetic order:

* Bastien Boussau
* Diego Darriba
* Tom&aacute;&scaron; Flouri
* Mark Holder
* Alexandros Stamatakis
