#!/bin/bash
../src/fastdate --method_nodeprior --tree_file ../example/small.tre --prior_file ../example/node_prior.txt --out_file test2 --bd_mu 2 --bd_lambda 3 --max_age 100 --bd_rho 0.5 --show_tree --grid 100
