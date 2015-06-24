#!/usr/bin/env python
import dendropy
import sys
ERR = sys.stderr
OUT = sys.stdout
VERSION = '0.0.0'
NAME = 'slowdate'
DESCRIPTION = 'testing implementation of speed dating'
def error(msg):
  ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))

def fatal(msg):
  error(msg)
  sys.exit(1)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(NAME, description=DESCRIPTION)
  parser.add_argument('--version',
                      default=False,
                      action='store_true',
                      help='display program version')
  parser.add_argument('--show_tree',
                      default=False,
                      action='store_true',
                      help='write an ASCII version of the tree and time estimates to standard output')
  parser.add_argument('--method_sampled',
                      default=False,
                      action='store_true',
                      help='consider the sampling to be incomplete')
  parser.add_argument('--quite',
                      default=False,
                      action='store_true',
                      help='only emit warnings and fatal errors')
  parser.add_argument('--tree_file',
                    type=str,
                    required=False,
                    help='filepath to the input newick tree file.')
  parser.add_argument('--out_file',
                    type=str,
                    required=False,
                    help='filepath to the output for the annotated newick tree file.')
  parser.add_argument('--out_form',
                    type=str,
                    choices=('dated', 'ultrametric'),
                    default='ultrametric',
                    help='format of the outputree tree file. The default is ultrametric.')
  parser.add_argument('--grid',
                    type=int,
                    default=100,
                    help='number of discrete time bins in the approximation')
  parser.add_argument('--bd_lambda',
                    type=float,
                    default=1.0,
                    help='birth rate of the birth/death prior')
  parser.add_argument('--bd_mu',
                    type=float,
                    default=0.5,
                    help='death rate of the birth/death prior')
  parser.add_argument('--bd_rho',
                    type=float,
                    default=0.5,
                    help='the sampling rate (probability) for extant tips')
  parser.add_argument('--rate_mean',
                    type=float,
                    default=1.0,
                    help='the mean of the gamma distribution used as the prior on the rate of character evolution')
  parser.add_argument('--rate_variance',
                    type=float,
                    default=1.0,
                    help='the variance of the gamma distribution used as the prior on the rate of character evolution')
  args = parser.parse_args()
  if args.version:
    ERR.write('{n} v{v}\n'.format(n=NAME, v=VERSION))
    sys.exit(0)
  if args.tree_file is None:
    fatal('an input tree file is required.')