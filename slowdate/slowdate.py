#!/usr/bin/env python
import math
import sys
try:
  import dendropy
except ImportError:
  sys.exit('dendropy 4 is required. Read the file that is helpfully called "README.md"\n')

ERR = sys.stderr
OUT = sys.stdout
VERSION = '0.0.0'
NAME = 'slowdate'
DESCRIPTION = 'testing implementation of speed dating'
VERBOSE = False
QUIET = False

def error(msg):
  ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))

def debug(msg):
  if VERBOSE:
    ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))
def info(msg):
  if not QUIET:
    ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))
def fatal(msg):
  error(msg)
  sys.exit(1)

def calc_stadler_c1(bd_lambda, bd_mu, bd_psi):
  f = bd_lambda - bd_mu - bd_psi
  s = 4.0 * bd_lambda * bd_psi
  return math.sqrt(f*f - s)

def calc_stadler_constants(bd_lambda, bd_mu, bd_rho, bd_psi):
  STADLER_C1 = calc_stadler_c1(bd_lambda, bd_mu, bd_psi)
  debug('STADLER_C1 = {s}'.format(s=STADLER_C1))
  c2_numerator = bd_lambda - bd_mu - 2 * bd_lambda * bd_rho - bd_psi
  STADLER_C2 = c2_numerator / STADLER_C1
  debug('STADLER_C2 = {s}'.format(s=STADLER_C2))
  return STADLER_C1, STADLER_C2

class StadlerFactors(object):
  def __init__(self,
               num_leaves_sampled,
               num_extinct_leaves_sampled,
               num_extinct_internals_sampled,
               bd_lambda,
               bd_mu,
               bd_rho,
               bd_psi):
    self.n = num_leaves_sampled
    self.m = num_extinct_leaves_sampled
    self.k = num_extinct_internals_sampled
    if self.m > 0 or self.k > 0:
      raise NotImplementedError('sampling of extinct species, is not supported yet')
    assert bd_rho > 0.0
    assert bd_rho <= 1.0
    self._c1, self._c2 = calc_stadler_constants(bd_lambda, bd_mu, bd_rho, bd_psi)
    self._ln_lambda = math.log(bd_lambda)
    self._ln_4rho = math.log(4.0*bd_rho)
    # store first term of denom of p1
    self._df = 2*(1.0 - self._c2 * self._c2) 
    # store coeff of second term of denom of p1
    s = 1 - self._c2
    self._ds = s*s  
    # store coeff of second term of denom of p1
    t = 1 + self._c2
    self._dt = t*t
    assert self.m == 0 and self.k == 0
    _root_ln_numerator = math.log(4*self.n * bd_rho * bd_lambda)
    _root_ln_denom_const = math.log(self._c1*(1.0 + self._c2))
    self._root_ln_const = _root_ln_numerator - _root_ln_denom_const
  def ln_root_factor(self, t):
    '''Calculates the log of the factor associated with the root being at time `t`
    assumes that the root has already been treated like an internal node.
    '''
    ln_var = math.log(1 - self._c2 + (1 + self._c2)*math.exp(self._c1*t))
    return self._root_ln_const - ln_var
  def ln_internal_node_factor(self, t):
    '''Calculates the log of the factor associated with an internal node begin at time t
    '''
    return self._ln_lambda + self._ln_p1(t)
  def _ln_p1(self, t):
    emag = self._c1 * t
    denominator = self._df + self._ds * math.exp(-emag) + self._dt * math.exp(emag)
    return self._ln_4rho - math.log(denominator)

def main(args):
  global VERBOSE, QUIET
  if args.quiet:
    QUIET = True
  if args.verbose:
    QUIET = False
    VERBOSE = True
  assert args.grid > 0, 'number of grid points must be positive'
  assert args.bd_lambda > 0.0, 'bd_lambda must be positive'
  assert args.bd_mu >= 0.0, 'bd_mu must be non-negative'
  assert args.bd_lambda > args.bd_mu, 'bd_lambda must be greater that bd_mu'
  assert args.bd_rho > 0.0, 'bd_rho must be positive'
  assert args.bd_rho <= 1.0, 'bd_rho cannot be greater than 1'
  assert args.rate_mean > 0.0, 'rate_mean must be positive'
  assert args.rate_variance >= 0.0, 'rate_variance must be non-negative'
  tree = dendropy.Tree.get(path=args.tree_file, schema='newick', rooting='force-rooted')
  for edge in tree.inorder_edge_iter():
    if (edge.tail_node is not None) and (edge.length is None):
      raise ValueError('Every branch in the input tree must have a branch length')
  # fill in min_height leaves at 0 (contemporaneous leaves assumption)
  for nd in tree.postorder_node_iter():
    if nd.is_leaf():
      nd.min_grid_idx = 0
    else:
      nd.min_grid_idx = 1 + max([c.min_grid_idx for c in nd.child_nodes()])
  if tree.seed_node.min_grid_idx >= args.grid:
    raise ValueError('Grid too small! A grid of at least {g} is required'.format(g=tree.seed_node.min_grid_idx + 1))
  # fill in min_height leaves at 0 (contemporaneous leaves assumption)
  for nd in tree.preorder_node_iter():
    if nd.is_leaf():
      nd.min_grid_idx = 0
    else:
      if nd is tree.seed_node:
        nd.max_grid_idx = args.grid - 1
      else:
        nd.max_grid_idx = nd.parent_node.max_grid_idx - 1
      assert nd.max_grid_idx >= nd.min_grid_idx
  num_leaves = len(tree.leaf_nodes())
  num_extinct_leaves_sampled = 0
  num_extinct_internals_sampled = 0
  bd_prior = StadlerFactors(num_leaves,
                            num_extinct_leaves_sampled,
                            num_extinct_internals_sampled,
                            args.bd_lambda,
                            args.bd_mu,
                            args.bd_rho,
                            args.bd_psi)

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
  parser.add_argument('--quiet',
                      default=False,
                      action='store_true',
                      help='only emit warnings and fatal errors')
  parser.add_argument('--verbose',
                      default=False,
                      action='store_true',
                      help='emit debug level messages')
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
                    default=1000,
                    help='number of discrete time bins in the approximation')
  parser.add_argument('--bd_lambda',
                    type=float,
                    default=2.0,
                    help='birth rate of the birth/death prior')
  parser.add_argument('--bd_mu',
                    type=float,
                    default=0.0,
                    help='death rate of the birth/death prior')
  parser.add_argument('--bd_rho',
                    type=float,
                    default=0.5,
                    help='the sampling rate (probability) for extant tips')
  parser.add_argument('--bd_psi',
                    type=float,
                    default=0.0,
                    help='the sampling rate (probability) for extinct tips')
  parser.add_argument('--rate_mean',
                    type=float,
                    default=5.0,
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
  try:
    main(args)
  except Exception as x:
    fatal('An error occurred when validating the constraints on the command line options.\n' \
          'The exception raised should give you some hints about what went wrong:\n' + str(x))
