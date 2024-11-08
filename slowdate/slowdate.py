#!/usr/bin/env python
import math
import sys
import re
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
NEG_INF = float('-inf')

def error(msg):
  ERR.write('{n} ERROR: {m}\n'.format(n=NAME, m=msg))
def warn(msg):
  ERR.write('{n} WARNING: {m}\n'.format(n=NAME, m=msg))

def debug(msg):
  if VERBOSE:
    ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))
def info(msg):
  if not QUIET:
    ERR.write('{n}: {m}\n'.format(n=NAME, m=msg))
def fatal(msg):
  error(msg)
  sys.exit(1)

# Adding the following attr to DendroPy nodes:
#   min_grid_idx lowest legal bin index for this node
#   max_grid_idx highest legal bin index for this node
#   idx2DP an integer indexable collection guaranteed to be 
#     indexable in the range [min_grid_idx, max_grid_idx]
#     holding None if not filled in or (relative log density, traceback)
#     where `traceback` is the pair daughter indices that yielded the 
#     highest density for this node+height combination, or None in the case
#     of leaves.
#    map_idx 
# The tree gets an new attr:
#   root_map_est which is a tuple (relative_ln_posterior_density, root_bin_index)

def _max_relative_date_map_for_des(par_age, par_bin_idx, des_node, rate_prior_calc, grid):
  highest_ln_density_idx = None
  highest_ln_density = NEG_INF
  #debug('  par_bin_idx = {p} des_node.max_grid_idx = {d}'.format(p=par_bin_idx, d=des_node.max_grid_idx))
  if des_node.is_leaf():
    end_idx = 1
  else:
    assert par_bin_idx <= 1 + des_node.max_grid_idx
    end_idx = par_bin_idx
  #debug('des_node.min_grid_idx = {}, end_idx = {}'.format(des_node.min_grid_idx, end_idx))
  for i_d in xrange(des_node.min_grid_idx, end_idx):
    t_d = grid.to_abs_age(i_d)
    duration = par_age - t_d
    rate_of_mol_evol = des_node.edge.length / duration
    ln_rate_prior = rate_prior_calc(rate_of_mol_evol)
    assert (ln_rate_prior > NEG_INF)
    d_dp = des_node.idx2DP[i_d]
    assert (d_dp is not None)
    dp_ln_d = d_dp[0]
    if not (dp_ln_d > NEG_INF):
      fatal('-inf in des idx2DP table for index {}.\n{}'.format(i_d, des_node.idx2DP))
    ln_d = ln_rate_prior + dp_ln_d
    assert ln_d > NEG_INF
    if ln_d > highest_ln_density:
      highest_ln_density = ln_d
      highest_ln_density_idx = i_d
  assert highest_ln_density_idx is not None
  return highest_ln_density, highest_ln_density_idx

def ln_of_sum(prev_ln_sum, new_ln_term):
  '''Avoid underflow/overflow, but return:
      log(x + y)
  where x = exp(prev_ln_sum) and y = exp(new_ln_term)
  '''
  if (prev_ln_sum is None) or (new_ln_term - prev_ln_sum > 30):
    return new_ln_term # first term, or previous is lost in rounding error
  if prev_ln_sum - new_ln_term > 30:
    return prev_ln_sum # new term is lost in rounding error
  ln_m = min(prev_ln_sum, new_ln_term)
  prev_ln_sum -= ln_m
  new_ln_term -= ln_m
  # one is now 0 and the other is no greater than 30
  exp_sum = math.exp(prev_ln_sum) + math.exp(new_ln_term)
  return math.log(exp_sum) + ln_m


def _sum_date_table_for_des(par_age, par_bin_idx, des_node, rate_prior_calc, grid):
  if des_node.is_leaf():
    end_idx = 1
  else:
    assert par_bin_idx <= 1 + des_node.max_grid_idx
    end_idx = par_bin_idx
  des_ln_prob_entries = {}
  des_ln_sum = None
  for i_d in xrange(des_node.min_grid_idx, end_idx):
    t_d = grid.to_abs_age(i_d)
    duration = par_age - t_d
    rate_of_mol_evol = des_node.edge.length / duration
    ln_rate_prior = rate_prior_calc(rate_of_mol_evol)
    assert (ln_rate_prior > NEG_INF)
    d_dp = des_node.idx2DP[i_d]
    assert (d_dp is not None)
    dp_ln_d = d_dp[0]
    if not (dp_ln_d > NEG_INF):
      fatal('-inf in des idx2DP table for index {}.\n{}'.format(i_d, des_node.idx2DP))
    ln_d = ln_rate_prior + dp_ln_d
    des_ln_prob_entries[i_d] = ln_d
    des_ln_sum = ln_of_sum(des_ln_sum, ln_d)
  return des_ln_sum, des_ln_prob_entries


def calc_date_credible_interval(tree, age_prior_calc, rate_prior_calc, grid, interval_prob):
  '''Uses a dynamic programming approach to approximate a `interval_prob` marginal credible
  interval for the age of each non-leaf node.
  '''
  for u in tree.postorder_node_iter():
    debug('node u = {n}'.format(n=u._as_newick_string()))
    if u.is_leaf():
      assert not u.calibrations, 'tip calibrations and interval estimates are not jointly supported (yet)'
    else:
      u.leftTable = {}
      u.rightTable = {}
      u.leftNormConstants = {}
      u.rightNormConstants = {}
      v, w = u.child_nodes()
      for i_u in xrange(u.min_grid_idx, 1 + u.max_grid_idx):
        assert u.idx2DP[i_u] is None
        t_u = grid.to_abs_age(i_u)
        
        ln_sum_left, table_left = _sum_date_table_for_des(t_u, i_u, v, rate_prior_calc, grid)
        u.leftTable[i_u] = table_left
        u.leftNormConstants[i_u] = ln_sum_left
        
        ln_sum_right, table_right = _sum_date_table_for_des(t_u, i_u, w, rate_prior_calc, grid)
        u.rightTable[i_u] = table_right
        u.rightNormConstants[i_u] = ln_sum_right
        
        ln_age_factor = age_prior_calc.ln_internal_node_factor(t_u)
        for calib_dist in u.calibrations:
          ln_age_factor += calib_dist.ln_calibration_density(t_u)
        d_terms = [ln_sum_left, ln_sum_right, ln_age_factor]
        if u is tree.seed_node:
          d_terms.append(age_prior_calc.ln_root_factor(t_u))
        d_u = sum(d_terms)
        if not (d_u > NEG_INF):
          fatal('-inf ln density at {} time={}. terms={}'.format(i_u, t_u, d_terms))
        u.idx2DP[i_u] = (d_u, None) # store density and traceback info
  root = tree.seed_node
  ln_norm_constant = None
  for i_u in xrange(root.min_grid_idx, 1 + root.max_grid_idx):
    lnt = root.idx2DP[i_u][0]
    ln_norm_constant = ln_of_sum(ln_norm_constant, lnt)
  root.age_prob = {}
  root.ln_age_prob = {}
  for i_u in xrange(root.min_grid_idx, 1 + root.max_grid_idx):
    ln_unnorm_prob = root.idx2DP[i_u][0]
    ln_norm_prob = ln_unnorm_prob - ln_norm_constant
    age_prob = math.exp(ln_norm_prob)
    print 'root i_u =', i_u, ' Pr(age) = ', age_prob
    root.ln_age_prob[i_u] = ln_norm_prob
    root.age_prob[i_u] = age_prob

  for u in tree.preorder_node_iter():
    if u.is_leaf():
      continue
    create_child_age_probs(u)
    v, w = u.child_nodes()
    calc_cred_interval_of_age(u, grid, interval_prob)

def create_child_age_probs(node):
  assert not node.is_leaf()
  v, w = node.child_nodes()
  v.ln_age_prob = {}
  w.ln_age_prob = {}
  for i_u in xrange(node.min_grid_idx, 1 + node.max_grid_idx):
    ln_u_age_prob = node.ln_age_prob[i_u]
    if not v.is_leaf():
      add_to_age_prob(v.ln_age_prob, ln_u_age_prob, node.leftTable[i_u], node.leftNormConstants[i_u])
    if not w.is_leaf():
      add_to_age_prob(w.ln_age_prob, ln_u_age_prob, node.rightTable[i_u], node.rightNormConstants[i_u])
  if not v.is_leaf():
    fill_age_prob_from_ln_age_prob(v)
  if not w.is_leaf():
    fill_age_prob_from_ln_age_prob(w)

def fill_age_prob_from_ln_age_prob(node):
  p = 0.0
  node.age_prob = {}
  for i, ln_ap in node.ln_age_prob.items():
    ap = math.exp(ln_ap)
    #print 'age bin=', i, ' ap = ', ap, ' ln_ap =', ln_ap
    p += ap
    node.age_prob[i] = ap
  
def calc_cred_interval_of_age(node, grid, interval_prob):
  bins = node.age_prob.keys()
  min_i = min(bins)
  max_i = max(bins)
  high_tail_prob = 0.0
  tail_p = (1 - interval_prob)/2.0
  for i in range(max_i, min_i - 1, -1):
    hp = node.age_prob[i]
    high_tail_prob += hp
    if high_tail_prob > tail_p:
      cred_upper = 1 + i
      break
  low_tail_prob = 0.0
  for i in range(min_i, max_i + 1):
    hp = node.age_prob[i]
    low_tail_prob += hp
    if low_tail_prob > tail_p:
      cred_lower = i - 1
      break
  node.lower_age_limit_bin = cred_lower
  node.upper_age_limit_bin = cred_upper
  node.age_credible_interval_prob = interval_prob
  node.lower_age_limit = grid.to_abs_age(node.lower_age_limit_bin)
  node.upper_age_limit = grid.to_abs_age(node.upper_age_limit_bin)


def add_to_age_prob(ln_age_prob, ln_par_age_prob, ln_unnorm_prob, ln_norm_constant):
  for i, lup in ln_unnorm_prob.items():
    ln_norm_prob_conditional_on_par = lup - ln_norm_constant
    ln_norm_joint = ln_norm_prob_conditional_on_par + ln_par_age_prob
    ln_age_prob[i] = ln_of_sum(ln_age_prob.get(i), ln_norm_joint)

def calc_date_map_estimate(tree, age_prior_calc, rate_prior_calc, grid):
  '''Uses a dynamic programming approach to approximate a MAP estimate
  for the grid_idx for each non-leaf node.
  '''
  for u in tree.postorder_node_iter():
    debug('node u = {n}'.format(n=u._as_newick_string()))
    if u.is_leaf():
      if u.calibrations:
        for i_u in xrange(u.min_grid_idx, 1 + u.max_grid_idx):
          t_u = grid.to_abs_age(i_u)
          ln_age_factor = 0.0
          for calib_dist in u.calibrations:
            ln_age_factor += calib_dist.ln_calibration_density(t_u)
          assert ln_age_factor > NEG_INF
          u.idx2DP[i_u] = (ln_age_factor, None) # store density and traceback info
    else:
      v, w = u.child_nodes()
      for i_u in xrange(u.min_grid_idx, 1 + u.max_grid_idx):
        assert u.idx2DP[i_u] is None
        t_u = grid.to_abs_age(i_u)
        d_v, i_v = _max_relative_date_map_for_des(t_u, i_u, v, rate_prior_calc, grid)
        d_w, i_w = _max_relative_date_map_for_des(t_u, i_u, w, rate_prior_calc, grid)
        ln_age_factor = age_prior_calc.ln_internal_node_factor(t_u)
        for calib_dist in u.calibrations:
          ln_age_factor += calib_dist.ln_calibration_density(t_u)
        d_terms = [d_v, d_w, ln_age_factor]
        if u is tree.seed_node:
          d_terms.append(age_prior_calc.ln_root_factor(t_u))
        d_u = sum(d_terms)
        if not (d_u > NEG_INF):
          fatal('-inf ln density at {} time={}. terms={}'.format(i_u, t_u, d_terms))
        u.idx2DP[i_u] = (d_u, (i_v, i_w)) # store density and traceback info
  root = tree.seed_node
  highest_ln_density = NEG_INF
  highest_ln_density_idx = None
  for i_r in  xrange(root.min_grid_idx, 1 + root.max_grid_idx):
    ln_density = root.idx2DP[i_r][0]
    if ln_density > highest_ln_density:
      highest_ln_density = ln_density
      highest_ln_density_idx = i_r
  assert highest_ln_density_idx is not None
  tree.root_map_est = (highest_ln_density, highest_ln_density_idx)
  tree.seed_node.map_idx = tree.root_map_est[1]
  tree.seed_node.map_age = grid.to_abs_age(tree.seed_node.map_idx)
  for node in tree.preorder_internal_node_iter():
    dp_node = node.idx2DP[node.map_idx]
    traceback = dp_node[1]
    assert(dp_node is not None)
    v, w = node.child_nodes()
    v.map_idx = traceback[0]
    v.map_age = grid.to_abs_age(v.map_idx)
    w.map_idx = traceback[1]
    w.map_age = grid.to_abs_age(w.map_idx)

def calc_stadler_c1(bd_lambda, bd_mu, bd_psi):
  f = bd_lambda - bd_mu - bd_psi
  s = 4.0 * bd_lambda * bd_psi
  return math.sqrt(f*f + s)

def calc_stadler_constants(bd_lambda, bd_mu, bd_rho, bd_psi):
  STADLER_C1 = calc_stadler_c1(bd_lambda, bd_mu, bd_psi)
  debug('STADLER_C1 = {s}'.format(s=STADLER_C1))
  c2_numerator = bd_lambda - bd_mu - 2 * bd_lambda * bd_rho - bd_psi
  STADLER_C2 = - c2_numerator / STADLER_C1
  debug('STADLER_C2 = {s}'.format(s=STADLER_C2))
  return STADLER_C1, STADLER_C2

_E_200 = math.exp(200)
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
    debug('StadlerFactors n={} m={} k={}'.format(self.n, self.m, self.k))
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
    self._root_ln_denom_const_factor = self._c1*(1.0 + self._c2)
    self._root_ln_const = _root_ln_numerator
  def ln_root_factor(self, t):
    '''Calculates the log of the factor associated with the root being at time `t`
    assumes that the root has already been treated like an internal node.
    '''
    emag = self._c1*t
    f = 1 - self._c2
    offset = 0.0
    #debug('   emag = {e} {f} {s} {t}'.format(e=emag, f=self._df, s=self._ds, t=self._dt))
    if emag > 200:
      while emag > 200:
        f = f / _E_200
        emag -= 200
        offset += 200
    var_factor = f + (1 + self._c2)*math.exp(emag)
    ln_var =offset + math.log(var_factor*self._root_ln_denom_const_factor)
    return self._root_ln_const - ln_var
  def ln_internal_node_factor(self, t):
    '''Calculates the log of the factor associated with an internal node begin at time t
    '''
    return self._ln_lambda + self._ln_p1(t)
  def _ln_p1(self, t):
    emag = self._c1 * t
    #debug('   emag = {e} {f} {s} {t}'.format(e=emag, f=self._df, s=self._ds, t=self._dt))
    if emag > 200:
      f = self._df / _E_200
      emag -= 200
      offset = 200
      while emag > 200:
        f = f / _E_200
        emag -= 200
        offset += 200
      ln_denom = offset + math.log(f + self._dt*math.exp(emag))
    else:
      assert emag > 0
      denominator = self._df + self._ds * math.exp(-emag) + self._dt * math.exp(emag)
      ln_denom = math.log(denominator)
    return self._ln_4rho - ln_denom

class LnRelUncorrelatedGammaRatePrior(object):
  '''Functor for the log of the relative prob density at a time.'''
  def __init__(self, mean, variance):
    assert mean > 0.0
    assert variance > 0.0
    self.beta = mean / variance
    self.alpha = mean * self.beta
  def __call__(self, r):
    assert(r > 0.0)
    return (self.alpha - 1)*math.log(r) - (self.beta * r)

class Grid(object):
  def __init__(self, num_bins, max_age):
    self.max_age = max_age
    self.num_bins = num_bins
    self.max_idx = self.num_bins - 1
    self.mill_years_per_bin = float(max_age)/num_bins
  def to_abs_age(self, bin_index):
    assert bin_index < self.num_bins
    return bin_index * self.mill_years_per_bin
  def min_age_to_min_grid_idx(self, age, calib):
    assert age < self.max_age
    b = int(math.ceil(age/self.mill_years_per_bin))
    assert b < self.num_bins
    if (self.to_abs_age(b) == age) and calib.has_zero_density_at_zero:
      assert calib.ln_calibration_density(age + self.mill_years_per_bin) > NEG_INF
      return b + 1 # not all of our priors have
    return b

class OffsetDistribution(object):
  def __init__(self, mean, offset):
    self.mean = mean
    assert mean > 0.0
    self.offset = offset
    assert offset >= 0.0
  def ln_calibration_density(self, t):
    '''Return the log density for time t'''
    variate = t - self.offset
    if variate < 0.0:
      return NEG_INF
    return self.ln_density_from_offset(variate)

class OffsetLogNormalDistribution(OffsetDistribution):
  def __init__(self, mean, variance, offset):
    OffsetDistribution.__init__(self, mean , offset)
    self.variance = variance
    assert variance > 0.0
    self.constants = -0.5 * math.log(self.variance*2.0 * math.pi)
    self.coefficient = -1.0/(2*variance)
    self.has_zero_density_at_zero = True
  def ln_density_from_offset(self, v):
    try:
      lv = math.log(v)
      d = lv - self.mean
      ln_density =  -lv + self.constants + self.coefficient * d * d
      info(' ln v={v} ln_density = {d}'.format(v=v, d=ln_density))
      return ln_density
    except:
      warn('underflow at v = {}'.format(v))
      return NEG_INF

class OffsetExponentialDistribution(OffsetDistribution):
  def __init__(self, mean, offset):
    OffsetDistribution.__init__(self, mean , offset)
    self.hazard = 1.0/mean
    self.constant = math.log(self.hazard)
    self.has_zero_density_at_zero = False
  def ln_density_from_offset(self, v):
    return self.constant + self.hazard * v

_CALIBRATION_PAT = re.compile(r'^(.*) (.*)\s*[(](.*)[)]$')
_DIST_PAT = re.compile(r'^(\S.+)\s*[(](.*)[)]$')
def distribution_string_to_distribution(dist_str, dist_params, min_age):
  v = {}
  dist_name_to_type = {
    'ln': OffsetLogNormalDistribution,
    'exp': OffsetExponentialDistribution
  }
  try:
    dn = dist_str.lower().strip()
    dc = dist_name_to_type[dn]
    args = list(dist_params)
    v = [float(i) for i in args]
    if dn == 'ln':
      assert(len(v) == 2)
      d = dc(v[0], v[1], min_age)
    else:
      assert(len(v) == 1)
      d = dc(v[0], min_age)
  except Exception as x:
    error(str(x))
    raise RuntimeError('Could not parse the distribution string "{}"'.format(dist_str))
  return d

def parse_node_calibration_line(line):
  chunks = _CALIBRATION_PAT.match(line)
  if not chunks:
    raise RuntimeError('Expecting calibration format to be:\n<designator> [designator] {exp|ln} (mean,[variance],offset)')
  dist_str = chunks.group(2)
  param_str_list = [i for i in chunks.group(3).strip().split(',')]
  mrca_of = [i.strip() for i in chunks.group(1).strip().split(' ') if i.strip()]
  if len(mrca_of) == 1:
    mrca_of = mrca_of[0]
  try:
    min_age = float(param_str_list[-1])
  except:
    raise RuntimeError('Expecting minimum age to be a floating point number; found "{}" \n'.format(age_str))
  dist_params = param_str_list[:-1]
  distribution = distribution_string_to_distribution(dist_str, dist_params, min_age)
  
  return {'mrca_of': mrca_of,
          'min_age': min_age,
          'distribution': distribution}

def parse_node_calibrations(fn):
  if fn is None:
    return []
  calibrations = []
  with open(fn, 'rU') as inp:
    for line in inp:
      ls = line.strip()
      if ls:
        calibrations.append(parse_node_calibration_line(ls))
  return calibrations

def add_node_calibrations(tree, grid, node_calibrations):
  for nc in node_calibrations:
    mrca_of = nc['mrca_of']
    try:
      mrca = tree.mrca(taxon_labels=mrca_of)
    except Exception as x:
      error(str(x))
      raise RuntimeError('MRCA could not be found for "{}"'.format('", "'.join(mrca_of)))
    debug('Setting node calibration for mrca_of = {}'.format(str(mrca_of)))
    min_age = nc['min_age']
    prob_dist = nc['distribution']
    try:
      calibration_min_grid_idx = grid.min_age_to_min_grid_idx(min_age, prob_dist)
    except:
      raise RuntimeError('Calibration age of {m} exceeds grid max age of {a}'.format(m=min_age, a=args.max_age))
    mrca.min_grid_idx = max(mrca.min_grid_idx, calibration_min_grid_idx)
    try:
      mrca.calibrations.append(prob_dist)
    except:
      mrca.calibrations = [prob_dist]

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
  assert args.interval is None or (args.interval > 0.0 and args.interval < 1.0)
  
  tree = dendropy.Tree.get(path=args.tree_file, schema='newick', rooting='force-rooted')
  for edge in tree.inorder_edge_iter():
    if (edge.tail_node is not None) and (edge.length is None):
      raise ValueError('Every branch in the input tree must have a branch length')
  # let's do some node punching...
  for nd in tree.preorder_node_iter():
    nd.calibrations = []
    nd.min_grid_idx = 0

  # Add 1 because we don't count the 0 bin (for the tips) in the the UI
  grid = Grid(1 + args.grid, max_age=args.max_age)
  # read node calibrations and allow then to adjust min_age
  node_calibrations = parse_node_calibrations(args.prior_file)
  tip_calibrations = parse_node_calibrations(args.age_file)
  node_calibrations.extend(tip_calibrations)

  add_node_calibrations(tree, grid, node_calibrations)

  for nd in tree.postorder_node_iter():
    if not nd.is_leaf():
      child_based_idx = 1 + max([c.min_grid_idx for c in nd.child_nodes()])
      nd.min_grid_idx = max(nd.min_grid_idx, child_based_idx)

  if tree.seed_node.min_grid_idx >= grid.max_idx:
    raise ValueError('Grid too small! A grid of at least {g} is required'.format(g=tree.seed_node.min_grid_idx + 1))
  # fill in min_height leaves at 0 (contemporaneous leaves assumption)
  for nd in tree.preorder_node_iter():
    if nd.is_leaf() and (not nd.calibrations):
      nd.max_grid_idx = 0
      nd.idx2DP = [(1.0, None)]
    else:
      if nd is tree.seed_node:
        nd.max_grid_idx = grid.max_idx
      else:
        nd.max_grid_idx = nd.parent_node.max_grid_idx - 1
      if nd.max_grid_idx < nd.min_grid_idx:
        raise ValueError('Grid too small! (or max_age is too small to comfortably deal with all of the calibrations)')
      nd.idx2DP = [None] * (1 + nd.max_grid_idx)
    debug('node {}: min {} max = {}'.format(nd._as_newick_string(), nd.min_grid_idx, nd.max_grid_idx))
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
  rate_prior_calc = LnRelUncorrelatedGammaRatePrior(mean=args.rate_mean, variance=args.rate_variance)
  if args.interval is not None:
    calc_date_credible_interval(tree,
                                age_prior_calc=bd_prior,
                                rate_prior_calc=rate_prior_calc,
                                grid=grid,
                                interval_prob=args.interval)
    for nd in tree.preorder_node_iter():
      if not nd.is_leaf():
        nd.annotations.add_new(name='lower_age', value=nd.lower_age_limit)
        nd.annotations.add_new(name='upper_age', value=nd.upper_age_limit)
  else:
    calc_date_map_estimate(tree,
                           age_prior_calc=bd_prior,
                           rate_prior_calc=rate_prior_calc,
                           grid=grid)
    for nd in tree.preorder_node_iter():
      if not nd.is_leaf():
        nd.annotations.add_new(name='age', value=nd.map_age)

  tree.write(file=sys.stdout,
             schema='newick',
             suppress_annotations=False,
             annotations_as_nhx=True)



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
  parser.add_argument('--method_relative',
                      default=False,
                      action='store_true',
                      help='produce a MAP estimate of the relative divergence times using the sampled birth-death process (Stadler, 2010)')
  parser.add_argument('--method_nodeprior',
                      default=False,
                      action='store_true',
                      help='Expect a --prior_file=FN arg')
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
  parser.add_argument('--max_age',
                      type=float,
                      default=1.0,
                      help='the absolute age of the maximum bin of the dynamic programming table.')
  parser.add_argument('--prior_file',
                      type=str,
                      default=None,
                      help='Filepath to list of node calibrations.')
  parser.add_argument('--age_file',
                      type=str,
                      default=None,
                      help='Filepath to list of tip calibrations.')
  parser.add_argument('--interval',
                      type=float,
                      required=False,
                      default=None,
                      help='If supplied, it should be a sub-probability used to defined the probability of the credible interval')
  
  args = parser.parse_args()
  if args.version:
    ERR.write('{n} v{v}\n'.format(n=NAME, v=VERSION))
    sys.exit(0)
  if args.tree_file is None:
    fatal('an input tree file is required.')
  try:
    main(args)
  except Exception as x:
    raise
    fatal('An error occurred when validating the constraints on the command line options.\n' \
          'The exception raised should give you some hints about what went wrong:\n' + str(x))
