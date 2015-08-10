/*
    Copyright (C) 2015 Tomas Flouri, Emily Jane McTavish

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "fastdate.h"


static long inner_entries = 0;
static double interval_age = 0;
static int tree_height = 0;

/* this resets the node heights for method_nodeprior */
static void reset_node_heights(tree_node_t * node)
{
  assert(node);
  
  double offset = 0;
  int min_height = 0;

  /* number of lines required above the new min height which will
     be selected for the current node */
  int required_space = tree_height - node->height;

  /* compute minimum possible interval line given a node has calibration info */
  if (node->prior == NODEPRIOR_EXP)
    offset = ((exp_params_t *)(node->prior_params))->offset; 
  else if (node->prior  == NODEPRIOR_LN)
    offset = ((ln_params_t *)(node->prior_params))->offset;
  min_height = lrint(ceil(offset / interval_age));

  if (!node->left)
  {
    node->height = min_height;
    if (node->height + required_space >= opt_grid_intervals)
      fatal("Error: Calibration info on line %d of file %s disrupts grid", 
            node->prior_lineno, opt_priorfile);
    return;
  }

  reset_node_heights(node->left);
  reset_node_heights(node->right);

  node->height = (node->left->height > node->right->height) ? 
                       node->left->height + 1 : node->right->height + 1;

  if (min_height > node->height)
    node->height = min_height;
  
  if (node->height + required_space >= opt_grid_intervals)
    fatal("Error: Calibration info on line %d of file %s disrupts grid",
          node->prior_lineno, opt_priorfile);
}

static void alloc_node_entries(tree_node_t * node)
{
  int entries;

  assert(node);

  /* tip case */
  if (!node->left)
  {
    if (!node->prior)
    {
      /* tips are always placed on the first grid line */
      node->entries      = 1;
      node->matrix       = (double *)xmalloc(sizeof(double));
      node->matrix_left  = (int *)xmalloc(sizeof(int));
      node->matrix_right = (int *)xmalloc(sizeof(int));
    }
    else
    {
      entries = node->parent->entries + node->parent->height - node->height - 1;
      node->entries      = entries;
      node->matrix       = (double *)xmalloc((size_t)entries * sizeof(double));
      node->matrix_left  = (int *)xmalloc((size_t)entries * sizeof(int));
      node->matrix_right = (int *)xmalloc((size_t)entries * sizeof(int));
    }
    return;
  }

  if (!node->parent)
    entries = opt_grid_intervals - node->height;
  else
    entries = node->parent->entries + node->parent->height - node->height - 1;

  /* allocate storage space for placement information at each node */
  node->entries      = entries;
  node->matrix       = (double *)xmalloc((size_t)entries * sizeof(double));
  node->matrix_left  = (int *)xmalloc((size_t)entries * sizeof(int));
  node->matrix_right = (int *)xmalloc((size_t)entries * sizeof(int));
  node->matrix_PP       = (double *)xmalloc((size_t)entries * sizeof(double));

  /* for progress bar indication */
  inner_entries += entries;

  /* allocate the space for its two subtrees */
  alloc_node_entries(node->left);
  alloc_node_entries(node->right);

}

void dp_recurse(tree_node_t * node, int root_height)
{
  static long sum_entries = 0;

  int i,j,k;
  int jmax, kmax;
  unsigned int low;
  unsigned int left_low;
  unsigned int right_low;

  tree_node_t * left;
  tree_node_t * right;

  double rel_age_node, rel_age_left, rel_age_right, age_diff, abs_age_node;
  double prob_rate_left, prob_rate_right;
  double dist_logprob;
  double score, PPscore;

  if (!node) return;

  /* leaves case */
  if (!node->left)
  {
    if (!node->prior)
    {
      node->matrix[0] = 0.0;    /* 0 and not 1 because of the log-scale */
      return;
    }

    /* tip fossil case */
    for (i = 0; i < node->entries; ++i)
    {
      rel_age_node = (1.0 / opt_grid_intervals) * (node->height+i);
      abs_age_node = (node->height+i)*interval_age;
      /* if estimating absolute ages check for node priors and compute PDF */
      dist_logprob = 0;
      if (node->prior == NODEPRIOR_EXP)
      {
        exp_params_t * params = (exp_params_t *)(node->prior_params);
        dist_logprob = exp_dist_logpdf(1/params->mean, 
                                       abs_age_node - params->offset);
      }
      else if (node->prior == NODEPRIOR_LN)
      {
        ln_params_t * params = (ln_params_t *)(node->prior_params);
        dist_logprob = ln_dist_logpdf(params->mean,
                                      params->stdev,
                                      abs_age_node - params->offset);
      }
      else assert(0);

      node->matrix[i] = dist_logprob + bd_tipdates_prod_tip(rel_age_node);
    }

    return;
  }
  
  /* inner nodes case */
  left  = node->left;
  right = node->right;

  dp_recurse(left,  root_height);
  dp_recurse(right, root_height);

  low = node->height;
  left_low = left->height;
  right_low = right->height;

  /* run DP */

  /*
         
                        o  (rel_age_node)
                       / \
                      /   \
      (rel_age_left) o     \
                            o (rel_age_right)
                    
  */

  for (i = 0; i < node->entries; ++i)
  {
    rel_age_node = (1.0 / opt_grid_intervals) * (low+i);
    abs_age_node = (low+i)*interval_age;

    /* check the ages of the left child */
    if (!left->left)
      jmax = 1;
    else
      jmax = (node->height + i - left->height);
    
    assert(jmax <= left->entries);

    int jbest = -1;
    double jbest_score = -__DBL_MAX__;
    double jsum_score = 0;
    for (j = 0; j < jmax; ++j)
    {
      assert(j+left->height < node->height + i);
      rel_age_left = (1.0 / opt_grid_intervals) * (left_low+j);

      prob_rate_left = gamma_dist_logpdf(left->length / 
                                         (rel_age_node - rel_age_left));
      score = left->matrix[j] + prob_rate_left;
      jsum_score = -log(exp(-jsum_score)+exp(-score)); /*TODO need to think trhough if summing the probabilities is actually correct...*/
      if (score  > jbest_score)
      {
        jbest = j;
        jbest_score = score;
      }
    }

    /* check the ages of right child */
    if (!right->left)
      kmax = 1;
    else
      kmax = (node->height + i - right->height);

    assert(kmax <= right->entries);

    int kbest = -1;
    double kbest_score = -__DBL_MAX__;
    double ksum_score = 0;
    for (k = 0; k < kmax; ++k)
    {
      assert(k+right->height < node->height + i);
      rel_age_right = (1.0 / opt_grid_intervals) * (right_low+k);

      age_diff = rel_age_node - rel_age_right;
      prob_rate_right = gamma_dist_logpdf(right->length / age_diff);

      score = right->matrix[k] + prob_rate_right;
      ksum_score = -log(exp(-ksum_score)+exp(-score)); /*TODO need to think trhough if summing the probabilities is actually correct...*/
      if (score > kbest_score)
      {
        kbest = k;
        kbest_score = score;
      }
    }

    assert(jbest > -1);
    assert(kbest > -1);

    double bd_term = 0;
    if (opt_method_relative || opt_method_nodeprior)
      bd_term = bd_relative_prod(rel_age_node);
    else if (opt_method_tipdates)
      bd_term = bd_tipdates_prod_inner(rel_age_node);
    else assert(0);


    /* if estimating absolute ages check for node priors and compute PDF */
    dist_logprob = 0;
    if (node->prior == NODEPRIOR_EXP)
    {
      exp_params_t * params = (exp_params_t *)(node->prior_params);
      dist_logprob = exp_dist_logpdf(1/params->mean, 
                                     abs_age_node - params->offset);
    }
    else if (node->prior == NODEPRIOR_LN)
    {
      ln_params_t * params = (ln_params_t *)(node->prior_params);
      dist_logprob = ln_dist_logpdf(params->mean,
                                    params->stdev,
                                    abs_age_node - params->offset);
    }
    else if (node->prior) assert(0);

    score = bd_term + jbest_score + kbest_score + dist_logprob;
    PPscore = bd_term + jsum_score + ksum_score + dist_logprob;
    /* if it's the root add one more term */
    if (node->height == root_height)
    {
      if (opt_method_relative || opt_method_nodeprior)
        score += bd_relative_root(node->leaves,
                                  rel_age_node);
      else if (opt_method_tipdates)
        score += bd_tipdates_root(node->leaves,
                                   rel_age_node);
      else assert(0);
    }

    /* store best placement of children and likelihood for interval line i */
    node->matrix[i] = score;
    node->matrix_PP[i] = PPscore;
    node->matrix_left[i] = jbest;
    node->matrix_right[i] = kbest;
  }
  sum_entries += node->entries;
  progress_update(sum_entries);
}


static void dp_backtrack_recursive(tree_node_t * node, int best_entry)
{
  if (!node) return;
  node->interval_line = node->height + best_entry;

  if (node->left)
  {
    dp_backtrack_recursive(node->left, node->matrix_left[best_entry]);
    dp_backtrack_recursive(node->right, node->matrix_right[best_entry]);
  }
}

static void dp_backtrack(tree_node_t * root)
{
  int i;
  int best_entry = 0;
  double max_score = -__DBL_MAX__;
  double score;

  for (i = 0; i < root->entries; ++i)
  {
    score = root->matrix[i];

    if (score > max_score)
    {
      max_score = score;
      best_entry = i;
    }
  }
  root->interval_line = root->height + best_entry;

  dp_backtrack_recursive(root, best_entry);
}

void dp(tree_node_t * tree)
{
  assert(opt_grid_intervals > tree->height);
  tree_height = tree->height;

  gamma_dist_init();
  bd_init(0,0);

  /* compute the absolute age of an interval line */
  if (opt_method_nodeprior|| opt_method_tipdates)
    interval_age = opt_max_age / (opt_grid_intervals - 1);
  
  /* reset node heights according to calibrations if nodeprior method is used */
  if (opt_method_nodeprior || opt_method_tipdates)
    reset_node_heights(tree);

  /* allocate space for node entries */
  alloc_node_entries(tree);
  
  progress_init("Running DP...", inner_entries);
  dp_recurse(tree,tree->height);
  progress_done();

  if (!opt_quiet)
    printf ("Backtracking...\n");

  dp_backtrack(tree);

}
