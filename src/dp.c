/*
    Copyright (C) 2015 Tomas Flouri

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

static void alloc_node_entries(tree_node_t * node)
{
  int entries;

  if (!node) return;

  /* tip case */
  if (!node->left)
  {
    /* tips are always placed on the first grid line */
    node->entries      = 1;
    node->matrix       = (double *)xmalloc(sizeof(double));
    node->matrix_left  = (int *)xmalloc(sizeof(int));
    node->matrix_right = (int *)xmalloc(sizeof(int));
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

  double rel_age_parent, rel_age_left, rel_age_right, age_diff;
  double prob_rate_left, prob_rate_right;
  double score;

  if (!node) return;

  /* leaves case */
  if (!node->left)
  {
    /* this should be 0 not 1, because of the log-scale */
    node->matrix[0] = 0.0;
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
         
                        o  (rel_age_parent)
                       / \
                      /   \
      (rel_age_left) o     \
                            o (rel_age_right)
                    
  */

  for (i = 0; i < node->entries; ++i)
  {
    rel_age_parent = (1.0 / opt_grid_intervals) * (low+i);

    /* check the ages of the left child */
    if (!left->left)
      jmax = 1;
    else
      jmax = (node->height + i - left->height);
    
    assert(jmax <= left->entries);

    int jbest = -1;
    double jbest_sum = -__DBL_MAX__;
    for (j = 0; j < jmax; ++j)
    {
      assert(j+left->height < node->height + i);
      rel_age_left = (1.0 / opt_grid_intervals) * (left_low+j);

      prob_rate_left = gamma_dist_logpdf(left->length / 
                                         (rel_age_parent - rel_age_left));

      if (left->matrix[j] + prob_rate_left > jbest_sum)
      {
        jbest = j;
        jbest_sum = left->matrix[j] + prob_rate_left;
      }
    }

    /* check the ages of right child */
    if (!right->left)
      kmax = 1;
    else
      kmax = (node->height + i - right->height);

    assert(kmax <= right->entries);

    int kbest = -1;
    double kbest_sum = -__DBL_MAX__;
    for (k = 0; k < kmax; ++k)
    {
      assert(k+right->height < node->height + i);
      rel_age_right = (1.0 / opt_grid_intervals) * (right_low+k);

      age_diff = rel_age_parent - rel_age_right;
      prob_rate_right = gamma_dist_logpdf(right->length / age_diff);

      if (right->matrix[k] + prob_rate_right > kbest_sum)
      {
        kbest = k;
        kbest_sum = right->matrix[k] + prob_rate_right;
      }
    }

    assert(jbest > -1);
    assert(kbest > -1);

    double node_term = bd_nofossil_prod(rel_age_parent);
    score = node_term + jbest_sum + kbest_sum;

    /* if it's the root add one more term */
    if (node->height == root_height)
      score += bd_nofossil_root(node->leaves, 
                                rel_age_parent);

    /* store best placement of children and likelihood for interval line i */
    node->matrix[i] = score;
    node->matrix_left[i] = jbest;
    node->matrix_right[i] = kbest;
  }
  sum_entries += node->entries;
  progress_update(sum_entries);
}

static void dp_backtrack_recursive(tree_node_t * node, int best_entry)
{
  if (!node) return;
  if (!node->left || !node->right)
  {
    node->interval_line = 0;
  }
  else
  {
    node->interval_line = node->height + best_entry;
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

  gamma_dist_init();
  bd_init();

  /* allocate space for node entries */
  alloc_node_entries(tree);
  
  progress_init("Running DP...", inner_entries);
  dp_recurse(tree,tree->height);
  progress_done();

  if (!opt_quiet)
    printf ("Backtracking...\n");

  dp_backtrack(tree);
}
