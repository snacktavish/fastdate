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

static double * matrix;
static int * matrix_left;
static int * matrix_right;
static int inner_node_entries;
static int matrix_entries;


static int set_matrix_chunks(tree_node_t * node, int start_offset)
{
  int offset;

  if (!node) return start_offset;
  if (!node->left || !node->right)
  {
    node->matrix = matrix + start_offset;
    node->matrix_left = matrix_left + start_offset;
    node->matrix_right = matrix_right + start_offset;
    node->entries = 1;

    return start_offset+1;
  }
  
  offset = set_matrix_chunks(node->right, 
                               set_matrix_chunks(node->left, start_offset));

  node->matrix       = matrix       + offset;
  node->matrix_left  = matrix_left  + offset;
  node->matrix_right = matrix_right + offset;
  node->entries = inner_node_entries;

  return (offset + inner_node_entries);
}

static void dp_init(tree_node_t * root)
{
  int inner_nodes = (root->leaves - 2);
  set_node_heights(root);

  inner_node_entries = (opt_grid_intervals - root->height + 1);
  matrix_entries = inner_node_entries*inner_nodes + root->leaves + 1;


  /* allocate matrix space */
  matrix = (double *)xmalloc((size_t)matrix_entries * sizeof(double));
  matrix_left = (int *)xmalloc((size_t)matrix_entries * sizeof(int));
  matrix_right = (int *)xmalloc((size_t)matrix_entries * sizeof(int));

  /* another tree traversal to set the matrix portions per node */
  set_matrix_chunks(root,0);
  root->entries = 1;
}

static void dp_kill(void)
{
  free(matrix);
  free(matrix_left);
  free(matrix_right);
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

  double prob_rate_left;
  double prob_rate_right;
  double prob_node_time;

  double itime, jtime, ktime;

  if (!node) return;

  /* leaves case */
  if (!(node->left) || !(node->right))
  {
    node->matrix[0] = 1.0;
    return;
  }

  /* inner nodes case */
  left  = node->left;
  right = node->right;

  dp_recurse(left,  root_height);
  dp_recurse(right, root_height);

  /* initialize discretization domain for the three nodes */
  if (node->height == root_height)
    low = opt_grid_intervals;
  else
    low = node->height - 1;

  left_low = left->height - 1;
  right_low = right->height - 1;

  
  /* run DP */
  for (i = 0; i < node->entries; ++i)
  {
    node->matrix[i] = -__DBL_MAX__;
    itime = (1.0 / opt_grid_intervals) * (low+i);
    
    jmax = (node->height == root_height) ? 
           left->entries : MIN(i+1,left->entries);
    
    int jbest = -1;
    double jbest_sum = -__DBL_MAX__;
    for (j = 0; j < jmax; ++j)
    {
      jtime = (1.0 / opt_grid_intervals) * (left_low+j);
      prob_rate_left = gamma_dist_logpdf(left->length / (itime - jtime));
      if (left->matrix[j] + prob_rate_left > jbest_sum)
      {
        jbest = j;
        jbest_sum = left->matrix[j] + prob_rate_left;
      }
    }

    kmax = (node->height == root_height) ? 
           right->entries : MIN(i+1,right->entries);
    int kbest = -1;
    double kbest_sum = -__DBL_MAX__;
    for (k = 0; k < kmax; ++k)
    {
      ktime = (1.0 / opt_grid_intervals) * (right_low+k);
      prob_rate_right = gamma_dist_logpdf(right->length / (itime - ktime));
      if (right->matrix[k] + prob_rate_right > kbest_sum)
      {
        kbest = k;
        kbest_sum = right->matrix[k] + prob_rate_right;
      }
    }

    assert(jbest > -1);
    assert(kbest > -1);

    prob_node_time = bd_prob(node->leaves, itime) + jbest_sum + kbest_sum;
    
    /* store best placement of children and likelihood for interval line i */
    if (prob_node_time > node->matrix[i])
    {
      node->matrix[i] = prob_node_time;
      node->matrix_left[i] = jbest;
      node->matrix_right[i] = kbest;
    }
  }
  sum_entries += node->entries;
  progress_update(sum_entries);
}

static void dp_backtrack(tree_node_t * node, int best_entry)
{
  if (!node) return;
  if (!node->left || !node->right)
  {
    node->interval_line = 0;
  }
  else
  {
    node->interval_line = node->height - 1 + best_entry;
    dp_backtrack(node->left, node->matrix_left[best_entry]);
    dp_backtrack(node->right, node->matrix_right[best_entry]);
  }
}

void dp(tree_node_t * tree)
{
  assert(opt_grid_intervals > tree->height);

  gamma_dist_init(opt_edgerate_mean, opt_edgerate_var);
  bd_init(opt_birth_rate,opt_death_rate);

  /* init */
  dp_init(tree);
  
  progress_init("Running DP...", matrix_entries);
  dp_recurse(tree,tree->height);
  progress_done();

  if (!opt_quiet)
    printf ("Backtracking...\n");
  dp_backtrack(tree,0);
  tree->interval_line = opt_grid_intervals;

  dp_kill();
}
