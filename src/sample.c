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

static double * cdf_vector;
static double interval_age;

static long sample_gridline(double * vector, long count) 
{
  long i;

  /* draw a random number between [0,1] */
  double r =  drand48();

  /* find the cdf/bin/gridline in which the random number falls in */
  for (i = 0; i < count; i++)
    if (r < vector[i])
      return i;

  assert(0);
}

static void recompute_scores(tree_node_t * node,
                             double rel_age_parent,
                             double * vector,
                             long count,
                             double * maxscore)
{
  long i;
  long node_low = node->height;
  double score;
  double rel_age_node;
  double prob_rate_node;

  *maxscore = -__DBL_MAX__;

  for (i = 0; i < count; ++i)
  {
    rel_age_node = (1.0 / opt_grid_intervals) * (node_low+i);

    /* TODO: Point of conflict */
    prob_rate_node = gamma_dist_logpdf(node->length /
                                       (rel_age_parent - rel_age_node));

    score = node->matrix[i] + prob_rate_node;
    if (score > *maxscore)
      *maxscore = score;
    vector[i] = score;
  }
}

static void normalize_cdf(double * vector,
                          long count,
                          double maxscore)
{
  long i;

  /*Subtract max from all*/
  for (i = 0; i < count; ++i)
    vector[i] = vector[i] - maxscore;

  /* Compute threshold for precision */
  double thresh = log(10e-16) - log(count);
  double sumscore = 0;

  /* discarding all scores below threshold and convert the rest to their
     exponentiations. Store the sum of exponentiations */
  for (i = 0; i < count; ++i)
  {
    if (vector[i] < thresh)
      vector[i] = 0;
    else
    {
      vector[i] = exp(vector[i]);
      sumscore = vector[i] + sumscore;
    }
  }
  assert (sumscore > 0);

  /* Normalize and store cumulative probability vector */
  vector[0] /= sumscore;
  for (i = 1; i < count; ++i)
  {
    vector[i] = (vector[i] / sumscore) + vector[i-1];
  }

  assert(0.9999 < vector[count-1]);
  assert(vector[count-1] < 1.0001);
}

static void dp_backtrack_sampling_recursive(tree_node_t * node)
{
  double maxscore;

  if (!node) return;
  assert(node->parent);

  tree_node_t * parent = node->parent;

  double rel_age_parent = (1.0 / opt_grid_intervals) *
                          (parent->sampled_gridline);

  /* number of grid-line entries on which we can place node */
  long entries = parent->sampled_gridline - node->height;

  /* in case we have a tip node without a prior, then it can have at most
     one entry */
  if (entries > node->entries)
    entries = node->entries;

  recompute_scores(node, rel_age_parent, cdf_vector, entries, &maxscore);
  normalize_cdf(cdf_vector, entries, maxscore);
  node->sampled_gridline = node->height + 
                           sample_gridline(cdf_vector, entries);

  dp_backtrack_sampling_recursive(node->left);
  dp_backtrack_sampling_recursive(node->right);
}

static void dp_backtrack_sampling(tree_node_t * root)
{
  int i;
  double maxscore = -__DBL_MAX__;

  cdf_vector = (double *)xmalloc((size_t)opt_grid_intervals * sizeof(double));

  /* copy logscores to cdf vector and get max */
  memcpy(cdf_vector, root->matrix, (size_t)(root->entries) * sizeof(double));
  for (i = 0; i < root->entries; ++i)
    if (cdf_vector[i] > maxscore)
      maxscore = cdf_vector[i];

  /* normalize logscores according to maxscore */
  normalize_cdf(cdf_vector, root->entries, maxscore);

  /* select interval line for root */
  root->sampled_gridline = root->height + 
                           sample_gridline(cdf_vector, root->entries);

  dp_backtrack_sampling_recursive(root->left);
  dp_backtrack_sampling_recursive(root->right);

  free(cdf_vector);
}

static void output_sample_dated_tree_recursive(tree_node_t * node,
                                               FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s[&age=%f]:%f",
            node->label, node->sampled_gridline * interval_age, node->length);
  else
  {
    fprintf(fp_out, "(");
    output_sample_dated_tree_recursive(node->left, fp_out);
    fprintf(fp_out, ",");
    output_sample_dated_tree_recursive(node->right, fp_out);
    fprintf(fp_out, ")%s[&age=%f]:%f", node->label ? node->label : "",
                    node->sampled_gridline * interval_age, node->length);
  }
}

static void output_sample_um_tree_recursive(tree_node_t * node, 
                                             FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s:%f", node->label, 
                    (node->parent->sampled_gridline - node->sampled_gridline) *
                      interval_age);
  else
  {
    fprintf(fp_out, "(");
    output_sample_um_tree_recursive(node->left, fp_out);
    fprintf(fp_out, ",");
    output_sample_um_tree_recursive(node->right, fp_out);
    if (node->parent)
     {
        fprintf(fp_out, ")%s:%f", node->label ? node->label : "", 
                      (node->parent->sampled_gridline - node->sampled_gridline)*
                         interval_age);
      }
   else
     {
        fprintf(fp_out, ")%s", node->label ? node->label : "");
      }
  }
}
void sample(tree_node_t * root)
{
  long i;

  interval_age = opt_max_age / (opt_grid_intervals - 1);

  char * filename = (char *)xmalloc((strlen(opt_outfile)+9)*sizeof(char));

  strcpy(filename, opt_outfile);
  strcat(filename,".sampled");

  FILE * fp_out = fopen(filename, "w");

  progress_init("Sampling...", (unsigned long)opt_sample);
  for (i = 0; i < opt_sample; ++i)
  {
    dp_backtrack_sampling(root);
    if (opt_outform == OUTPUT_ULTRAMETRIC)
      output_sample_um_tree_recursive(root, fp_out);
    else if (opt_outform == OUTPUT_DATED)
      output_sample_dated_tree_recursive(root, fp_out);
    else
      fatal("Internal error while selecting output format");
    fprintf(fp_out, ";\n");
    progress_update((unsigned long)i);
  }
  progress_done();

  fclose(fp_out);
  free(filename);
}

static double logsum(double old_logsum, double logterm)
{
   /*Avoid underflow/overflow, but return:
         log(x + y)
   where x = exp(prev_ln_sum) and y = exp(new_ln_term)
   */

   double ln_m;
   double exp_sum;

   /* check if old_logsum is lost in rounding error */
   if (logterm - old_logsum > 30)
     return logterm;

   /* check if logterm is lost in rounding error */
   if (old_logsum - logterm > 30)
     return old_logsum;

   ln_m = (old_logsum < logterm) ? old_logsum : logterm;

   old_logsum -= ln_m;
   logterm -= ln_m;

   exp_sum = exp(old_logsum) + exp(logterm);

   return (log(exp_sum) + ln_m);
}

static void cr_update(tree_node_t * node, double * logprob)
{
  long i;
  double tailsum = 0;
  long cred_upper = 0;
  long cred_lower = 0;

  double thresh = (1 - opt_cred_interval)/2;

  /* start from the high side */
  for (i = node->entries - 1; i >= 0; --i)
  {
    tailsum += exp(logprob[i]);
    if (tailsum > thresh)
    {
      cred_upper = node->height + i;
      break;
    }
  }

  /* start from the low side */
  tailsum = 0;
  for (i = 0; i < node->entries; ++i)
  {
    tailsum += exp(logprob[i]);
    if (tailsum > thresh)
    {
      cred_lower = node->height+i;
      break;
    }
  }

  node->cr_minage = cred_lower * interval_age;
  node->cr_maxage = cred_upper * interval_age;
}

static double * create_child_probs(tree_node_t * node, double * parent_prob)
{
  long i,j;
  tree_node_t * parent;
  double abs_age_node, abs_age_parent;
  double prob_rate;
  double table_sum;

  assert(node->parent);

  /* if tip (and not fossil) return */
  if (!node->left && !node->prior) return NULL;

  parent = node->parent;

  double * table = (double *)xmalloc((size_t)(node->entries) * sizeof(double));
  double * logprob = (double *)xmalloc((size_t)(node->entries) * sizeof(double));
  int * mask = (int *)xmalloc((size_t)(node->entries) * sizeof(int));
  memset(mask,0,node->entries * sizeof(int));

  /* go through the entries of the parent */
  for (i = 0; i < parent->entries; ++i)
  {
    abs_age_parent = (parent->height+i)*interval_age;

    /* get number of allowed entries for child node given i (line of parent) */
    long jmax = (parent->height + i - node->height);
    if (jmax > node->entries)
      jmax = node->entries;

    /* traverse allowed entries of child and obtain unnormalized probabilities
       for each line, and their sum */
    for (j = 0; j < jmax; ++j)
    {
      abs_age_node = (node->height+j)*interval_age;
      prob_rate = gamma_dist_logpdf(node->length /
                                    (abs_age_parent - abs_age_node));
      table[j] = node->matrix_sum[j] + prob_rate;

      if (j)
        table_sum = logsum(table_sum, node->matrix_sum[j] + prob_rate);
      else
        table_sum = table[j];
    }

    /* go through the table, normalize and multiply with parent probability */
    for (j = 0; j < jmax; ++j)
    {
      if (mask[j])
        logprob[j] = logsum(logprob[j], table[j] - table_sum + parent_prob[i]);
      else
        logprob[j] = table[j] - table_sum + parent_prob[i];

      /* mark that entry j of table and logprob is initialized to a value,
         required for hte logsum function */
      mask[j] = 1;
    }
  }

  free(mask);
  free(table);

  return logprob;
}

static void credible_recursive(tree_node_t * node, double * node_logprob)
{
  double * child_logprob;

  if (!node) return;
  if (!node->left && !node->prior) return;

  /* update credible interval for node */
  cr_update(node, node_logprob);

  /* compute grid line probabilities for the left child and recursively process
     left subtree */
  child_logprob = create_child_probs(node->left, node_logprob);
  credible_recursive(node->left, child_logprob);

  free(child_logprob);

  /* compute grid line probabilities for the right child and recursively process
     right subtree */
  child_logprob = create_child_probs(node->right, node_logprob);
  credible_recursive(node->right, child_logprob);

  free(child_logprob);
}

void credible(tree_node_t * root)
{
  double norm = 0;
  long i;

  double * logprob = (double *)xmalloc((size_t)(root->entries) * sizeof(double));

  interval_age = opt_max_age / (opt_grid_intervals - 1);

  /* compute normalization constant for root */
  norm = root->matrix_sum[0];
  for (i = 1; i < root->entries; ++i)
  {
    norm = logsum(norm, root->matrix_sum[i]);
  }

  /* normalize log probabilities for root */
  for (i = 0; i < root->entries; ++i)
    logprob[i] = root->matrix_sum[i] - norm;

  /* process tree starting from root */
  credible_recursive(root, logprob);

  free(logprob);
}
