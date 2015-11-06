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
    double score;
    double rel_age_node;
    double prob_rate_node;
    long node_low = node->height;

    *maxscore = -__DBL_MAX__;

    for (i = 0; i < count; ++i)
    {
      rel_age_node = (1.0 / opt_grid_intervals) * (node_low+i);

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
                                               double interval_age,
                                               FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s[&age=%f]:%f",
            node->label, node->sampled_gridline * interval_age, node->length);
  else
  {
    fprintf(fp_out, "(");
    output_sample_dated_tree_recursive(node->left, interval_age, fp_out);
    fprintf(fp_out, ",");
    output_sample_dated_tree_recursive(node->right, interval_age, fp_out);
    fprintf(fp_out, ")%s[&age=%f]:%f", node->label ? node->label : "",
                    node->sampled_gridline * interval_age, node->length);
  }
}

static void output_sample_um_tree_recursive(tree_node_t * node, 
                                             double interval_age,
                                             FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s:%f", node->label, 
                    (node->parent->sampled_gridline - node->sampled_gridline) *
                      interval_age);
  else
  {
    fprintf(fp_out, "(");
    output_sample_um_tree_recursive(node->left, interval_age, fp_out);
    fprintf(fp_out, ",");
    output_sample_um_tree_recursive(node->right, interval_age, fp_out);
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
  double interval_age;

  char * filename = (char *)xmalloc((strlen(opt_outfile)+9)*sizeof(char));

  if (opt_max_age)
    interval_age = opt_max_age / (opt_grid_intervals - 1);
  else
    interval_age = 1.0 / (opt_grid_intervals - 1);

  strcpy(filename, opt_outfile);
  strcat(filename,".sampled");

  FILE * fp_out = fopen(filename, "w");

  progress_init("Sampling...", (unsigned long)opt_sample);
  for (i = 0; i < opt_sample; ++i)
  {
    dp_backtrack_sampling(root);
    if (opt_outform == OUTPUT_ULTRAMETRIC)
      output_sample_um_tree_recursive(root, interval_age, fp_out);
    else if (opt_outform == OUTPUT_DATED)
      output_sample_dated_tree_recursive(root, interval_age, fp_out);
    else
      fatal("Internal error while selecting output format");
    fprintf(fp_out, ";\n");
    progress_update((unsigned long)i);
  }
  progress_done();

  fclose(fp_out);
  free(filename);
}
