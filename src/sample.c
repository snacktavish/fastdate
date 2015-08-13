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

static int sample_gridline(double * vector, int count) 
{
  int i;

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
                             int count,
                             double * sumscore)
{
    int i;
    double score;
    double rel_age_node;
    double prob_rate_node;
    unsigned int node_low = node->height;

    *sumscore = 0;/*sumscore and score in probability space, NOT LOG*/

    for (i = 0; i < count; ++i)
    {
      rel_age_node = (1.0 / opt_grid_intervals) * (node_low+i);

      prob_rate_node = gamma_dist_logpdf(node->length / 
                                         (rel_age_parent - rel_age_node));

      score = node->matrix_PP[i] + exp(prob_rate_node);
      *sumscore += score;
      assert (sumscore > 0);
      vector[i] = score;
    }
}

static void normalize_cdf(tree_node_t * node,
                          double * vector,
                          int count,
                          double sumscore)
{
    int i;

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
  double sumscore;

  if (!node) return;
  assert(node->parent);

  tree_node_t * parent = node->parent;

  double rel_age_parent = (1.0 / opt_grid_intervals) *
                          (parent->sampled_gridline);

  int entries = (node->left) ? 
                (parent->sampled_gridline - node->height) : 1;

  assert(entries <= node->entries);

  recompute_scores(node, rel_age_parent, cdf_vector, entries, &sumscore); /*HALP what is in the cdf vector here?! Doesn't matter bc it is overwritten...*/
  normalize_cdf(node, cdf_vector, entries, sumscore);
  node->sampled_gridline = node->height + 
                           sample_gridline(cdf_vector, entries);

  dp_backtrack_sampling_recursive(node->left);
  dp_backtrack_sampling_recursive(node->right);
}

static void dp_backtrack_sampling(tree_node_t * root)
{
  int i;
  double sumscore = 0;

  cdf_vector = (double *)xmalloc(opt_grid_intervals * sizeof(double));

  /* copy logscores to cdf vector and get max */
  memcpy(cdf_vector, root->matrix_PP, root->entries * sizeof(double));
  for (i = 0; i < root->entries; ++i)
    sumscore += cdf_vector[i];

  /* normalize logscores according to maxscore */
  normalize_cdf(root, cdf_vector, root->entries, sumscore);

  /* select interval line for root */
  root->sampled_gridline = root->height + 
                           sample_gridline(cdf_vector, root->entries);

  dp_backtrack_sampling_recursive(root->left);
  dp_backtrack_sampling_recursive(root->right);

  free(cdf_vector);
}

static void output_sample_tree_recursive(tree_node_t * node,
                                         double interval_age,
                                         FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s[&age=%f]:%f",
            node->label, node->sampled_gridline * interval_age, node->length);
  else
  {
    fprintf(fp_out, "(");
    output_sample_tree_recursive(node->left, interval_age, fp_out);
    fprintf(fp_out, ",");
    output_sample_tree_recursive(node->right, interval_age, fp_out);
    fprintf(fp_out, ")%s[&age=%f]:%f", node->label ? node->label : "",
                    node->sampled_gridline * interval_age, node->length);
  }
}

void sample(tree_node_t * root)
{
  int i;

  char * filename = (char *)xmalloc((strlen(opt_outfile)+9)*sizeof(char));
  double interval_age = opt_max_age / (opt_grid_intervals - 1);
  if (opt_method_relative)
      interval_age = 1;
  strcpy(filename, opt_outfile);
  strcat(filename,".sampled");

  FILE * fp_out = fopen(filename, "w");

  progress_init("Sampling...", opt_sample);
  for (i = 0; i < opt_sample; ++i)
  {
    dp_backtrack_sampling(root);
    output_sample_tree_recursive(root, interval_age, fp_out);
    fprintf(fp_out, ";\n");
    progress_update(i);
  }
  progress_done();

  fclose(fp_out);
  free(filename);
}
