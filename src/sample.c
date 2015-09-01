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
    int i; /*TODO: Is this supposed to be a long or an int?*/
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

  long entries = parent->sampled_gridline - node->height;
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
  long i;

  char * filename = (char *)xmalloc((strlen(opt_outfile)+9)*sizeof(char));
  double interval_age = opt_max_age / (opt_grid_intervals - 1);
  if (opt_method_relative)
    interval_age =1;
  strcpy(filename, opt_outfile);
  strcat(filename,".sampled");

  FILE * fp_out = fopen(filename, "w");

  progress_init("Sampling...", (unsigned long)opt_sample);
  for (i = 0; i < opt_sample; ++i)
  {
    dp_backtrack_sampling(root);
    output_sample_tree_recursive(root, interval_age, fp_out);
    fprintf(fp_out, ";\n");
    progress_update((unsigned long)i);
  }
  progress_done();

  fclose(fp_out);
  free(filename);
}



static void dp_backtrack_interval_recursive(tree_node_t * node)
{
  double maxscore;

  if (!node) return;

  assert(node->parent);

  tree_node_t * parent = node->parent;
  int par_lowerbound = parent->lowerbound;
  int par_upperbound = parent->upperbound;
  int i,z,x;
  double * mult_vector = (double *)xmalloc((size_t)opt_grid_intervals * sizeof(double));


  long poss_entries = (par_upperbound + parent->height) - node->height; /*the maximum range of locations for the child*/
  if (poss_entries > node->entries)
    poss_entries = node->entries;

  for (z = 0; z < poss_entries; ++z) /*TODO is this entirely uneccessary? is vector allocated as 0's*/
          mult_vector[z] = 0;

  for (i = par_lowerbound; i < par_upperbound; ++i) /*loops through 95% comfidence intevral for parent*/
  { 
      double rel_age_parent = (1.0 / opt_grid_intervals) * (i + parent->height);
      double par_factor = parent->interval_weights[i - par_lowerbound];
      long entries = (i + parent->height) - node->height;
      if (entries > node->entries)
           entries = node->entries; /*Should probbaly special case the tips, because this will loop through all parent pos even though there si only one possible position..;.*/
      recompute_scores(node, rel_age_parent, cdf_vector, entries, &maxscore); /*compute the scores given this parental position*/
      normalize_cdf(cdf_vector, entries, maxscore);  /*rescale */
      for (x = 0; x < entries; ++x)
        {
          assert(x < poss_entries);
/*          printf("node %s, parent is %s, parent at %i, child at %i, par weight is %f, child weight is %f\n", node->label, parent->label, i, z, parent->interval_weights[i - par_lowerbound],  cdf_vector[z]);*/
          assert(cdf_vector[x] >= 0);
          double local_prob;
          if (x > 0)
            local_prob = cdf_vector[x] - cdf_vector[x-1]; /*de-cumulate*/
          else
          {
            local_prob = cdf_vector[x];
          }
          mult_vector[x] = ((local_prob * par_factor) +  mult_vector[x]); /*sum the probability weight fo those positions, geven the weight of the parent being at that position, this is not a cumulative probability density!*/
          assert(mult_vector[x]>=0);
        } 
    }

    /*Now re-normalize the vector across probabilities TODO risk of underflow?*/
    double mult_vector_sum = 0;
    for (i = 0; i < poss_entries; ++i)
    {
       mult_vector_sum += mult_vector[i];
    }
    if (poss_entries == 1) /*should probably handle tips better...*/
      {
        mult_vector[poss_entries-1] = 1;
      }
    else
       { mult_vector[0] /= mult_vector_sum;
        for (i = 1; i < poss_entries; ++i)
        {
          mult_vector[i] = (mult_vector[i] / mult_vector_sum) + mult_vector[i-1];
        }
      }
    assert(0.9999 < mult_vector[poss_entries-1]);
    assert(mult_vector[poss_entries-1] < 1.0001);

   node->lowerbound = -1;
   node->upperbound = -1;
   double upper_thresh = mult_vector[node->interval_line - node->height] + (opt_conf_interval/2); /*cumulative probability density of the upper conf*/
   double lower_thresh = mult_vector[node->interval_line - node->height] - (opt_conf_interval/2); /*cumulative probability density of the lower conf*/
   if (upper_thresh > 1)
        {
          node->upperbound = poss_entries; /* TODO Confidence interval of node label %s exceeded the upper bound, how should this be handled?*/
          printf("Confidence intervals hitting bound, \n");
        }
   if (lower_thresh < 0)
        {
          node->lowerbound = 0; /* TODO Confidence interval of node label %s exceeded the lower bound, how should this be handled?*/
        }
   for (i = 0; i < poss_entries; ++i)
    {   
        if (node->entries == 1)
        {
            node->lowerbound = 0;
            node->upperbound = 0;
        }
        if ((mult_vector[i] > lower_thresh) && (node->lowerbound == -1))
          node->lowerbound = i;
        if ((mult_vector[i] > upper_thresh) && (node->upperbound == -1))
          node->upperbound = i;
    }
    assert(node->lowerbound != -1);
    assert(node->upperbound != -1);
  node->interval_weights[0] = 0;
  for (i = node->lowerbound + 1; i < node->upperbound; ++i)
  {  
     node->interval_weights[i - node->lowerbound] = (mult_vector[i]-mult_vector[i-1])/opt_conf_interval;   /*Making it non-cumulative, TODO off by one?!*/
  }

  dp_backtrack_interval_recursive(node->left);
  dp_backtrack_interval_recursive(node->right);
  free(mult_vector);

}



static void dp_backtrack_interval(tree_node_t * root)
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

  root->lowerbound = -1;
  root->upperbound = -1;
  double upper_thresh = cdf_vector[root->interval_line - root->height] + (opt_conf_interval/2); /*cumulative probability density of the upper conf, centered on best estimate...*/
  double lower_thresh = cdf_vector[root->interval_line - root->height] - (opt_conf_interval/2); /*cumulative probability density of the lower conf, centered on best estimate..*/
  if (upper_thresh > 1)
        {
          root->upperbound = root->entries; /* TODO Confidence interval of node label %s exceeded the upper bound, how should this be handled?*/
          printf("Root confidence interval hit maximum interval line value\n");
        }
  if (lower_thresh < 0)
        {
          root->lowerbound = 0; /* TODO Confidence interval of node label %s exceeded the lower bound, how should this be handled?*/
          printf("Root confidence interval hit minimum interval line value\n");
        }
  for (i = 0; i < root->entries; ++i)
    { 
        if ((cdf_vector[i] > lower_thresh) && (root->lowerbound == -1))
        {
          root->lowerbound = i;
        }
        if ((cdf_vector[i] > upper_thresh) && (root->upperbound == -1))
        {
          root->upperbound = i;
        }
    }  
    assert(root->lowerbound != -1);
    assert(root->upperbound != -1);
  root->interval_weights[0] = 0;
  for (i = root->lowerbound + 1; i < root->upperbound; ++i)
  {  
    root->interval_weights[i - root->lowerbound] = (cdf_vector[i]-cdf_vector[i-1])/opt_conf_interval;   /*Making it non-cumulative, TODO off by one?!*/
    assert(root->interval_weights[i - root->lowerbound] >= 0);
  }

  dp_backtrack_interval_recursive(root->left);
  dp_backtrack_interval_recursive(root->right);

  free(cdf_vector);
}


/*static void output_intervals_tree_recursive(tree_node_t * node,
                                         double interval_age,
                                         FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s[&age=%f-%f]:%f",
            node->label, (node->lowerbound + node->height) * interval_age, (node->upperbound + node->height) * interval_age, node->length);
  else
  {
    fprintf(fp_out, "(");
    output_intervals_tree_recursive(node->left, interval_age, fp_out);
    fprintf(fp_out, ",");
    output_intervals_tree_recursive(node->right, interval_age, fp_out);
    fprintf(fp_out, ")%s[&age=%f-%f]:%f", node->label ? node->label : "",
                    (node->lowerbound + node->height) * interval_age, (node->upperbound + node->height) * interval_age, node->length);
  }
}*/


void interval(tree_node_t * root)
{
/*  char * filename = (char *)xmalloc((strlen(opt_outfile)+9)*sizeof(char));
  double interval_age = opt_max_age / (opt_grid_intervals - 1);
  if (opt_method_relative)
    interval_age =1;
  strcpy(filename, opt_outfile)
  strcat(filename,".intervals");

  FILE * fp_out = fopen(filename, "w");*/
  dp_backtrack_interval(root);
/*  output_intervals_tree_recursive(root, interval_age, fp_out);*/

/*  fclose(fp_out);
  free(filename);*/
}

 