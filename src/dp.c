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

/* L-BFGS-B optimization factor
 * (relative to the machine epsilon)
 */
static const double opt_factor = 1e12;
/* gradient tolerance */
static const double opt_pgtol = 0.1;
/* epsilon for the iterative optimization */
static const double opt_epsilon = 1.0e-2;

static long inner_entries = 0;
static double interval_age = 0;
static long tree_height = 0;

static pthread_attr_t attr;

typedef struct thread_info_s
{
  pthread_t thread;
  pthread_mutex_t mutex;
  pthread_cond_t cond;
  int work;
  tree_node_t * node;
  long entry_first;
  long entry_count;
} thread_info_t;

static thread_info_t * ti;

static double age_prior_logpdf(tree_node_t * node, double abs_age_node)
{
  double dist_logprob;

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
  else if (node->prior == NODEPRIOR_UNI)
  {
    uni_params_t * params = (uni_params_t *)(node->prior_params);
    dist_logprob = uni_dist_logpdf(params->min_age,
                                   params->max_age,
                                   abs_age_node);
  }
  else if (node->prior == NODEPRIOR_NORM)
  {
    norm_params_t * params = (norm_params_t *)(node->prior_params);
    dist_logprob = norm_dist_logpdf(params->mean,
                                    params->variance,
                                    abs_age_node - params->offset);
  }
  else if (node->prior)
    fatal("Unknown prior used");

  return dist_logprob;
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

static void fill_tip_table(tree_node_t * node)
{

  long i;
  double abs_age_node;
  double dist_logprob;

  if (!node->prior)
  {
    node->matrix[0] = 0.0;    /* 0 and not 1 because of the log-scale */
    node->matrix_sum[0] = 0;
    return;
  }

  /* tip fossil case*/
  for (i = 0; i < node->entries; ++i)
  {
    abs_age_node = (node->height+i)*interval_age;

    /* if estimating absolute ages check for node priors and compute PDF */
    dist_logprob = age_prior_logpdf(node, abs_age_node);
    
    node->matrix[i] = dist_logprob + bd_tipdates_prod_tip(abs_age_node);
    node->matrix_sum[i] = 0;
  }
}

static void fill_inner_table(tree_node_t * node, long entry_first, long entry_last)
{
  long i,j,k;
  long jmax, kmax;
  long low;
  long left_low;
  long right_low;
  
  tree_node_t * left;
  tree_node_t * right;

  double age_diff;
  double abs_age_node;
  double abs_age_left;
  double abs_age_right;
  double prob_rate_left;
  double prob_rate_right;
  double dist_logprob;
  double score;
  double score_lsum;
  double score_rsum;

  left  = node->left;
  right = node->right;

  low = node->height;
  left_low = left->height;
  right_low = right->height;

  /* run DP */

  /*

                        o  (abs_age_node)
                       / \
                      /   \
      (abs_age_left) o     \
                            o (abs_age_right)
  */

  for (i = entry_first; i < entry_last; ++i)
  {
    abs_age_node = (low+i)*interval_age;
    score_lsum = score_rsum = 0;

    /* check the ages of the left child */
    if (!left->left)
      jmax = 1;
    else
      jmax = (node->height + i - left->height);

    /* in case we had a uniform prior. TODO: Check whether this should also
       replace 'jmax = 1' three lines above for tip fossils */
    if (jmax > left->entries) jmax = left->entries;

    assert(jmax <= left->entries);

    long jbest = -1;
    double jbest_score = -__DBL_MAX__;
    for (j = 0; j < jmax; ++j)
    {
      assert(j+left->height < node->height + i);
      abs_age_left = (left_low+j)*interval_age;

      prob_rate_left = gamma_dist_logpdf(left->length / 
                                         (abs_age_node - abs_age_left));
      score = left->matrix[j] + prob_rate_left;

      if (j)
        score_lsum = logsum(score_lsum, left->matrix_sum[j] + prob_rate_left);
      else
        score_lsum = left->matrix_sum[j] + prob_rate_left;

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

    /* in case we had a uniform prior. TODO: Check whether this should also
       replace 'kmax = 1' three lines above for tip fossils */
    if (kmax > right->entries) kmax = right->entries;

    assert(kmax <= right->entries);

    long kbest = -1;
    double kbest_score = -__DBL_MAX__;
    for (k = 0; k < kmax; ++k)
    {
      assert(k+right->height < node->height + i);
      abs_age_right = (right_low+k)*interval_age;

      age_diff = abs_age_node - abs_age_right;
      prob_rate_right = gamma_dist_logpdf(right->length / age_diff);

      score = right->matrix[k] + prob_rate_right;

      if (k)
        score_rsum = logsum(score_rsum, right->matrix_sum[k] + prob_rate_right);
      else
        score_rsum = right->matrix_sum[k] + prob_rate_right;
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
      bd_term = bd_relative_prod(abs_age_node);
    else if (opt_method_tipdates)
      bd_term = bd_tipdates_prod_inner(abs_age_node);
    else assert(0);


    dist_logprob = age_prior_logpdf(node, abs_age_node);

    score = bd_term + jbest_score + kbest_score + dist_logprob;
    double score_sum = bd_term + score_lsum + score_rsum + dist_logprob;

    /* if it's the root add one more term */
    if (node->height == tree_height)
    {
      if (opt_method_relative || opt_method_nodeprior)
      {
        score += bd_relative_root(node->leaves,
                                  abs_age_node);
        score_sum += bd_relative_root(node->leaves,
                                      abs_age_node);
      }
      else if (opt_method_tipdates)
      {
        score += bd_tipdates_root(node->leaves,
                                  abs_age_node);
        score_sum += bd_tipdates_root(node->leaves,
                                      abs_age_node);
      }
      else assert(0);
    }

    /* store best placement of children and likelihood for interval line i */
    node->matrix[i] = score;
    node->matrix_sum[i] = score_sum;
    node->matrix_left[i] = jbest;
    node->matrix_right[i] = kbest;
  }
}

static inline void dp_recurse_worker(long t)
{
  thread_info_t * tip = ti + t;
  tree_node_t * node  = tip->node;

  fill_inner_table(node, tip->entry_first, tip->entry_first + tip->entry_count);
  
}

static void * threads_worker(void *vp)
{
  long t = (long) vp;
  thread_info_t * tip = ti + t;

  pthread_mutex_lock(&tip->mutex);

  /* loop until signalled to quit */
  while (tip->work >= 0)
  {
    /* wait for work available */
    if (tip->work == 0)
      pthread_cond_wait(&tip->cond, &tip->mutex);
    if (tip->work > 0)
    {
      dp_recurse_worker(t);
      tip->work = 0;
      pthread_cond_signal(&tip->cond);
    }
  }
  pthread_mutex_unlock(&tip->mutex);

  pthread_exit(NULL);
}

static void threads_init()
{
  long t;

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* allocate memory for thread info */
  ti = (thread_info_t *)xmalloc((size_t)opt_threads * sizeof(thread_info_t));

  /* init and create worker threads */
  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;
    tip->work = 0;
    pthread_mutex_init(&tip->mutex,NULL);
    pthread_cond_init(&tip->cond,NULL);
    if (pthread_create(&tip->thread, &attr, threads_worker, (void *)(long)t))
      fatal("Cannot create thread");
  }
}

static void threads_wakeup(tree_node_t * node)
{
  long t;
  long entries = node->entries;
  long threads = entries > opt_threads ? opt_threads : entries;
  long entries_rest = entries;
  long entry_next = 0;


  /* precompute total workload as the sum of operations required for each
     gridline node can be placed on */
  double term = 2 * node->height + node->left->height + node->right->height;
  unsigned long workload = (unsigned long)entries * (unsigned long)term -
                           (unsigned long)entries +
                           (unsigned long)(entries*entries);
  double wl_thread = workload / (double)threads;

  entries_rest = entries;
  entry_next = 0;


  /* tell the threads there is work to do */
  for (t = 0; t < threads; ++t)
  {
    thread_info_t * tip = ti + t;

    tip->node = node;
    tip->entry_first = entry_next;

    /* compute the workload for thread t by solving the quadratic equation */
    if (t != threads-1)
      tip->entry_count = (long)((sqrt((term-1)*(term-1) + 4.0*wl_thread*(t+1)) +
                               1 - term ) / 2 - tip->entry_first);
    else
      tip->entry_count = entries_rest;

    entries_rest -= tip->entry_count;
    entry_next += tip->entry_count;

    pthread_mutex_lock(&tip->mutex);
    tip->work = 1;
    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);
  }

  /* wait for threads to finish their work */
  for (t = 0; t < threads; ++t)
  {
    thread_info_t * tip = ti+t;
    pthread_mutex_lock(&tip->mutex);
    while(tip->work > 0)
      pthread_cond_wait(&tip->cond, &tip->mutex);
    pthread_mutex_unlock(&tip->mutex);
  }
}

static void threads_exit()
{
  long t;

  for (t = 0; t < opt_threads; ++t)
  {
    thread_info_t * tip = ti + t;

    /* tell worker to quit */
    pthread_mutex_lock(&tip->mutex);
    tip->work = -1;
    pthread_cond_signal(&tip->cond);
    pthread_mutex_unlock(&tip->mutex);

    /* wait for worker to quit */
    if (pthread_join(tip->thread, 0))
      fatal("Cannot join thread");

    pthread_cond_destroy(&tip->cond);
    pthread_mutex_destroy(&tip->mutex);
  }

  free(ti);
  pthread_attr_destroy(&attr);
}

static void dp_recurse_parallel(tree_node_t * node)
{
  static long sum_entries = 0;


  if (!node) return;

  /* leaves case */
  if (!node->left)
  {
    fill_tip_table(node);
    return;
  }

  /* inner nodes case */
  dp_recurse_parallel(node->left);
  dp_recurse_parallel(node->right);

  /* schedule work for threads for current node */
  threads_wakeup(node);

  sum_entries += node->entries;
  progress_update((unsigned long)sum_entries);
}

/* this resets the node heights for method_nodeprior or method_tipdates */
static void reset_node_heights(tree_node_t * node)
{
  assert(node);

  double offset = 0;
  long min_height = 0;

  /* number of lines required above the new min height which will
     be selected for the current node */
  long required_space = tree_height - node->height;

  /* compute minimum possible interval line given a node has calibration info */
  if (node->prior == NODEPRIOR_EXP)
    offset = ((exp_params_t *)(node->prior_params))->offset;
  else if (node->prior  == NODEPRIOR_LN)
    offset = ((ln_params_t *)(node->prior_params))->offset;
  else if (node->prior == NODEPRIOR_NORM)
    offset = ((norm_params_t *)(node->prior_params))->offset;
  else if (node->prior == NODEPRIOR_UNI)
    offset = ((uni_params_t *)(node->prior_params))->min_age;
  min_height = lrint(ceil(offset / interval_age));

  /* tip case */
  if (!node->left)
  {
    node->height = min_height;
    if (node->height + required_space >= opt_grid_intervals)
      fatal("Error: Calibration info on line %d of file %s disrupts grid",
            node->prior_lineno, opt_priorfile);
    return;
  }

  /* inner node case */
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
  long entries;

  assert(node);

  /* tip case */
  if (!node->left)
  {
    if (!node->prior)
    {
      /* tips are always placed on the first grid line */
      node->entries      = 1;
      node->matrix       = (double *)xmalloc(sizeof(double));
      node->matrix_sum   = (double *)xmalloc(sizeof(double));
      node->matrix_left  = (long *)xmalloc(sizeof(long));
      node->matrix_right = (long *)xmalloc(sizeof(long));
      node->interval_weights = (double *)xmalloc(sizeof(double));
    }
    else
    {
      entries = node->parent->entries + node->parent->height - node->height - 1;

      /* if uniform distribution we need to check for max age */
      if (node->prior == NODEPRIOR_UNI)
      {
        double max_age = ((uni_params_t *)(node->prior_params))->max_age;
        long max_grid_line = lrint(ceil(max_age / interval_age)); 

        assert(max_grid_line >= node->height);

        if (node->height + entries - 1 > max_grid_line)
          entries = max_grid_line - node->height + 1;
      }
      node->entries      = entries;
      node->matrix       = (double *)xmalloc((size_t)entries * sizeof(double));
      node->matrix_sum   = (double *)xmalloc((size_t)entries * sizeof(double));
      node->matrix_left  = (long *)xmalloc((size_t)entries * sizeof(long));
      node->matrix_right = (long *)xmalloc((size_t)entries * sizeof(long));
      node->interval_weights = (double *)xmalloc((size_t)entries * sizeof(double));
    }
    return;
  }

  if (!node->parent)
    entries = opt_grid_intervals - node->height;
  else
    entries = node->parent->entries + node->parent->height - node->height - 1;

  /* if uniform distribution we need to check for max age */
  if (node->prior == NODEPRIOR_UNI)
  {
    double max_age = ((uni_params_t *)(node->prior_params))->max_age;
    long max_grid_line = lrint(ceil(max_age / interval_age)); 

    assert(max_grid_line >= node->height);

    if (node->height + entries - 1 > max_grid_line)
      entries = max_grid_line - node->height + 1;
  }

  /* allocate storage space for placement information at each node */
  node->entries      = entries;
  node->matrix       = (double *)xmalloc((size_t)entries * sizeof(double));
  node->matrix_sum   = (double *)xmalloc((size_t)entries * sizeof(double));
  node->matrix_left  = (long *)xmalloc((size_t)entries * sizeof(long));
  node->matrix_right = (long *)xmalloc((size_t)entries * sizeof(long));
  node->interval_weights = (double *)xmalloc((size_t)entries * sizeof(double));

  /* for progress bar indication */
  inner_entries += entries;

  /* allocate the space for its two subtrees */
  alloc_node_entries(node->left);
  alloc_node_entries(node->right);

}

static void dp_recurse_serial(tree_node_t * node)
{
  static long sum_entries = 0;

  tree_node_t * left;
  tree_node_t * right;

  if (!node) return;

  /* leaves case */
  if (!node->left)
  {
     fill_tip_table(node);
     return;
  }

  /* inner nodes case */
  left  = node->left;
  right = node->right;

  dp_recurse_serial(left);
  dp_recurse_serial(right);

  fill_inner_table(node,0,node->entries);
  sum_entries += node->entries;
  progress_update((unsigned long)sum_entries);
}

static void dp_backtrack_recursive(tree_node_t * node, long best_entry)
{
  if (!node) return;
  node->interval_line = node->height + best_entry;

  if (node->left)
  {
    dp_backtrack_recursive(node->left, node->matrix_left[best_entry]);
    dp_backtrack_recursive(node->right, node->matrix_right[best_entry]);
  }
}

static double dp_backtrack(tree_node_t * root)
{
  long i;
  long best_entry = 0;
  double max_score = -__DBL_MAX__;
  double score = 0.0;

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

  return max_score;
}

double dp_evaluate(tree_node_t * tree)
{
  double score = 0.0;

  gamma_dist_init();
  bd_init(0, 0);

  if (opt_threads > 1)
  {
    //TODO: assert that pthreads is initialized
    dp_recurse_parallel(tree);
  }
  else
    dp_recurse_serial(tree);

  score = dp_backtrack(tree);

  return score;
}

void dp(tree_node_t * tree)
{
  double score;

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

  progress_init("Running DP...", (unsigned long)inner_entries);
  if (opt_threads > 1)
  {
    threads_init();
    dp_recurse_parallel(tree);
    threads_exit();
  }
  else
    dp_recurse_serial(tree);

  progress_done();

  if (!opt_quiet)
    printf ("Backtracking...\n");

  score = dp_backtrack(tree);
  
  /* optimize parameters */
  if (opt_parameters_bitv)
  {
    printf("\n\n*** Parameters optimization ***\n\n");

    opt_quiet = 1;

    printf("Starting parameters:\n");
    if (opt_parameters_bitv & PARAM_LAMBDA)
      printf(" lambda:    %6.4f\n", opt_lambda);
    if (opt_parameters_bitv & PARAM_MU)
      printf(" mu:        %6.4f\n", opt_mu);
    if (opt_parameters_bitv & PARAM_RHO)
      printf(" rho:       %6.4f\n", opt_rho);
    if (opt_parameters_bitv & PARAM_PSI)
      printf(" psi:       %6.4f\n", opt_psi);
    if (opt_parameters_bitv & PARAM_RATE_MEAN)
    {
      printf(" rate mean: %6.4f\n", opt_rate_mean);
      if (opt_fixgamma)
        printf(" rate var:  %6.4f\n", opt_rate_var);
    }
    if (opt_parameters_bitv & PARAM_RATE_VAR)
    {
      printf(" rate var:  %6.4f\n", opt_rate_var);
    }
    printf("\n");

    printf("Starting score: %f\n\n", score);

    /* single param iterative optimization */
    if (opt_threads > 1)
      threads_init();
    double cur_score = score + 1;
    int i = 0;
    while (fabs(cur_score - score) > opt_epsilon)
    {
      printf("[Iteration %d]\n", ++i);
      cur_score = score;
      //opt_parameters(tree, PARAM_PSI, opt_factor, opt_pgtol);
      if (opt_parameters_bitv & PARAM_LAMBDA)
      {
        score = opt_parameters(tree, PARAM_LAMBDA, opt_factor, opt_pgtol, cur_score);
        printf("%15.4f lambda:    %6.4f\n", score, opt_lambda);
      }
      if (opt_parameters_bitv & PARAM_MU)
      {
        score = opt_parameters(tree, PARAM_MU, opt_factor, opt_pgtol, cur_score);
        printf("%15.4f mu:        %6.4f\n", score, opt_mu);
      }
      if (opt_parameters_bitv & PARAM_PSI)
      {
        score = opt_parameters(tree, PARAM_PSI, opt_factor, opt_pgtol, cur_score);
        printf("%15.4f psi:       %6.4f\n", score, opt_psi);
      }
      if (opt_parameters_bitv & PARAM_RHO)
      {
        score = opt_parameters(tree, PARAM_RHO, opt_factor, opt_pgtol, cur_score);
        printf("%15.4f rho:       %6.4f\n", score, opt_rho);
      }
      if (opt_parameters_bitv & PARAM_RATE_MEAN)
      {
        score = opt_parameters (tree, PARAM_RATE_MEAN, opt_factor, opt_pgtol, cur_score);
        printf ("%15.4f rate mean: %6.4f\n", score, opt_rate_mean);
      }
      if (opt_parameters_bitv & PARAM_RATE_VAR)
      {
        score = opt_parameters(tree, PARAM_RATE_VAR, opt_factor, opt_pgtol, cur_score);
        printf("%15.4f rate var:  %6.4f\n", score, opt_rate_var);
      }
    }
    if (opt_threads > 1)
      threads_exit();

    printf("\nFinal parameters:\n");
    if (opt_parameters_bitv & PARAM_LAMBDA)
      printf(" lambda:    %6.4f\n", opt_lambda);
    if (opt_parameters_bitv & PARAM_MU)
      printf(" mu:        %6.4f\n", opt_mu);
    if (opt_parameters_bitv & PARAM_RHO)
      printf(" rho:       %6.4f\n", opt_rho);
    if (opt_parameters_bitv & PARAM_PSI)
      printf(" psi:       %6.4f\n", opt_psi);
    if (opt_parameters_bitv & PARAM_RATE_MEAN)
    {
      printf (" rate mean: %6.4f\n", opt_rate_mean);
    }
    if (opt_parameters_bitv & PARAM_RATE_VAR)
    {
      printf (" rate var:  %6.4f\n", opt_rate_var);
    }

    printf("\nScore after optimization %f\n", score);
    printf("\n*** ***** ************ ***\n\n");
  }
}
