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



double check_prior(tree_node_t * node, double abs_age_node)
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
    else if (node->prior)
      assert(0);
    return dist_logprob;
}

static void calc_tip(tree_node_t * node)
{

  long i;
  double rel_age_node;
  double abs_age_node;
  double dist_logprob;
  if (!node->prior)
  {
    node->matrix[0] = 0.0;    /* 0 and not 1 because of the log-scale */
    node->matrix_PP[0] = 0.0;
    return;
  }

  /* tip fossil case*/
  for (i = 0; i < node->entries; ++i)
  {
    rel_age_node = (1.0 / opt_grid_intervals) * (node->height+i);
    abs_age_node = (node->height+i)*interval_age;
    /* if estimating absolute ages check for node priors and compute PDF */
    dist_logprob = 0;
    dist_logprob = check_prior(node, abs_age_node);
    
    node->matrix[i] = dist_logprob + bd_tipdates_prod_tip(rel_age_node);
    node->matrix_PP[i] = dist_logprob + bd_tipdates_prod_tip(rel_age_node);
  }
}

double calc_log_sum(double  prev_ln_sum, double new_ln_term) /*swiped from slowdate.py*/
 /*Avoid underflow/overflow, but return:
      log(x + y)
  where x = exp(prev_ln_sum) and y = exp(new_ln_term)
  */
{
  double ln_m, exp_sum;
  if ((new_ln_term - prev_ln_sum) > 30)
    return new_ln_term; /* first term, or previous is lost in rounding error*/
  if ( prev_ln_sum - new_ln_term > 30)
    return prev_ln_sum; /* new term is lost in rounding error*/
  if (prev_ln_sum < new_ln_term)
    ln_m = prev_ln_sum;
  else
    ln_m = new_ln_term;
  prev_ln_sum -= ln_m;
  new_ln_term -= ln_m;
  /* one is now 0 and the other is no greater than 30*/
  exp_sum = exp(prev_ln_sum) + exp(new_ln_term);
  return (log(exp_sum) + ln_m);
}


static void calc_interval_scores(tree_node_t * node, long i_min, long i_max)
{
  long i,j,k;
  long jmax, kmax;
  long low;
  long left_low;
  long right_low;
  
  tree_node_t * left;
  tree_node_t * right;

  double rel_age_node;
  double rel_age_left;
  double rel_age_right;
  double age_diff;
  double abs_age_node;
  double abs_age_left;
  double abs_age_right;
  double prob_rate_left;
  double prob_rate_right;
  double dist_logprob;
  double dist_logprob_left;
  double dist_logprob_right;
  double score, PPscore;

  left  = node->left;
  right = node->right;

  low = node->height;
  left_low = left->height;
  right_low = right->height;

  for (i = i_min; i < i_max; ++i)
  {
    rel_age_node = (1.0 / opt_grid_intervals) * (low+i);
    abs_age_node = (low+i)*interval_age;

    /* check the ages of the left child */
    if (!left->left)
      jmax = 1;
    else
      jmax = (node->height + i - left->height);

    assert(jmax <= left->entries);

    long jbest = -1;
    double jbest_score = -__DBL_MAX__;
    double jsum_score = 0;
    for (j = 0; j < jmax; ++j)
    {
      assert(j+left->height < node->height + i);
      rel_age_left = (1.0 / opt_grid_intervals) * (left_low+j);
      abs_age_left = (left_low+j)*interval_age;
      age_diff = rel_age_node - rel_age_left;
      prob_rate_left = gamma_dist_logpdf(left->length / age_diff);
      dist_logprob_left = check_prior(left, abs_age_left);
      score = left->matrix[j] + prob_rate_left + dist_logprob_left;
      PPscore = left->matrix_PP[j] + prob_rate_left + dist_logprob_left; 
      jsum_score = calc_log_sum(jsum_score, PPscore);
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

    long kbest = -1;
    double kbest_score = -__DBL_MAX__;
    double ksum_score = 0;
    for (k = 0; k < kmax; ++k)
    {
      assert(k+right->height < node->height + i);
      rel_age_right = (1.0 / opt_grid_intervals) * (right_low+k);
      abs_age_right = (right_low+k)*interval_age;
      age_diff = rel_age_node - rel_age_right;
      prob_rate_right = gamma_dist_logpdf(right->length / age_diff);
      dist_logprob_right = check_prior(right, abs_age_right);
      score = right->matrix[k] + prob_rate_right + dist_logprob_right;
      PPscore = right->matrix_PP[k] + prob_rate_right + dist_logprob_right;
      ksum_score = calc_log_sum(ksum_score, PPscore); 
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


    score = bd_term + jbest_score + kbest_score;/*TODO: the prior is now only added at the children and when you hit the root. Double check this is correct!*/
    PPscore = bd_term + jsum_score + ksum_score;
        /* if it's the root add one more term */
    if (node->height == tree_height)
    {
      if (opt_method_relative || opt_method_nodeprior)
        {
        score += bd_relative_root(node->leaves,
                                  rel_age_node);
        PPscore += bd_relative_root(node->leaves,
                                  rel_age_node); 
        }
      else if (opt_method_tipdates)
        {
        score += bd_tipdates_root(node->leaves,
                                   rel_age_node);
        PPscore += bd_tipdates_root(node->leaves,
                                  rel_age_node);
       }
      else assert(0);
      dist_logprob = check_prior(node, abs_age_node);
      score += dist_logprob;
      PPscore += dist_logprob;
    }

    /* store best placement of children and likelihood for interval line i */
    node->matrix[i] = score;
    node->matrix_PP[i] = PPscore;
    node->matrix_left[i] = jbest;
    node->matrix_right[i] = kbest;
  }

}


static thread_info_t * ti;

static inline void dp_recurse_worker(long t)
{

  thread_info_t * tip = ti + t;
  tree_node_t * node  = tip->node;

  calc_interval_scores(node, tip->entry_first, tip->entry_first + tip->entry_count);
  
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
    calc_tip(node);
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
      node->matrix_left  = (long *)xmalloc(sizeof(long));
      node->matrix_right = (long *)xmalloc(sizeof(long));
      node->matrix_PP       = (double *)xmalloc(sizeof(double));
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
      node->matrix_left  = (long *)xmalloc((size_t)entries * sizeof(long));
      node->matrix_right = (long *)xmalloc((size_t)entries * sizeof(long));
      node->matrix_PP       = (double *)xmalloc((size_t)entries * sizeof(double));
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
  node->matrix_left  = (long *)xmalloc((size_t)entries * sizeof(long));
  node->matrix_right = (long *)xmalloc((size_t)entries * sizeof(long));
  node->matrix_PP       = (double *)xmalloc((size_t)entries * sizeof(double));
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
     calc_tip(node);
     return;
  }

  /* inner nodes case */
  left  = node->left;
  right = node->right;

  dp_recurse_serial(left);
  dp_recurse_serial(right);

  /* run DP */

  /*

                        o  (rel_age_node)
                       / \
                      /   \
      (rel_age_left) o     \
                            o (rel_age_right)

  */
  calc_interval_scores(node,0,node->entries);

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

  assert(score < 0);

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
    printf("\n\n*** RATES OPTIMIZATION ***\n\n");
    printf("Starting score: %f\n\n", score);

    opt_quiet = 1;

    printf("Starting parameters:\n");
    if (opt_parameters_bitv & PARAM_LAMBDA)
      printf(" lambda: %6.4f\n", opt_lambda);
    if (opt_parameters_bitv & PARAM_MU)
      printf(" mu:     %6.4f\n", opt_mu);
    if (opt_parameters_bitv & PARAM_RHO)
      printf(" rho:    %6.4f\n", opt_rho);
    if (opt_parameters_bitv & PARAM_PSI)
      printf(" psi:    %6.4f\n", opt_psi);
    printf("\n");

    /* single param iterative optimization */
    if (opt_threads > 1)
      threads_init();
    double cur_score = score + 1;
    int i = 0;
    while (fabs(cur_score - score) > opt_epsilon)
    {
      printf("[%d]\n", i++);
      cur_score = score;
      //opt_parameters(tree, PARAM_PSI, opt_factor, opt_pgtol);
      if (opt_parameters_bitv & PARAM_LAMBDA)
      {
        score = opt_parameters(tree, PARAM_LAMBDA, opt_factor, opt_pgtol);
        printf("%15.4f lambda: %6.4f\n", score, opt_lambda);
      }
      if (opt_parameters_bitv & PARAM_MU)
      {
        score = opt_parameters(tree, PARAM_MU, opt_factor, opt_pgtol);
        printf("%15.4f mu:     %6.4f\n", score, opt_mu);
      }
      if (opt_parameters_bitv & PARAM_PSI)
      {
        score = opt_parameters(tree, PARAM_PSI, opt_factor, opt_pgtol);
        printf("%15.4f psi:    %6.4f\n", score, opt_psi);
      }
      if (opt_parameters_bitv & PARAM_RHO)
      {
        score = opt_parameters(tree, PARAM_RHO, opt_factor, opt_pgtol);
        printf("%15.4f rho:    %6.4f\n", score, opt_rho);
      }
    }
    if (opt_threads > 1)
      threads_exit();

    printf("\nFinal parameters:\n");
    if (opt_parameters_bitv & PARAM_LAMBDA)
      printf(" lambda: %6.4f\n", opt_lambda);
    if (opt_parameters_bitv & PARAM_MU)
      printf(" mu:     %6.4f\n", opt_mu);
    if (opt_parameters_bitv & PARAM_RHO)
      printf(" rho:    %6.4f\n", opt_rho);
    if (opt_parameters_bitv & PARAM_PSI)
      printf(" psi:    %6.4f\n", opt_psi);

    printf("\nScore after optimization %f\n", score);
    printf("\n*** ***** ************ ***\n\n");
  }
}
