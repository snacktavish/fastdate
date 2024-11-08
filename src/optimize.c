/*
 Copyright (C) 2015 Diego Darriba

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "lbfgsb/lbfgsb.h"
#include "fastdate.h"

/* L-BFGS-B bound type */

#define LBFGSB_BOUND_NONE  0
#define LBFGSB_BOUND_LOWER 1
#define LBFGSB_BOUND_BOTH  2
#define LBFGSB_BOUND_UPPER 3

#define ERROR_X 1.0e-4

/* information for parameter optimization */
typedef struct
{
  /* which parameter to optimize */
  int which_parameters;
  int num_variables;

  tree_node_t * tree;

  /* optimization level */
  double factr;
  double pgtol;
} optimize_options_t;

static int count_bits(int i)
{
  i = i - ((i >> 1) & 0x55555555);
  i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
  return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

//static int n_comp = 0;

static double compute_score(optimize_options_t * params,
                            double * x,
                            long extinct_leaves_count,
                            long fossils_count)
{
  //printf("BDG compute %d\n", n_comp++);
  int i;
  double score;

  /* set parameters */
  /* set variables and bounds */
  i = 0;
  if (params->which_parameters & PARAM_LAMBDA)
  {
    if (opt_mu_scale)
    {
      /* scale mu accordingly */
      opt_mu *= x[i]/opt_lambda;
#ifdef DEBUG
      printf("x = %f (mu=%f, lambda=%f)\n", x[i], opt_mu, opt_lambda);
#endif
    }
#ifdef DEBUG
    else
    {
      printf("lambda=%f\n", x[i]);
    }
#endif
    opt_lambda = x[i];
    i++;
  }
  if (params->which_parameters & PARAM_MU)
  {
    opt_mu = x[i] * opt_lambda;
#ifdef DEBUG
    printf("x = %f (mu=%f)\n", x[i], opt_mu);
#endif
    i++;
  }
  if (params->which_parameters & PARAM_PSI)
  {
    opt_psi = x[i];
    i++;
  }
  if (params->which_parameters & PARAM_RHO)
  {
    opt_rho = x[i];
    i++;
  }
  if (params->which_parameters & PARAM_RATE_MEAN)
  {
    opt_rate_mean = x[i];
    if (opt_fixgamma)
    {
      opt_rate_var = 1/x[i];
#ifdef DEBUG
      printf("rate_mean = rate_var = %f\n", opt_rate_mean);
#endif
    }
#ifdef DEBUG
    else
    {
      printf("rate_mean = %f\n", opt_rate_mean);
    }
#endif
    i++;
  }
  if (params->which_parameters & PARAM_RATE_VAR)
  {
    assert (!opt_fixgamma);
    opt_rate_var = x[i];
#ifdef DEBUG
    printf("rate_var = %f\n", opt_rate_var);
#endif
    i++;
  }

  /* assert that all free variables are set */
  assert(i == params->num_variables);

  assert (opt_lambda > opt_mu);

  /* evaluate proposal */
  score = dp_evaluate(params->tree,extinct_leaves_count,fossils_count) * -1;

#ifdef DEBUG
  printf(" score = %f\n", score);
#endif
  return score;
}

static void compute_gradients(optimize_options_t * params,
                              double * x,
                              double * g,
                              int n,
                              double score,
                              long extinct_leaves_count,
                              long fossils_count)
{
  int i;
  double h, temp;

  for (i = 0; i < n; i++)
  {
    temp = x[i];
    h = ERROR_X * fabs(temp);
    if (h < 1e-12)
      h = ERROR_X;

    x[i] = temp + h;
    h = x[i] - temp;
    double lnderiv = compute_score(params,
                                   x,
                                   extinct_leaves_count,
                                   fossils_count);

    /* compute gradient */
    g[i] = (lnderiv - score) / h;

    /* reset variable */
    x[i] = temp;
  }
}

static double opt_parameters_lbfgsb(tree_node_t * tree,
                                    int which,
                                    double factr,
                                    double pgtol,
                                    double initial_score,
                                    long extinct_leaves_count,
                                    long fossils_count)
{
  int i;
  int continue_opt;
  int iprint = -1;

  /* L-BFGS-B */
  int max_corrections = 10;
  double score = 0;
  double *start_x, *x, *g, *lower_bounds, *upper_bounds, *wa;
  int *bound_type, *iwa;
  int taskValue;
  int *task = &taskValue;
  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  optimize_options_t * params = (optimize_options_t *) calloc (
      1, sizeof(optimize_options_t));

  params->tree = tree;
  params->factr = factr;
  params->pgtol = pgtol;
  params->which_parameters = which;
  params->num_variables = count_bits (params->which_parameters);
  if (!params->num_variables)
  {
    /* nothing to optimize */
    return 0;
  }

  /* memory allocation */
  x = (double *) calloc ((size_t) params->num_variables, sizeof(double));
  start_x = (double *) calloc ((size_t) params->num_variables, sizeof(double));
  g = (double *) calloc ((size_t) params->num_variables, sizeof(double));
  lower_bounds = (double *) calloc ((size_t) params->num_variables,
                                    sizeof(double));
  upper_bounds = (double *) calloc ((size_t) params->num_variables,
                                    sizeof(double));
  bound_type = (int *) calloc ((size_t) params->num_variables, sizeof(int));
  iwa = (int *) calloc (3 * (size_t) params->num_variables, sizeof(int));
  wa = (double *) calloc (
      (2 * (size_t) max_corrections + 5) * (size_t) params->num_variables
          + 12 * (size_t) max_corrections * ((size_t) max_corrections + 1),
      sizeof(double));

  /* set variables and bounds */
  i = 0;
  if (params->which_parameters & PARAM_LAMBDA)
  {
    x[i] = (opt_lambda > ERROR_X) ? opt_lambda : ERROR_X;
    lower_bounds[i] = opt_mu + 2 * ERROR_X;
    bound_type[i] = LBFGSB_BOUND_LOWER;
    i++;
  }
  if (params->which_parameters & PARAM_MU)
  {
    x[i] = opt_mu/opt_lambda;
    if (x[i] < ERROR_X) x[i] = ERROR_X;
    if (x[i] > (1.0 - ERROR_X)) x[i] = (1.0 - ERROR_X);
    lower_bounds[i] = ERROR_X;
    upper_bounds[i] = 1.0 - ERROR_X;
    bound_type[i] = LBFGSB_BOUND_BOTH;
    i++;
  }
  if (params->which_parameters & PARAM_PSI)
  {
    x[i] = (opt_psi > ERROR_X) ? opt_psi : ERROR_X;
    lower_bounds[i] = ERROR_X;
    upper_bounds[i] = MAX_PSI;
    bound_type[i] = LBFGSB_BOUND_BOTH;
    i++;
  }
  if (params->which_parameters & PARAM_RHO)
  {
    x[i] = (opt_rho > ERROR_X) ? opt_rho : ERROR_X;
    lower_bounds[i] = ERROR_X;
    upper_bounds[i] = MAX_RHO;
    bound_type[i] = LBFGSB_BOUND_BOTH;
    i++;
  }
  if (params->which_parameters & PARAM_RATE_MEAN)
  {
    x[i] = (opt_rate_mean > ERROR_X) ? opt_rate_mean : ERROR_X;
    lower_bounds[i] = ERROR_X;
    // upper_bounds[i] = MAX_RATE_MEAN;
    bound_type[i] = LBFGSB_BOUND_LOWER;
    i++;
  }
  if (params->which_parameters & PARAM_RATE_VAR)
  {
    x[i] = (opt_rate_var > ERROR_X) ? opt_rate_var : ERROR_X;
    lower_bounds[i] = ERROR_X;
    bound_type[i] = LBFGSB_BOUND_LOWER;
    i++;
  }

  /* assert that all free variables are set */
  assert(i == params->num_variables);
  memcpy(start_x, x, params->num_variables * sizeof(double));

  /* start the iteration by initializing task */
  *task = (int) START;

  continue_opt = 1;
  while (continue_opt)
  {
    /* this is the call to the L-BFGS-B code */
    setulb (&params->num_variables, &max_corrections, x, lower_bounds,
            upper_bounds, bound_type, &score, g, &(params->factr),
            &(params->pgtol), wa, iwa, task, &iprint, csave, lsave, isave,
            dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

      /* reset bounds */
      if ((params->which_parameters & PARAM_LAMBDA)
          && (params->which_parameters & PARAM_MU))
      {
        opt_lambda = x[0];
        opt_mu = x[1];
        if (opt_lambda < opt_mu)
        {
          x[0] = opt_mu;
          x[1] = opt_lambda;
          opt_lambda = x[0];
          opt_mu = x[1];
        }
        lower_bounds[0] = opt_mu + ERROR_X;
        upper_bounds[1] = opt_lambda - ERROR_X;
      }

      score = compute_score(params,x,extinct_leaves_count,fossils_count);

      compute_gradients(params,
                        x,
                        g,
                        params->num_variables,
                        score,
                        extinct_leaves_count,
                        fossils_count);
    }
    else if (*task != NEW_X)
    {
      continue_opt = 0;
    }
  }

  score = compute_score(params,x,extinct_leaves_count,fossils_count);

  /* check if score has improved */
  if (score + initial_score > 1e-5)
  {
    /* reset starting values */
    memcpy(x, start_x, params->num_variables * sizeof(double));
    score = compute_score(params,x,extinct_leaves_count,fossils_count);
  }

  free (iwa);
  free (wa);
  free (x);
  free (start_x);
  free (g);
  free (lower_bounds);
  free (upper_bounds);
  free (bound_type);

  free (params);

  return -1 * score;
}

double opt_parameters(tree_node_t * tree,
                      int which,
                      double factr,
                      double pgtol,
                      double cur_score,
                      long extinct_leaves_count,
                      long fossils_count)
{
  double score = opt_parameters_lbfgsb(tree,
                                       which,
                                       factr,
                                       pgtol,
                                       cur_score,
                                       extinct_leaves_count,
                                       fossils_count);

  return score;
} /* optimize_parameters */
