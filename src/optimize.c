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

static int count_bits(int i)
{
     i = i - ((i >> 1) & 0x55555555);
     i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
     return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

static double compute_score( optimize_options_t * params,
                             double * x
                           )
{
  double next_score = 0.0;

  /* TODO: implement */

  return next_score;
}

static void compute_gradients( optimize_options_t * params,
                               double * x,
                               double * g,
                               int n
                             )
{
  int i;
  double h, temp;
  double local_score;

  local_score = compute_score (params, x);

  for (i = 0; i < n; i++)
  {
    double ERROR_X = 1.0e-4;
    temp = x[i];
    h = ERROR_X * fabs (temp);
    if (h < 1e-12)
      h = ERROR_X;

    x[i] = temp + h;
    h = x[i] - temp;
    double lnderiv = compute_score (params, x);

    g[i] = (lnderiv - local_score) / h;

    /* reset variable */
    x[i] = temp;
  }
}

double optimize_parameters(optimize_options_t * params)
{
  int i;
  int continue_opt;
  int iprint = -1;

  /* L-BFGS-B */
  int max_corrections;
  int num_variables;
  double score = 0;
  double *x, *g, *lower_bounds, *upper_bounds, *wa;
  int *bound_type, *iwa;
  int taskValue;
  int *task = &taskValue;
  int csaveValue;
  int *csave = &csaveValue;
  double dsave[29];
  int isave[44];
  logical lsave[4];

  num_variables = count_bits(params->which_parameters);
  if (!num_variables)
  {
    /* nothing to optimize */
    return 0;
  }

  /* memory allocation */
  x = (double *) calloc ((size_t) num_variables, sizeof(double));
  g = (double *) calloc ((size_t) num_variables, sizeof(double));
  lower_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  upper_bounds = (double *) calloc ((size_t) num_variables, sizeof(double));
  bound_type = (int *) calloc ((size_t) num_variables, sizeof(int));
  iwa = (int *) calloc (3 * (size_t)num_variables, sizeof(int));
  wa = (double *) calloc (
      (2 * (size_t)max_corrections + 5) * (size_t)num_variables
      + 12 * (size_t)max_corrections * ((size_t)max_corrections + 1),
      sizeof(double));

  /* set variables and bounds */
  i=0;
  if (params->which_parameters & PARAM_LAMBDA)
    {
      x[i] = params->lambda;
      lower_bounds[i] = params->mu;
      bound_type[i] = LBFGSB_BOUND_LOWER;
      i++;
    }
  else if (params->which_parameters & PARAM_MU)
    {
      x[i] = params->mu;
      lower_bounds[i] = 1e-8;
      upper_bounds[i] = params->lambda;
      bound_type[i] = LBFGSB_BOUND_BOTH;
      i++;
    }
  if (params->which_parameters & PARAM_PSI)
    {
      x[i] = params->psi;
      lower_bounds[i] = 0;
      upper_bounds[i] = 1;
      bound_type[i] = LBFGSB_BOUND_BOTH;
      i++;
    }

  /* check all free variables are set */
  assert(i == num_variables);

  /* start the iteration by initializing task */
  *task = (int) START;

  continue_opt = 1;
  while (continue_opt)
  {
    /* this is the call to the L-BFGS-B code */
    setulb (&num_variables, &max_corrections, x, lower_bounds, upper_bounds,
            bound_type, &score, g, &(params->factr), &(params->pgtol), wa, iwa,
            task, &iprint, csave, lsave, isave, dsave);
    if (IS_FG(*task))
    {
      /*
       * the minimization routine has returned to request the
       * function f and gradient g values at the current x.
       * Compute function value f for the sample problem.
       */

       /* reset bounds */
       if ((params->which_parameters & PARAM_LAMBDA) && (params->which_parameters & PARAM_MU))
         {
           params->lambda = x[0];
           params->mu = x[1];
           lower_bounds[0] = params->mu;
           upper_bounds[1] = params->lambda;
         }

       /* TODO: compute score and gradient */
      score = compute_score (params, x);

      compute_gradients(params, x, g, num_variables);
    }
    else if (*task != NEW_X)
    {
      continue_opt = 0;
    }
  }

  free (iwa);
  free (wa);
  free (x);
  free (g);
  free (lower_bounds);
  free (upper_bounds);
  free (bound_type);

  return score;
} /* optimize_parameters */
