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

static double diff;
static double c1;
static double c2;
static long k,m;

void bd_init(long fossils_count, long extinct_leaves_count)
{
  diff = opt_mu - opt_lambda;
  c1 = sqrt( pow((opt_lambda - opt_mu - opt_psi),2) + 4*opt_lambda*opt_psi);
  c2 = - (opt_lambda - opt_mu - 2*opt_lambda*opt_rho - opt_psi) / c1;

  k = fossils_count;
  m = extinct_leaves_count;
}

/* Computes the first part (before the products) of the special case of
 * Equation (9) from (Stadler 2010) for which no fossil information is
 * available as part of the tree */
double bd_relative_root(long leaves, double t)
{
  double terma = leaves * (opt_lambda - opt_mu) * exp(diff * t);
  double termb = opt_rho*opt_lambda + (opt_lambda*(1 - opt_rho) - opt_mu) *
                                                                exp(diff*t);

  return log(terma / termb);
}

/* Computes the second part (products) of the special case of Equation (9)
 * from (Stadler 2010) for which no fossil information is available */
double bd_relative_prod(double t)
{
  double terma = (opt_lambda * opt_rho * diff * diff * exp(diff*t));
  double termb = (opt_rho*opt_lambda + (opt_lambda*(1 - opt_rho) - opt_mu) * 
                                                                exp(diff*t));

  return log(terma / (termb * termb));
}

double bd_tipdates_prod_inner(double t) /*pi(t) from Stadler 2010 eq.2*/
{
  double term;

  term = (2*(1 - c2*c2) + exp(-c1*t)*(1-c2)*(1-c2) + exp(c1*t)*(1+c2)*(1+c2));

  
  double p1 = 4*opt_rho / term;

  return log(opt_lambda * p1);
}

double bd_tipdates_prod_tip(double t)
{
  double term;

  term = (opt_lambda + opt_mu + opt_psi) + c1* (exp(-c1*t)*(1-c2)-(1+c2)) /
                                               (exp(-c1*t)*(1-c2)+(1+c2));
  
  double p0 = term / (2*opt_lambda);

  term = (2*(1 - c2*c2) + exp(-c1*t)*(1-c2)*(1-c2) + exp(c1*t)*(1+c2)*(1+c2));
  double p1 = 4*opt_rho / term;
  
  return log(p0/p1);
}

double bd_tipdates_root(long leaves, double t)
{
  double terma = log(4 * leaves * opt_rho) + (k+m)*log(opt_psi);
  double termb = log(c1*(c2+1)*(1 - c2 + (1+c2)*exp(c1*t)));

  return terma - termb;
}
