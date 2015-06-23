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

static double diff;

void bd_init(void)
{
  diff = opt_mu - opt_lambda;
}

/* Computes the first part (before the products) of the special case of
 * Equation (9) from (Stadler 2010) for which no fossil information is
 * available */
double bd_nofossil_root(int leaves, double t)
{
  double terma = leaves * (opt_lambda - opt_mu) * exp(diff * t);
  double termb = opt_rho*opt_lambda + (opt_lambda*(1 - opt_rho) - opt_mu) *
                                                                exp(diff*t);

  return log(terma / termb);
}

/* Computes the second part (products) of the special case of Equation (9)
 * from (Stadler 2010) for which no fossil information is available */
double bd_nofossil_prod(double t)
{
  double terma = (opt_lambda * opt_rho * diff * diff * exp(diff*t));
  double termb = (opt_rho*opt_lambda + (opt_lambda*(1 - opt_rho) - opt_mu) * 
                                                                exp(diff*t));

  return log(terma / (termb * termb));
}
