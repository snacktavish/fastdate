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

static double alpha = 0;
static double beta = 0;
static double constant = 0;

void gamma_dist_init(double mean, double variance)
{
  beta  = mean/variance;
  alpha = mean*beta;

  constant = alpha*log(beta) - lgamma(alpha);
}

double gamma_dist_logpdf(double x)
{
  assert(x > 0);
  if (x <= 0) return 0;

  return (constant + (alpha - 1)*log(x) - x*beta);
}


