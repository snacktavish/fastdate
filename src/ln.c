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

#ifndef M_SQRTPI
#define M_SQRTPI       1.77245385090551602729  /* square root of pi*/
#endif

#ifndef M_SQRT_2
#define M_SQRT_2       1.41421356237309504880 /* square root of 2 */
#endif




static double constant = M_SQRTPI * M_SQRT_2;

double ln_dist_pdf(double mean, double stdev, double x)
{
  assert(stdev > 0);
  assert(x > 0);

  return (1 / (x * stdev * constant)) *
         exp(-(log(x) - mean)*(log(x) - mean) / (2 * stdev * stdev));
}

double ln_dist_logpdf(double mean, double stdev, double x)
{
  assert(stdev > 0);
  assert(x >= 0);

  if (x == 0) return  -__DBL_MAX__ / 2;

  return (-log(x) - log(stdev) - log(constant) -
         (log(x) - mean)*(log(x) - mean) / (2 * stdev * stdev));

}


