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

double exp_dist_pdf(double lambda, double x)
{
  assert(lambda > 0);

  if (x < 0) return 0;

  return (lambda * exp(-lambda*x));
}

double exp_dist_logpdf(double lambda, double x)
{
  assert(lambda > 0);
  assert(x >= 0);

  if (x < 0) return -__DBL_MAX__ / 2;

  return log(lambda) - lambda*x;
}
