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


#ifndef M_PI
#define 3.14159265358979323846
#endif

double norm_dist_pdf(double mean, double stdev, double x)
{
  assert (stdev > 0);

  double terma = exp(-((x-mean)*(x-mean)/(2*stdev*stdev)));
  double termb = stdev*sqrt(2*M_PI);

  return (terma/termb);
}

double norm_dist_logpdf(double mean, double stdev, double x)
{
  assert (stdev> 0);

  double terma = -(x-mean)*(x-mean) / (2*stdev*stdev);
  double termb = log(sqrt(2*M_PI)*stdev);

  return terma-termb;
}
