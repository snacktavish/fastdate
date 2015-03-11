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

static double ut;

static double birth_rate;
static double death_rate;
static double diff;

void bd_init(double br, double dr)
{
  birth_rate = br;
  death_rate = dr;

  diff = death_rate - birth_rate;

  /* compute ut */

  if (death_rate == birth_rate)
  {
    assert ( death_rate > 0);
    //ut = death_rate / (1 + death_rate);
    ut = log(death_rate) - log(1 + death_rate);
  }
  else if (death_rate == 0)
  {
    assert((1 - exp(-birth_rate)) > 0);
    ut = log(1 - exp(-birth_rate));
  }
  else
  {
    assert(  birth_rate*(1 - exp(diff)) / (birth_rate - (death_rate * exp(diff))) > 0);
    ut = log( birth_rate*(1 - exp(diff)) / (birth_rate - (death_rate * exp(diff))) );
  }
}

double bd_prob(int leaves, double t)
{
  double pt;
  double e;

  assert(t != 0);

  double factor = (leaves-1) * birth_rate;
  assert(factor > 0);
  factor = log(factor);
  factor -= ut;
   
  if (diff == 0)
  {
    //return factor / ((1 + death_rate*t)*(1 + death_rate*t));
    return factor - log((1 + death_rate*t)*(1 + death_rate*t));
  }
  
  e = exp(diff*t);
  pt = -diff / (birth_rate - death_rate * e);
  assert(pt > 0);
  return factor + log(pt * pt * e);
}
