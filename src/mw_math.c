/* mw_math.c -- Simple mathematical functions applied to entire fields

   Copyright (C) 2008 Robin Hogan <r.j.hogan@reading.ac.uk> 

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>
#include "maxwell.h"

/* Subtract arg2 from arg1 and return the answer in ans */
int
mw_subtract(int nx, int ny, real **arg1, real **arg2, real **ans)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      ans[j][i] = arg1[j][i]-arg2[j][i];
    }
  }
  return MW_SUCCESS;
}

/* Scale the field "arg" by "factor" */
int
mw_scale(int nx, int ny, real **arg, real factor)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      arg[j][i] *= factor;
    }
  }
  return MW_SUCCESS;
}

/* Calculate the complex susceptibility from the complex refractive
   index using the Clausius-Mosotti relationship
   xi=3(n^2-1)/(n^2+2). The input ni may be positive or negative, but
   xii is always returned as positive. */
/*
int
mw_susceptibility(real nr, real ni, real *xir, real *xii)
{
  real er = nr*nr - ni*ni;
  real ei = 2.0*nr*ni;
  real denom = er*er + ei*ei + 4.0*er + 4.0;
  *xir = (er*er + ei*ei + er - 2.0) * 3.0 / denom;
  *xii = fabs(ei) * 9.0 / denom;
  return MW_SUCCESS;
}
*/

/* Return the complex susceptibility *not* using Clausius-Mosotti */
int
mw_susceptibility(real nr, real ni, real *xir, real *xii)
{
  *xir = nr-1.0;
  *xii = ni;
  return MW_SUCCESS;
}
