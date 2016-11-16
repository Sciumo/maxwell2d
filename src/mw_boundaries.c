/* mw_boundaries.c -- Find the boundaries of a dielectric constant field

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

#include <stdio.h>
#include "maxwell.h"

/* Return the maximum of the nine numbers entered, but incremented by
   1.0 if they are all the same */
real
get_max(real r1, real r2, real r3, real r4, real r5,
	real r6, real r7, real r8, real r9)
{
  real rmax = r1;
  if (r1 == r2 && r1 == r3 && r1 == r4 && r1 == r5 
      && r1 == r6 && r1 == r7 && r1 == r8 && r1 == r9) {
    return r1+1.0;
  }
  else {
    if (r2 > rmax) rmax = r2;
    if (r3 > rmax) rmax = r3;
    if (r4 > rmax) rmax = r4;
    if (r5 > rmax) rmax = r5;
    if (r6 > rmax) rmax = r6;
    if (r7 > rmax) rmax = r7;
    if (r8 > rmax) rmax = r8;
    if (r9 > rmax) rmax = r9;
  }
  return rmax;
}

/* Set the domain->boundaries matrix to contain ones and zeros
   demarking the edges of the regions of different dielectric
   constant */
int
mw_find_boundaries(mwDomain *domain)
{
  int i, j;
  /* Assume the domain is already set to zeros */

  /* Decide if each pixel should be set as a boundary */
  for (j = 1; j < domain->ny-2; j++) {
    for (i = 1; i < domain->nx-2; i++) {
      real emax = get_max(domain->epsilon[j-1][i-1],
			  domain->epsilon[j-1][i],
			  domain->epsilon[j-1][i+1],
			  domain->epsilon[j][i-1],
			  domain->epsilon[j][i],
			  domain->epsilon[j][i+1],
			  domain->epsilon[j+1][i-1],
			  domain->epsilon[j+1][i],
			  domain->epsilon[j+1][i+1]);
      if (domain->epsilon[j][i] == emax) {
	domain->boundaries[j][i] = 1.0;
      }
    }
  }

  return MW_SUCCESS;
}
