/* mw_print.c -- Functions to print 2D files as ASCII to a file

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

int
mw_print_field(FILE *file, real **field, int nx, int ny)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      fprintf(file, " %6g", field[j][i]);
    }
    fprintf(file, "\n");
  }
  return MW_SUCCESS;
}

int
mw_visualize_field(FILE *file, real **field, int nx, int ny)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      if (field[j][i] > 0.0) {
	fprintf(file, "+");
      }
      else if (field[j][i] < -0.0) {
	fprintf(file, "-");
      }
      else {
	fprintf(file, " ");
      }
    }
    fprintf(file, "\n");
  }
  return MW_SUCCESS;
}
