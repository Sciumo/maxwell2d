/* mw_alloc.c -- Memory allocation functions

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

#include <stdlib.h>
#include <math.h>

#include "maxwell.h"

/* Initialize a matrix of real numbers with a specified size and set
   every element to "value" */
int
mw_new_field(real ***field, int nx, int ny, int value)
{
  int i;
  *field = (real**) malloc(sizeof(real*)*ny);
  if (!*field) {
    return MW_FAILURE;
  }
  **field = (real*) malloc(sizeof(real)*ny*nx);
  if (!**field) {
    return MW_FAILURE;
  }
  for (i = 1; i < ny; i++) {
    (*field)[i] = (*field)[0] + i*nx;
  }
  mw_reset_field(*field, nx, ny, value);
  return MW_SUCCESS;
}

/* Initialize a new domain with a specified size, and use "mode" to
   indicate which fields are required */
int
mw_new_domain(mwDomain *domain, int nx, int ny, real dx, int mode)
{
  if (mode & MW_MODE_EXY) {
    mw_new_field(&domain->Ex, nx, ny, 0.0);
    mw_new_field(&domain->Ey, nx, ny, 0.0);
    mw_new_field(&domain->Bz, nx, ny, 0.0);
  }
  if (mode & MW_MODE_EZ) {
    mw_new_field(&domain->Ez, nx, ny, 0.0);
    mw_new_field(&domain->Bx, nx, ny, 0.0);
    mw_new_field(&domain->By, nx, ny, 0.0);
  }
  if (mode & MW_MODE_VACUUM) {
    if (mode & MW_MODE_EXY) {
      mw_new_field(&domain->Ex_vacuum, nx, ny, 0.0);
      mw_new_field(&domain->Ey_vacuum, nx, ny, 0.0);
      mw_new_field(&domain->Bz_vacuum, nx, ny, 0.0);
    }
    if (mode & MW_MODE_EZ) {
      mw_new_field(&domain->Ez_vacuum, nx, ny, 0.0);
      mw_new_field(&domain->Bx_vacuum, nx, ny, 0.0);
      mw_new_field(&domain->By_vacuum, nx, ny, 0.0);
    }
  }

  mw_new_field(&domain->epsilon, nx, ny, 1.0);
  mw_new_field(&domain->Edamping, nx, ny, 0.0);
  mw_new_field(&domain->Bdamping, nx, ny, 1.0);
  mw_new_field(&domain->forcingI, nx, ny, 0.0);
  mw_new_field(&domain->forcingQ, nx, ny, 0.0);
  mw_new_field(&domain->boundaries, nx, ny, 0.0);
  mw_new_field(&domain->Poynting_x, nx, ny, 0.0);
  mw_new_field(&domain->Poynting_y, nx, ny, 0.0);
  if (mode & MW_MODE_VACUUM) {
    mw_new_field(&domain->Poynting_x_scat, nx, ny, 0.0);
    mw_new_field(&domain->Poynting_y_scat, nx, ny, 0.0);
  }
  domain->Eprefix = NULL;

  domain->mode = mode;
  domain->nx = nx;
  domain->ny = ny;
  domain->dx = dx;
  domain->c = MW_C;
  domain->dt = 0.8 * dx / domain->c;
  domain->dt_dx = domain->dt/domain->dx;
  domain->Ez_forcingI = 1.0;
  domain->Ez_forcingQ = 0.0;
  domain->Ex_forcingI = domain->Ey_forcingI = 0.0;
  domain->Ex_forcingQ = domain->Ey_forcingQ = 0.0;
  domain->time = 0.0;
  domain->iframe = 0;
  domain->plot_E_max = 1.0e-8;
  domain->plot_B_max = domain->plot_E_max/domain->c;
  domain->plot_scat_ratio = 1.0;
  domain->mag = 1;
  domain->epsilon_plot_file = NULL;
  domain->frequencies = NULL;
  domain->nfrequencies = 0;
  return MW_SUCCESS;
}

/* Free the memory used to store a matrix field */
int
mw_free_field(real **field)
{
  if (field) {
    if (*field) {
      free(*field);
    }
    free(field);
  }
  return MW_SUCCESS;
}

/* Free the memory used to store an entire domain */
int
mw_free_domain(mwDomain *domain)
{
  mw_free_field(domain->Ex);
  mw_free_field(domain->Ey);
  mw_free_field(domain->Ez);
  mw_free_field(domain->Bx);
  mw_free_field(domain->By);
  mw_free_field(domain->Bz);
  mw_free_field(domain->epsilon);
  mw_free_field(domain->Edamping);
  mw_free_field(domain->Bdamping);
  mw_free_field(domain->Eprefix);
  domain->Ex = domain->Ey = domain->Ez = NULL;
  domain->Bx = domain->By = domain->Bz = NULL;
  domain->epsilon = domain->Edamping = domain->Bdamping
    = domain->Eprefix = NULL;
  return MW_SUCCESS;
}

/* Set all the elements of an existing field to a particular value */
int
mw_reset_field(real **field, int nx, int ny, real value)
{
  int i, j;
  for (j = 0; j < ny; j++) {
    for (i = 0; i < nx; i++) {
      field[j][i] = value;
    }
  }
  return MW_SUCCESS;
}
