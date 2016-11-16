/* mw_step.c -- Execute one timestep of the E-M wave simulator

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

/* Move the E and B fields forward one timestep. */
int
mw_step(mwDomain *domain)
{
  real dt_dx = 0.5*domain->dt/domain->dx;
  real Eprefix_vacuum = 0.5*domain->dt*domain->c*domain->c
    / domain->dx;
  int i, j;

  /* If this is the first call then create a convenience field that
     reduces the number of multiplications and divisions. */
  if (domain->Eprefix == NULL) {
    mw_new_field(&domain->Eprefix, domain->nx, domain->ny, 1.0);
    for (j = 0; j < domain->ny-1; j++) {
      for (i = 0; i < domain->nx-1; i++) {
	domain->Eprefix[j][i] = 0.5*domain->dt*domain->c*domain->c
	  /(domain->dx*domain->epsilon[j][i]);
	domain->Edamping[j][i] = domain->Bdamping[j][i]
	  * exp(-2.0*M_PI*domain->primary_frequency*domain->dt
		*domain->Edamping[j][i]/domain->epsilon[j][i]);
      }
    }
  }

  /* If wave has a horizontally polarized component... */
  if (domain->mode & MW_MODE_EXY) {
    /* Increment the Ex and Ey components. */
    for (j = 0; j < domain->ny-1; j++) {
      for (i = 0; i < domain->nx-1; i++) {
	domain->Ex[j][i] = domain->Edamping[j][i]*domain->Ex[j][i]
	  + domain->dt*(domain->forcingI[j][i]*domain->Ex_forcingI
			-domain->forcingQ[j][i]*domain->Ex_forcingQ)
	  + domain->Eprefix[j][i]*(domain->Bz[j+1][i+1] - domain->Bz[j][i+1]);
	domain->Ey[j][i] = domain->Edamping[j][i]*domain->Ey[j][i]
	  + domain->dt*(domain->forcingI[j][i]*domain->Ey_forcingI
			-domain->forcingQ[j][i]*domain->Ey_forcingQ)
	  + domain->Eprefix[j][i]*(domain->Bz[j+1][i] - domain->Bz[j+1][i+1]);
      }
    }
  }
  /* If wave has a vertically polarized component... */
  if (domain->mode & MW_MODE_EZ) {
    /* Increment the Ez component. */
    for (j = 1; j < domain->ny-1; j++) {
      for (i = 1; i < domain->nx-1; i++) {
	domain->Ez[j][i] = domain->Edamping[j][i]*domain->Ez[j][i]
	  + domain->dt*(domain->forcingI[j][i]*domain->Ez_forcingI
			-domain->forcingQ[j][i]*domain->Ez_forcingQ)
	  + domain->Eprefix[j][i]*(domain->By[j-1][i] - domain->By[j-1][i-1]
			   - domain->Bx[j][i-1] + domain->Bx[j-1][i-1]);
      }
    }
    /* Increment the Bx and By components. */
    for (j = 0; j < domain->ny-1; j++) {
      for (i = 0; i < domain->nx-1; i++) {
	domain->Bx[j][i] = domain->Bdamping[j][i]*domain->Bx[j][i]
	  - dt_dx*(domain->Ez[j+1][i+1] - domain->Ez[j][i+1]);
	domain->By[j][i] = domain->Bdamping[j][i]*domain->By[j][i]
	  - dt_dx*(domain->Ez[j+1][i] - domain->Ez[j+1][i+1]);
      }
    }
  }
  /* If wave has a horizontally polarized component. */
  if (domain->mode & MW_MODE_EXY) {
    /* Increment the Bz component. */
    for (j = 1; j < domain->ny-1; j++) {
      for (i = 1; i < domain->nx-1; i++) {
	domain->Bz[j][i] = domain->Bdamping[j][i]*domain->Bz[j][i]
	  - dt_dx*(domain->Ey[j-1][i] - domain->Ey[j-1][i-1]
		   - domain->Ex[j][i-1] + domain->Ex[j-1][i-1]);
      }
    }
  }

  /* Is a parallel calculation required for vacuum? */
  if (domain->mode & MW_MODE_VACUUM) {
    if (domain->mode & MW_MODE_EXY) {
      /* Increment the Ex and Ey components. */
      for (j = 0; j < domain->ny-1; j++) {
	for (i = 0; i < domain->nx-1; i++) {
	  domain->Ex_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->Ex_vacuum[j][i]
	    + domain->dt*(domain->forcingI[j][i]*domain->Ex_forcingI
			  -domain->forcingQ[j][i]*domain->Ex_forcingQ)
	    + Eprefix_vacuum
	    *(domain->Bz_vacuum[j+1][i+1] - domain->Bz_vacuum[j][i+1]);
	  domain->Ey_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->Ey_vacuum[j][i]
	    + domain->dt*(domain->forcingI[j][i]*domain->Ey_forcingI
			  -domain->forcingQ[j][i]*domain->Ey_forcingQ)
	    + Eprefix_vacuum
	    *(domain->Bz_vacuum[j+1][i] - domain->Bz_vacuum[j+1][i+1]);
	}
      }
    }
    /* If wave has a vertically polarized component... */
    if (domain->mode & MW_MODE_EZ) {
      /* Increment the Ez component. */
      for (j = 1; j < domain->ny-1; j++) {
	for (i = 1; i < domain->nx-1; i++) {
	  domain->Ez_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->Ez_vacuum[j][i]
	    + domain->dt*(domain->forcingI[j][i]*domain->Ez_forcingI
			  -domain->forcingQ[j][i]*domain->Ez_forcingQ)
	    + Eprefix_vacuum
	    *(domain->By_vacuum[j-1][i] - domain->By_vacuum[j-1][i-1]
	      - domain->Bx_vacuum[j][i-1] + domain->Bx_vacuum[j-1][i-1]);
	}
      }
      /* Increment the Bx and By components. */
      for (j = 0; j < domain->ny-1; j++) {
	for (i = 0; i < domain->nx-1; i++) {
	  domain->Bx_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->Bx_vacuum[j][i]
	    - dt_dx*(domain->Ez_vacuum[j+1][i+1] - domain->Ez_vacuum[j][i+1]);
	  domain->By_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->By_vacuum[j][i]
	    - dt_dx*(domain->Ez_vacuum[j+1][i] - domain->Ez_vacuum[j+1][i+1]);
	}
      }
    }
    /* If wave has a horizontally polarized component. */
    if (domain->mode & MW_MODE_EXY) {
      /* Increment the Bz component. */
      for (j = 1; j < domain->ny-1; j++) {
	for (i = 1; i < domain->nx-1; i++) {
	  domain->Bz_vacuum[j][i]
	    = domain->Bdamping[j][i]*domain->Bz_vacuum[j][i]
	    - dt_dx*(domain->Ey_vacuum[j-1][i] - domain->Ey_vacuum[j-1][i-1]
	     - domain->Ex_vacuum[j][i-1] + domain->Ex_vacuum[j-1][i-1]);
	}
      }
    }
  }
  domain->time += domain->dt;
  return MW_SUCCESS;
}
