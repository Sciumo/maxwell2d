/* mw_frame.c -- Simulate a "frame" of the simulation (several timesteps)

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

/* Run a "frame" of the simulation (usually 7 timesteps) and calculate
   the Poynting vector summation */
int
mw_frame(mwDomain *domain)
{
  int i, j, l;
  /* The Poynting vector calculation is multiplied by 0.5 because of
     the need to use B-field values at two points on the grid, and is
     divided by the magnetic constant to yield the correct units */
  real POYNTING_FACTOR = 0.5 / (4.0 * M_PI * 1.0e-7);

  /* Run simulation forward several timesteps, updating the "forcing"
     each time */
  for (l = 0; l < MW_MINOR_STEPS; l++) {
    if (domain->time*domain->primary_frequency < domain->cycles) {
      real oscillatorI = 0.0, oscillatorQ = 0.0;
      if (domain->nfrequencies == 0) {
	/* Only a primary_frequency has been assigned */
	oscillatorI = sin(domain->time*domain->primary_frequency*2.0*M_PI);
	oscillatorQ = cos(domain->time*domain->primary_frequency*2.0*M_PI);
      }
      else {
	/* The "frequencies" vector is treated in groups of three,
	   with the first element being the frequency, the second the
	   amplitude scaling and the third the phase offset in
	   degrees */
	int ifreq;
	for (ifreq = 0; ifreq < domain->nfrequencies; ifreq++) {
	  oscillatorI += domain->frequencies[ifreq*3+1]
	    *sin(2.0*M_PI*(domain->time*domain->frequencies[ifreq*3]
			   +domain->frequencies[ifreq*3+2]));
	  oscillatorQ += domain->frequencies[ifreq*3+1]
	    *cos(2.0*M_PI*(domain->time*domain->frequencies[ifreq*3]
			   +domain->frequencies[ifreq*3+2]));
	}
      }
      domain->Ex_forcingI = domain->Ex_amplitude*oscillatorI;
      domain->Ey_forcingI = domain->Ey_amplitude*oscillatorI;
      domain->Ez_forcingI = domain->Ez_amplitude*oscillatorI;
      domain->Ex_forcingQ = domain->Ex_amplitude*oscillatorQ;
      domain->Ey_forcingQ = domain->Ey_amplitude*oscillatorQ;
      domain->Ez_forcingQ = domain->Ez_amplitude*oscillatorQ;
    }
    else {
      domain->Ex_forcingI = 0.0;
      domain->Ey_forcingI = 0.0;
      domain->Ez_forcingI = 0.0;
      domain->Ex_forcingQ = 0.0;
      domain->Ey_forcingQ = 0.0;
      domain->Ez_forcingQ = 0.0;
    }
    MW_CHECK(mw_step(domain));
  }

  /* Calculate the Poynting vector */
  if (domain->mode & MW_MODE_EZ) {
    for (j = 1; j < domain->ny-1; j++) {
      for (i = 1; i < domain->nx-1; i++) {
	domain->Poynting_x[j][i] -= POYNTING_FACTOR*domain->Ez[j][i]
	  *(domain->By[j-1][i-1]+domain->By[j-1][i]);
	domain->Poynting_y[j][i] += POYNTING_FACTOR*domain->Ez[j][i]
	  *(domain->Bx[j-1][i-1]+domain->Bx[j][i-1]);
      }
    }
    if (domain->mode & MW_MODE_VACUUM) {
      for (j = 1; j < domain->ny-1; j++) {
	for (i = 1; i < domain->nx-1; i++) {
	  domain->Poynting_x_scat[j][i] -= POYNTING_FACTOR
	    *(domain->Ez[j][i]-domain->Ez_vacuum[j][i])
	    *(domain->By[j-1][i-1]-domain->By_vacuum[j-1][i-1]
	      +domain->By[j-1][i]-domain->By_vacuum[j-1][i]);
	  domain->Poynting_y_scat[j][i] += POYNTING_FACTOR
	    *(domain->Ez[j][i]-domain->Ez_vacuum[j][i])
	    *(domain->Bx[j-1][i-1]-domain->Bx_vacuum[j-1][i-1]
	      +domain->Bx[j][i-1]-domain->Bx_vacuum[j][i-1]);
	}
      }
    }
  }

  if (domain->mode & MW_MODE_EXY) {
    for (j = 0; j < domain->ny-1; j++) {
      for (i = 0; i < domain->nx-1; i++) {
	domain->Poynting_x[j][i] += POYNTING_FACTOR*domain->Ey[j][i]
	  *(domain->Bz[j+1][i]+domain->By[j+1][i+1]);
	domain->Poynting_y[j][i] -= POYNTING_FACTOR*domain->Ex[j][i]
	  *(domain->Bz[j][i+1]+domain->Bz[j+1][i+1]);
      }
    }
    if (domain->mode & MW_MODE_VACUUM) {
      for (j = 1; j < domain->ny-1; j++) {
	for (i = 1; i < domain->nx-1; i++) {
	  domain->Poynting_x_scat[j][i] += POYNTING_FACTOR
	    *(domain->Ey[j][i]-domain->Ey_vacuum[j][i])
	    *(domain->Bz[j+1][i]-domain->Bz_vacuum[j+1][i]
	      +domain->Bz[j+1][i+1]-domain->Bz_vacuum[j+1][i+1]);
	  domain->Poynting_y_scat[j][i] -= POYNTING_FACTOR
	    *(domain->Ex[j][i]-domain->Ex_vacuum[j][i])
	    *(domain->Bz[j][i+1]-domain->Bz_vacuum[j][i+1]
	      +domain->Bz[j+1][i+1]-domain->Bz_vacuum[j+1][i+1]);
	}
      }
    }
  }

  domain->iframe++;
  return MW_SUCCESS;
}
