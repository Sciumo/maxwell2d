/* mw_shape.c -- Functions to add various shapes to the dielectric
   constant field

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

/* Reset the magnetic-field damping with an absorbing border */
int
mw_reset_damping(mwDomain *domain, int borderwidth)
{
  int i, k;
  mw_reset_field(domain->Bdamping, domain->nx, domain->ny, 1.0);
  for (k = 0; k < borderwidth; k++) {
    for (i = k; i < domain->nx-k; i++) {
      domain->Bdamping[k][i] = domain->Bdamping[domain->ny-1-k][i]
	//	= (k+1.0)/(borderwidth+1);
	= sqrt((k+1.0)/(borderwidth+1));
    }
    for (i = k; i < domain->ny-k; i++) {
      domain->Bdamping[i][k] = domain->Bdamping[i][domain->nx-1-k]
	//	= (k+1.0)/(borderwidth+1);
	= sqrt((k+1.0)/(borderwidth+1));
    }
  }
  return MW_SUCCESS;
}

/* Add one or more circles as determined by the vector "var" of length
   "ivar", where each group of five elements corresponds to: (0) x of
   centre, (1) y of centre, (2) radius, (3) epsilon_r, and (4)
   epsilon_i */
int
mw_add_circle(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 4) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real radius = var[2]/domain->dx;
    int minx = x0-radius;
    int maxx = x0+radius+1;
    int miny = y0-radius;
    int maxy = y0+radius+1;
    int i, j;
    real radius2 = radius*radius;
    real xir, xii;
    mw_susceptibility(var[3], var[4], &xir, &xii);
    if (minx < 0) {
      minx = 0;
    }
    if (miny < 0) {
      miny = 0;
    }
    if (maxx > domain->nx-1) {
      maxx = domain->nx-1;
    }
    if (maxy > domain->ny-1) {
      maxy = domain->ny-1;
    }
    for (j = miny; j <= maxy; j++) {
      for (i = minx; i <= maxx; i++) {
	if ((x0-i)*(x0-i) + (y0-j)*(y0-j)
	    < radius2) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 5;
    var += 5;
  }
  
  return MW_SUCCESS;
}

/* Add sharp edge as determined by the vector "var" of length "ivar",
   where each group of five elements corresponds to: (0) x of point on
   the edge, (1) y of point on the edge, (2) angle of edge in degrees,
   (3) epsilon_r on one side of edge, and (4) epsilon_i on one side of
   edge */
int
mw_add_edge(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 4) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real angle = var[2]*M_PI/180.0;
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    int i, j;
    real xir, xii;
    mw_susceptibility(var[3], var[4], &xir, &xii);
    //    fprintf(stderr, "%g %g %g %g\n", var[3], var[4], xir, xii);
    for (j = 0; j < domain->ny; j++) {
      for (i = 0; i < domain->nx; i++) {
	if ((i-x0)*sin_angle+(j-y0)*cos_angle > 0.0) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 5;
    var += 5;
  }
  
  return MW_SUCCESS;
}

/* Add a gradient of dielectric constant as determined by the vector
   "var" of length "ivar", where each group of six elements corresponds
   to: (0) x of centre, (1) y of centre, (2) angle in degrees, (3)
   width factor, (4) epsilon_r offset, and (5) epsilon_r scaling */
int
mw_add_gradient(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 5) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real angle = var[2]*M_PI/180.0;
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    real xfactor = domain->dx/var[5];
    real dist;
    int i, j;
    for (j = 0; j < domain->ny; j++) {
      for (i = 0; i < domain->nx; i++) {
	dist = (i-x0)*sin_angle+(j-y0)*cos_angle;
	domain->epsilon[j][i]
	  += var[3] + var[4]/(1.0+exp(-dist*xfactor));
      }
    }
    nvar -= 6;
    var += 6;
  }
  
  return MW_SUCCESS;
}

/* Add one or more ripples as determined by the vector "var" of length
   "ivar", where each group of six elements corresponds to: (0) x of
   centre, (1) y of centre, (2) angle in degrees, (3) epsilon_r scale
   of ripples, (4) wavelength of ripples, (5) decay scale of ripples */
int
mw_add_ripple(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 5) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real angle = var[2]*M_PI/180.0;
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    real xfactor = domain->dx/var[5];
    real wavenumber = 2.0*M_PI*domain->dx/var[4];
    real dist;
    int i, j;
    for (j = 0; j < domain->ny; j++) {
      for (i = 0; i < domain->nx; i++) {
	dist = (i-x0)*sin_angle+(j-y0)*cos_angle;
	domain->epsilon[j][i]
	  += var[3]*sin(wavenumber*dist)
	  *exp(-pow(xfactor*dist, 4.0));
      }
    }
    nvar -= 6;
    var += 6;
  }
  
  return MW_SUCCESS;
}

/* Add one or more rotated rectangles as determined by the vector
   "var" of length "ivar", where each group of seven elements
   corresponds to: (0) x of centre, (1) y of centre, (2) angle in
   degrees, (3) first axis length, (4) second axis length, (5)
   epsilon_r, (6) epsilon_i */
int
mw_add_rotated_rectangle(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 6) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real angle = var[2]*M_PI/180.0;
    real halfwidth1 = 0.5*var[3]/domain->dx;
    real halfwidth2 = 0.5*var[4]/domain->dx;
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    real dist1, dist2;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[5], var[6], &xir, &xii);
    for (j = 0; j < domain->ny; j++) {
      for (i = 0; i < domain->nx; i++) {
	dist1 = (i-x0)*sin_angle+(j-y0)*cos_angle;
	dist2 = (i-x0)*cos_angle-(j-y0)*sin_angle;
	if (fabs(dist1) <= halfwidth1 && fabs(dist2) <= halfwidth2) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 7;
    var += 7;
  }
  
  return MW_SUCCESS;
}

/* Add one or more "wave packets" as determined by the vector
   "var" of length "ivar", where each group of eight elements
   corresponds to: (0) x of centre, (1) y of centre, (2) angle in
   degrees, (3) first axis length, (4) second axis length, (5)
   wavelength along first axis, (6) epsilon_r maximum, (7) 
   epsilon_i maximum */
int
mw_add_wave_packet(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 6) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real angle = var[2]*M_PI/180.0;
    real halfwidth1 = 0.5*var[3]/domain->dx;
    real halfwidth2 = 0.5*var[4]/domain->dx;
    real wavelength = var[5];
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    real dist1, dist2;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[6], var[7], &xir, &xii);
    for (j = 0; j < domain->ny; j++) {
      for (i = 0; i < domain->nx; i++) {
	dist1 = (i-x0)*sin_angle+(j-y0)*cos_angle;
	dist2 = (i-x0)*cos_angle-(j-y0)*sin_angle;
	if (fabs(dist1) <= halfwidth1 && fabs(dist2) <= halfwidth2) {
	  real tmp = sin(M_PI*(halfwidth1-dist1)/wavelength);
	  real amplitude = tmp*tmp;
	  domain->epsilon[j][i] += xir*amplitude;
	  domain->Edamping[j][i] += xii*amplitude;
	}
      }
    }
    nvar -= 7;
    var += 7;
  }
  
  return MW_SUCCESS;
}

/* Add one or more dish antennae as determined by the vector "var" of
   length "ivar", where each group of seven elements correspond to:
   (0) x of focus, (1) y of focus, (2) distance scale, (3) radius left,
   (4) radius right, (5) dish thickness, (6) epsilon_r, and (7) epsilon_i */
int
mw_add_dish(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 7) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real dist = var[2]/domain->dx;
    real radius1 = var[3]/domain->dx;
    real radius2 = var[4]/domain->dx;
    real thickness = var[5]/domain->dx;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[6], var[7], &xir, &xii);
    for (i = x0-radius1; i <= x0+radius2; i++) {
      int k;
      if (i < 0 || i >= domain->nx) {
	continue;
      }
      k = y0 + (0.25*(i-x0)*(i-x0)/dist - dist);
      for (j = k; j > k-thickness; j--) {
	if (j >= 0 && j < domain->ny) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 7;
    var += 7;
  }
  return MW_SUCCESS;
}

/* Add one or more rectangles as determined by the vector "var" of
   length "ivar", where each group of seven elements corresponds to:
   (0) x of bottom left, (1) y of bottom left, (2) x of top right, (3)
   y of top right, (4) epsilon_r, (5) epsilon_i */
int
mw_add_rectangle(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 5) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real x1 = var[2]/domain->dx + domain->nx/2.0;
    real y1 = var[3]/domain->dx + domain->ny/2.0;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[4], var[5], &xir, &xii);
    for (i = x0; i <= x1; i++) {
      if (i < 0 || i >= domain->nx) {
	continue;
      }
      for (j = y0; j <= y1; j++) {
	if (j >= 0 && j < domain->ny) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 6;
    var += 6;
  }
  return MW_SUCCESS;
}

/* Add one or more convex lenses as determined by the vector "var" of
   length "ivar", where each group of six elements corresponds to: (0)
   x of centre, (1) y of centre, (2) radius of curvature, (3) size,
   (4) epsilon_r, (5) epsilon_i */
int
mw_add_lens(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 5) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real radcurv = var[2]/domain->dx;
    real radius = var[3]/domain->dx;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[4], var[5], &xir, &xii);
    for (i = x0-radius; i <= x0+radius; i++) {
      int thickness = 2+radcurv - sqrt(radcurv*radcurv
				     -radius*radius+(i-x0)*(i-x0));
      if (i < 0 || i >= domain->nx) {
	continue;
      }
      for (j = y0-thickness; j < y0; j++) {
	if (j >= 0 && j < domain->ny) {
	  domain->epsilon[j][i] += xir;
	  domain->Edamping[j][i] += xii;
	}
      }
    }
    nvar -= 6;
    var += 6;
  }
  return MW_SUCCESS;
}

/* Add one or more cavities (a rectangular block with a circular
   section removed), as determined by the vector "var" of length
   "ivar", where each group of six elements corresponds to: (0) x of
   bottom left, (1) y of bottom left, (2) x of top right, (3) y of top
   right), (4) x of centre of circle, (5) y of centre of circle, (6)
   radius of curvature,(7) epsilon_r, (8) epsilon_i */
int
mw_add_cavity(mwDomain *domain, int nvar, real *var)
{
  while (nvar > 8) {
    real x0 = var[0]/domain->dx + domain->nx/2.0;
    real y0 = var[1]/domain->dx + domain->ny/2.0;
    real x1 = var[2]/domain->dx + domain->nx/2.0;
    real y1 = var[3]/domain->dx + domain->ny/2.0;
    real xc = var[4]/domain->dx + domain->nx/2.0;
    real yc = var[5]/domain->dx + domain->ny/2.0;
    real radius = var[6]/domain->dx;
    real radius2 = radius*radius;
    real xir, xii;
    int i, j;
    mw_susceptibility(var[7], var[8], &xir, &xii);
    for (i = x0; i <= x1; i++) {
      if (i < 0 || i >= domain->nx) {
	continue;
      }
      for (j = y0; j <= y1; j++) {
	if (j >= 0 && j < domain->ny) {
	  if ((xc-i)*(xc-i) + (yc-j)*(yc-j)
	      > radius2) {
	    domain->epsilon[j][i] += xir;
	    domain->Edamping[j][i] += xii;
	  }
	}
      }
    }
    nvar -= 9;
    var += 9;
  }
  return MW_SUCCESS;
}
