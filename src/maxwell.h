/* maxwell.h -- Header file for E-M wave simulator

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

#ifndef _MAXWELL_H
#define _MAXWELL_H 1

#ifdef __cplusplus
extern "C" {
#endif                          /* __cplusplus */

#include <stdio.h>
#include "readconfig.h"

#define MW_CHECK(a) if ((a) != MW_SUCCESS) { \
  fprintf(stderr, "Error at line %d of " __FILE__ "\n", __LINE__); \
  return MW_FAILURE; }

/* For double-precision, change this to double */
#define real float

/* Error codes reported by functions */
#define MW_SUCCESS 0
#define MW_FAILURE 1

/* Speed of light */
#define MW_C 2.997924583e8

/* Various domain modes: "VACUUM" means that a parallel simulation is
   done in vacuum such that the scattered field can be calculated by
   subtracting the field in a vacuum from the total field. The "EZ"
   and "EXY" determine whether the Z-component of the E-field is
   simulated and/or the X and Y components. */
#define MW_MODE_VACUUM (1L<<0)
#define MW_MODE_EZ (1L<<1)
#define MW_MODE_EXY (1L<<2)

/* The number of timesteps in a frame */
#define MW_MINOR_STEPS 7

/* The mwDomain structure */
  typedef struct {
    real **Ex;
    real **Ey;
    real **Ez;
    real **Bx;
    real **By;
    real **Bz;
    real **Ex_vacuum;
    real **Ey_vacuum;
    real **Ez_vacuum;
    real **Bx_vacuum;
    real **By_vacuum;
    real **Bz_vacuum;
    real **epsilon;
    real **Edamping;
    real **Bdamping;
    real **forcingI;
    real **forcingQ;
    real **Eprefix;
    real **scat_field;
    real **Poynting_x;
    real **Poynting_y;
    real **Poynting_x_scat;
    real **Poynting_y_scat;
    real **boundaries;
    real *frequencies;
    char *epsilon_plot_file;
    rc_data *config;
    real Ex_forcingI;
    real Ey_forcingI;
    real Ez_forcingI;
    real Ex_forcingQ;
    real Ey_forcingQ;
    real Ez_forcingQ;
    real Ex_amplitude;
    real Ey_amplitude;
    real Ez_amplitude;
    real dx;
    real dt;
    real dt_dx;
    real c;
    real primary_frequency;
    real time;
    real plot_E_max;
    real plot_B_max;
    real plot_scat_ratio;
    real duration;
    int nx;
    int ny;
    int mode;
    int borderwidth;
    int cycles;
    int iframe;
    int mag;
    int nfrequencies;
  } mwDomain;

  /* Functions */
  int mw_subtract(int nx, int ny, real **arg1, real **arg2, real **ans);
  int mw_scale(int nx, int ny, real **arg, real factor);

  int mw_start(int argc, char **argv, mwDomain *domain);
  int mw_frame(mwDomain *domain);

  int mw_new_field(real ***field, int nx, int ny, int value);
  int mw_free_field(real **field);

  int mw_new_domain(mwDomain *domain, int nx, int ny, real dx, int mode);
  int mw_free_domain(mwDomain *domain);

  int mw_reset_field(real **field, int nx, int ny, real value);
  int mw_reset_damping(mwDomain *domain, int borderwidth);
  int mw_add_circle(mwDomain *domain, int nvar, real *var);
  int mw_add_edge(mwDomain *domain, int nvar, real *var);
  int mw_add_gradient(mwDomain *domain, int nvar, real *var);
  int mw_add_ripple(mwDomain *domain, int nvar, real *var);
  int mw_add_dish(mwDomain *domain, int nvar, real *var);
  int mw_add_rectangle(mwDomain *domain, int nvar, real *var);
  int mw_add_rotated_rectangle(mwDomain *domain, int nvar, real *var);
  int mw_add_wave_packet(mwDomain *domain, int nvar, real *var);
  int mw_add_lens(mwDomain *domain, int nvar, real *var);
  int mw_add_cavity(mwDomain *domain, int nvar, real *var);

  int mw_print_field(FILE *file, real **field, int nx, int ny);

  int mw_step(mwDomain *domain);

  int mw_nc_init(char *filename, mwDomain *domain, int argc, char **argv);
  int mw_nc_write_frame(mwDomain *domain);
  int mw_nc_close();

  int mw_gif_init(char *filename, mwDomain *domain);
  int mw_gif_write_frame(mwDomain *domain);
  int mw_gif_close();
  int mw_gif_write_epsilon(char *filename, mwDomain *domain);

  int mw_susceptibility(real nr, real ni, real *xir, real *xii);
  int mw_find_boundaries(mwDomain *domain);


#ifdef __cplusplus
}                               /* extern "C" */
#endif                          /* __cplusplus */

#endif /* maxwell.h */
