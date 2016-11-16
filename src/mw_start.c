/* mw_start.c -- Initialize the domain for the simulation

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
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include "maxwell.h"
#include "readconfig.h"

/* Initialize the domain for the simulation based on the command-line
   arguments and any config files on standard input */
int
mw_start(int argc, char **argv, mwDomain *domain)
{
  int k;
  rc_data *config;

  int nx = 64;
  int ny = 64;
  real dx = 1.0;
  int borderwidth = 6.0;
  char *polarization = "z";
  int mode = 0;
  int vacuum = 0;
  real *line_osc;
  int n_line_osc;
  real *point_osc;
  real *var;
  int n_var;
  //  char *epsilon_plot_file = NULL;

  /* Find the first config file on the command line. */
  int ifile = rc_get_file(argc, argv);

  if (!ifile) {
    /* No file given - assume command-line arguments contain all the
       information. */
    config = rc_read(NULL, stderr);
  }
  else {
    /* Read configuration information from the file. */
    config = rc_read(argv[ifile], stderr);
  }
  if (!config) {
    fprintf(stderr, "Error initializing configuration information\n");
    return MW_FAILURE;
  }
  
  /* Supplement configuration information with command-line
     arguments. */
  rc_register_args(config, argc, argv);

  rc_assign_int(config, "x_pixels", &nx);
  rc_assign_int(config, "y_pixels", &ny);
  rc_assign_real(config, "pixel_spacing", &dx);
  rc_assign_int(config, "border_width", &borderwidth);
  rc_assign_string(config, "polarization", &polarization);
  if (strcasecmp(polarization, "xyz") == 0) {
    mode |= (MW_MODE_EXY | MW_MODE_EZ);
  }
  else if (strcasecmp(polarization, "xy") == 0) {
    mode |= (MW_MODE_EXY);
  }
  else if (strcasecmp(polarization, "z") == 0) {
    mode |= (MW_MODE_EZ);
  }
  else {
    fprintf(stderr, "Config variable \"polarization\" must be \"xy\", \"z\" or \"xyz\"\n");
    return MW_FAILURE;
  }

  vacuum = rc_get_boolean(config, "vacuum");
  if (vacuum) {
    mode |= MW_MODE_VACUUM;
    mw_new_field(&domain->scat_field, nx, ny, 0.0);
  }

  mw_new_domain(domain, nx, ny, dx, mode);
  mw_reset_damping(domain, borderwidth);

  domain->primary_frequency = 0.1*MW_C;
  domain->Ex_amplitude = 0.0;
  domain->Ey_amplitude = 0.0;
  domain->Ez_amplitude = 1.0;
  domain->cycles = 10;
  domain->duration = 200.0*MW_MINOR_STEPS*domain->dt;

  rc_assign_real(config, "frequency", &domain->primary_frequency);
  if ((domain->frequencies
       = rc_get_real_vector(config, "frequencies",
			    &n_var)) && n_var > 2) {
    domain->nfrequencies = n_var / 3;
    domain->primary_frequency = domain->frequencies[0];    
  }

  rc_assign_real(config, "x_amplitude", &domain->Ex_amplitude);
  rc_assign_real(config, "y_amplitude", &domain->Ey_amplitude);
  rc_assign_real(config, "z_amplitude", &domain->Ez_amplitude);
  rc_assign_int(config, "cycles", &domain->cycles);
  rc_assign_int(config, "mag", &domain->mag);
  rc_assign_real(config, "plot_E_max", &domain->plot_E_max);
  rc_assign_real(config, "plot_B_max", &domain->plot_B_max);
  rc_assign_real(config, "plot_scat_ratio", &domain->plot_scat_ratio);
  rc_assign_real(config, "dx", &domain->dx);
  rc_assign_real(config, "duration", &domain->duration);

  //  domain->Ez_forcing = 1.0;

  if ((line_osc = rc_get_real_vector(config, "line_oscillator",
				     &n_line_osc)) && n_line_osc > 1) {
    for (k = 0; k < domain->nx; k++) {
      domain->forcingI[borderwidth+1][k]
	= line_osc[0]*exp(-pow(((real)k-domain->nx/2.0)*2.0
			       /(line_osc[1]*domain->nx), 4.0));
    }
  }

  if ((var = rc_get_real_vector(config, "point_oscillator",
				&n_var))) {
    point_osc = var;
    while (n_var > 2) {
      real x0 = point_osc[1]/domain->dx + domain->nx/2.0;
      real y0 = point_osc[2]/domain->dx + domain->ny/2.0;

      if (x0 > 0 && x0 < nx-1 && y0 > 0 && y0 < ny-1) {
	domain->forcingI[(int)y0][(int)x0] = point_osc[0];
      }
      n_var -= 3;
      point_osc += 3;
    }
    free(var);
  }

  if ((var = rc_get_real_vector(config, "phased_point_oscillator",
				&n_var))) {
    point_osc = var;
    while (n_var > 3) {
      real x0 = point_osc[1]/domain->dx + domain->nx/2.0;
      real y0 = point_osc[2]/domain->dx + domain->ny/2.0;

      if (x0 > 0 && x0 < nx-1 && y0 > 0 && y0 < ny-1) {
	domain->forcingI[(int)y0][(int)x0] = point_osc[0] 
	  * cos(M_PI*point_osc[3]/180.0);
	domain->forcingQ[(int)y0][(int)x0] = point_osc[0]
	  * sin(M_PI*point_osc[3]/180.0);
      }
      n_var -= 4;
      point_osc += 4;
    }
    free(var);
  }

  if ((var = rc_get_real_vector(config, "circle",
				   &n_var)) && n_var > 3) {
    mw_add_circle(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "edge",
				   &n_var)) && n_var > 3) {
    mw_add_edge(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "ripple",
				   &n_var)) && n_var > 3) {
    mw_add_ripple(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "gradient",
				   &n_var)) && n_var > 3) {
    mw_add_gradient(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "dish",
				   &n_var)) && n_var > 3) {
    mw_add_dish(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "rectangle",
				   &n_var)) && n_var > 3) {
    mw_add_rectangle(domain, n_var, var);
    free(var);
  }


  if ((var = rc_get_real_vector(config, "rotated_rectangle",
				&n_var)) && n_var > 3) {
    mw_add_rotated_rectangle(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "wave_packet",
				&n_var)) && n_var > 3) {
    mw_add_wave_packet(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "lens",
				   &n_var)) && n_var > 3) {
    mw_add_lens(domain, n_var, var);
    free(var);
  }

  if ((var = rc_get_real_vector(config, "cavity",
				   &n_var)) && n_var > 3) {
    mw_add_cavity(domain, n_var, var);
    free(var);
  }


  /*
  if (epsilon_plot_file) {
    mw_gif_write_epsilon(epsilon_plot_file, domain);
    free(epsilon_plot_file);
  }
  */

  //  rc_clear(config);
  domain->config = config;

  mw_find_boundaries(domain);

  return MW_SUCCESS;
}
