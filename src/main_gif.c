/* main_gif.c -- Program code for maxwell2d_gif to create gif
   animations of electromagnetic waves

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
#include "maxwell.h"

int
main(int argc, char **argv)
{
  mwDomain domain;
  char *epsilon_plot_file = NULL;

  /* Initialize the domain based on command-line arguments and
     standard input */
  if (mw_start(argc, argv, &domain)) {
    exit(1);
  }

  /* Write out a plot of the dielectric constant if required */
  rc_assign_string(domain.config, "epsilon_gif_file", 
		   &epsilon_plot_file);
  if (epsilon_plot_file) {
    mw_gif_write_epsilon(epsilon_plot_file, &domain);
  }

  /* Initialize a gif file to be written to standard output */
  mw_gif_init(NULL, &domain);

  /* Continue simulation until the total required time has elapsed */
  while (domain.time < domain.duration) {
    /* Write a frame to a gif file */
    mw_gif_write_frame(&domain);
    fprintf(stderr, ".");

    /* Move the simulation forward one frame (7 timesteps) */
    mw_frame(&domain);
  }
  fprintf(stderr, "\n");
  mw_gif_close();
  exit(0);
}
