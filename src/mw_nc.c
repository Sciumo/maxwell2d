/* mw_nc.c -- Functions to write the contents of NetCDF files

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
#include <string.h>
#include <netcdf.h>
#include "maxwell.h"
#include "nctools.h"

static int ncstatus;
#define NC_CHECK(a) if ((ncstatus = (a)) != NC_NOERR) { \
  return MW_FAILURE; }

/* NetCDF variable IDs are stored as static variables within this
   file; that means that only one NetCDF file may safely be written at
   a time */
static int ncid;
static int xdimid, ydimid, timedimid;
static int timeid;
//static int Exid, Eyid;
//static int Bxid, Byid;
static int Ezid, Bzid;
//static int Exscatid, Eyscatid;
//static int Bxscatid, Byscatid;
static int Ezscatid, Bzscatid;
static int Sxid, Syid;
static int Sxscatid, Syscatid;
static int nc_skip = 0;

/* Add some standard attributes to a variable */
static
int
add_attributes(int ncid, int fieldid, char *units,
	       char *long_name, char *comment)
{
  if (units) {
    MW_CHECK(nct_add_string_attribute(ncid, fieldid, "units", units));
  }
  if (long_name) {
    MW_CHECK(nct_add_string_attribute(ncid, fieldid, "long_name", long_name));
  }
  if (comment) {
    MW_CHECK(nct_add_string_attribute(ncid, fieldid, "comment", comment));
  }
  return MW_SUCCESS;
}

/* Put a 2D field in the file */
static
int
put_field(int ncid, int fieldid, real **M, int nx, int ny)
{
  size_t start[2], count[2];

  start[0] = 0;
  start[1] = 0;
  count[0] = ny;
  count[1] = nx;

  NC_CHECK(nc_put_vara_float(ncid, fieldid, start, count, M[0]));

  return MW_SUCCESS;
}

/* Put one 2D slice of a 3D field into the file */
static
int
put_slice(int ncid, int fieldid, real **M, int nx, int ny, int iframe)
{
  size_t start[3], count[3];

  start[0] = iframe;
  start[1] = 0;
  start[2] = 0;
  count[0] = 1;
  count[1] = ny;
  count[2] = nx;

  NC_CHECK(nc_put_vara_float(ncid, fieldid, start, count, M[0]));

  return MW_SUCCESS;
}

/* Initialize the NetCDF file */
int
mw_nc_init(char *filename, mwDomain *domain, int argc, char **argv)
{
  char *confstring = NULL;
  char *title = NULL;
  int epsilon_r_id, epsilon_i_id;
  int dimids[3];
  
  /* If we are only interested in the Poynting vector and the
     dielectric constant then this option will result in the
     time-dependent fields not being stored */
  nc_skip = rc_get_boolean(domain->config,
			   "nc_skip_time_dependent_fields");
  /* Open new file */
  NC_CHECK(nc_create(filename, NC_CLOBBER, &ncid));

  /* Set the dimensions */
  if (!nc_skip) {
    NC_CHECK(nc_def_dim(ncid, "time", NC_UNLIMITED, &timedimid));
  }
  NC_CHECK(nc_def_dim(ncid, "y", domain->ny, &ydimid));
  NC_CHECK(nc_def_dim(ncid, "x", domain->nx, &xdimid));

  /* Define the variables */
  dimids[0] = timedimid;
  dimids[1] = ydimid;
  dimids[2] = xdimid;

  if (!nc_skip) {
    NC_CHECK(nc_def_var(ncid, "time", NC_FLOAT,
			1, dimids, &timeid));
  }
  NC_CHECK(nc_def_var(ncid, "epsilon_r", NC_FLOAT, 
		      2, &dimids[1], &epsilon_r_id));
  NC_CHECK(nc_def_var(ncid, "epsilon_i", NC_FLOAT, 
		      2, &dimids[1], &epsilon_i_id));
  if (!nc_skip) {
    NC_CHECK(add_attributes(ncid, timeid, "s", 
			    "Time since start of simulation", NULL));
  }
  NC_CHECK(add_attributes(ncid, epsilon_r_id, "1", 
			  "Real part of the dielectric constant", NULL));
  NC_CHECK(add_attributes(ncid, epsilon_i_id, "1", 
			  "Imaginary part of the dielectric constant", 
			  "Note that this field is positive for ordinary materials and the full dielectric constant is given by epsilon_r-i*epsilon_i"));

  NC_CHECK(nc_def_var(ncid, "Sx", NC_FLOAT, 2, &dimids[1], &Sxid));
  NC_CHECK(add_attributes(ncid, Sxid, "W m-2",
	  "Mean x-component of Poynting vector for total field", NULL));
  NC_CHECK(nc_def_var(ncid, "Sy", NC_FLOAT, 2, &dimids[1], &Syid));
  NC_CHECK(add_attributes(ncid, Syid, "W m-2",
	  "Mean y-component of Poynting vector for total field", NULL));

  if (domain->mode & MW_MODE_VACUUM) {
    NC_CHECK(nc_def_var(ncid, "Sx_scat", NC_FLOAT, 2, &dimids[1], 
			&Sxscatid));
    NC_CHECK(add_attributes(ncid, Sxscatid, "W m-2",
    "Mean x-component of Poynting vector for scattered field", NULL));
    NC_CHECK(nc_def_var(ncid, "Sy_scat", NC_FLOAT, 2, &dimids[1], 
			&Syscatid));
    NC_CHECK(add_attributes(ncid, Syscatid, "W m-2",
    "Mean y-component of Poynting vector for scattered field", NULL));
  }

  if (!nc_skip) {
    if (domain->mode & MW_MODE_EZ) {
      NC_CHECK(nc_def_var(ncid, "Ez", NC_FLOAT, 3, dimids, &Ezid));
      NC_CHECK(add_attributes(ncid, Ezid, "V m-1",
	      "Z-component of the total electric field", NULL));
    }
    if (domain->mode & MW_MODE_EXY) {
      NC_CHECK(nc_def_var(ncid, "Bz", NC_FLOAT, 3, dimids, &Bzid));
      NC_CHECK(add_attributes(ncid, Bzid, "T",
	      "Z-component of the total magnetic field", NULL));
    }
    if (domain->mode & MW_MODE_EZ && domain->mode & MW_MODE_VACUUM) {
      NC_CHECK(nc_def_var(ncid, "Ez_scat", NC_FLOAT, 3, dimids, &Ezscatid));
      NC_CHECK(add_attributes(ncid, Ezscatid, "V m-1",
	      "Z-component of the scattered electric field",
	      "This field is simply the total electric field minus the electric field that would have occurred if the same electromagnetic wave had occurred in a vacuum"));
      
    }
    if (domain->mode & MW_MODE_EXY && domain->mode & MW_MODE_VACUUM) {
      NC_CHECK(nc_def_var(ncid, "Bz_scat", NC_FLOAT, 3, dimids, &Bzscatid));
      NC_CHECK(add_attributes(ncid, Bzscatid, "T",
	      "Z-component of the scattered magnetic field",
	      "This field is simply the total magnetic field minus the magnetic field that would have occurred if the same electromagnetic wave had occurred in a vacuum"));
    }
  }

  /* Define some global attributes */
  rc_assign_string(domain->config, "title", &title);
  if (title) {
    NC_CHECK(nct_add_string_attribute(ncid, NC_GLOBAL, "title", title));
    free(title);
  }

  NC_CHECK(nct_add_command_line(ncid, argc, argv));
  NC_CHECK(nct_add_history(ncid, "Maxwell2D simulation performed", NULL))
  confstring = rc_sprint(domain->config);
  if (confstring) {
    NC_CHECK(nc_put_att_text(ncid, NC_GLOBAL, "config",
			     strlen(confstring),
			     confstring));
    free(confstring);
  }

  /* End define mode */
  NC_CHECK(nc_enddef(ncid));

  /* Write the time-independent fields */
  NC_CHECK(put_field(ncid, epsilon_r_id, domain->epsilon,
		     domain->nx, domain->ny));
  NC_CHECK(put_field(ncid, epsilon_i_id, domain->Edamping,
		     domain->nx, domain->ny));
  return MW_SUCCESS;
}

/* Write a frame of the NetCDF file */
int
mw_nc_write_frame(mwDomain *domain)
{
  size_t index = domain->iframe;

  if (nc_skip) {
    return MW_SUCCESS;
  }

  NC_CHECK(nc_put_var1_float(ncid, timeid, &index, &domain->time));

  if (domain->mode & MW_MODE_EZ) {
    NC_CHECK(put_slice(ncid, Ezid, domain->Ez,
		       domain->nx, domain->ny, domain->iframe));
  }
  if (domain->mode & MW_MODE_EXY) {
    NC_CHECK(put_slice(ncid, Bzid, domain->Bz,
		       domain->nx, domain->ny, domain->iframe));
  }
  if (domain->mode & MW_MODE_EZ && domain->mode & MW_MODE_VACUUM) {
    mw_subtract(domain->nx, domain->ny, domain->Ez, domain->Ez_vacuum,
		domain->scat_field);
    NC_CHECK(put_slice(ncid, Ezscatid, domain->scat_field,
		       domain->nx, domain->ny, domain->iframe));
  }
  if (domain->mode & MW_MODE_EXY && domain->mode & MW_MODE_VACUUM) {
    mw_subtract(domain->nx, domain->ny, domain->Bz, domain->Bz_vacuum,
		domain->scat_field);
    NC_CHECK(put_slice(ncid, Bzscatid, domain->scat_field,
		       domain->nx, domain->ny, domain->iframe));
  }
  return MW_SUCCESS;
}

/* Write the mean Poynting vector to the file then close it */
int
mw_nc_close(mwDomain *domain)
{
  /* Currently the Poynting vector contains the sum of the values from
     each frame, so it needs to be scaled to obtain the mean */
  mw_scale(domain->nx, domain->ny, domain->Poynting_x, 1.0/domain->iframe);
  mw_scale(domain->nx, domain->ny, domain->Poynting_y, 1.0/domain->iframe);
  NC_CHECK(put_field(ncid, Sxid, domain->Poynting_x,
		     domain->nx, domain->ny));
  NC_CHECK(put_field(ncid, Syid, domain->Poynting_y,
		     domain->nx, domain->ny));
  if (domain->mode & MW_MODE_VACUUM) {
    mw_scale(domain->nx, domain->ny, domain->Poynting_x_scat,
	     1.0/domain->iframe);
    mw_scale(domain->nx, domain->ny, domain->Poynting_y_scat,
	     1.0/domain->iframe);
    NC_CHECK(put_field(ncid, Sxscatid, domain->Poynting_x_scat,
		       domain->nx, domain->ny));
    NC_CHECK(put_field(ncid, Syscatid, domain->Poynting_y_scat,
		       domain->nx, domain->ny));
  }

  NC_CHECK(nc_close(ncid));
  return MW_SUCCESS;
}
