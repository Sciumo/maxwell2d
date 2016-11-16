/* nctools.h -- useful NetCDF functions

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

/* See nctools.c to see what each function does */

#ifdef __cplusplus
extern "C" {
#endif

int nct_prepend_attribute(int ncid, int varid, char *attname,
			  char *string, char *separator);
int nct_append_attribute(int ncid, int varid, char *attname,
			 char *string, char *separator);
int nct_strip_parentheses(int ncid, int varid, char *attname);
int nct_add_history(int ncid, char *action, char *user);
int nct_add_command_line(int ncid, int argc, char **argv);
int nct_add_string_attribute(int ncid, int fieldid,
			     char *attname, char *value);
int nct_add_standard_attributes(int ncid, int fieldid, char *units,
				char *long_name, char *standard_name,
				char *comment);
int nct_add_missing_value(int ncid, int fieldid, float value);
int nct_add_missing_value_type(int ncid, int varid, int type, float value);
int nct_copy_attributes(int ncid_src, int vardid_src,
			int ncid_dest, int varid_dest);
int nct_copy_attributes_by_name(int ncid_src, char *varname_src,
				int ncid_dest, char *varname_dest);
int nct_add_plot_attributes(int ncid, int varid, char *units_html,
			    float plot_min, float plot_max, int is_log_scale);

#ifdef __cplusplus
}
#endif
