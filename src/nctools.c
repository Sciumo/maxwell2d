/* nctools.c -- useful NetCDF functions 

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
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>

#include <netcdf.h>

#define MAX_STRING_LENGTH 512

static int ncstatus;
#define CHECK(cmd) if ((ncstatus = (cmd)) != NC_NOERR) { return ncstatus; }

/* Prepend "string" to the attribute "attname" of the variable with ID
 * "varid" (where NC_GLOBAL indicates a global attribute) of NetCDF
 * dataset with ID "ncid", using separator "separator". */
int 
nct_prepend_attribute(int ncid, int varid, char *attname, char *string, char *separator)
{ 
  size_t oldlen, newlen = strlen(string);
  int status;

  status = nc_inq_attlen(ncid, varid, attname, &oldlen);

  if (status == NC_ENOTATT) {
    /* New attribute */
    return nc_put_att_text(ncid, varid, attname, newlen, string);
  }
  else if (status == NC_NOERR) {
    /* Existing attribute */
    size_t fulllen = oldlen + newlen + strlen(separator) + 1;
    char *newstring = malloc(fulllen * sizeof(char));
    if (!newstring) {
      return NC_ENOMEM;
    }
    sprintf(newstring, "%s%s", string, separator);
    status = nc_get_att_text(ncid, varid, attname, newstring+strlen(newstring));
    if (status != NC_NOERR) {
      free(newstring);
      return status;
    }
    newstring[fulllen-1] = '\0';
    status = nc_put_att_text(ncid, varid, attname, strlen(newstring), newstring);
    free(newstring);
    return status;
  }
  else {
    /* An error occurred */
    return status;
  }
}

/* Append "string" to the attribute "attname" of the variable with ID
 * "varid" (where NC_GLOBAL indicates a global attribute) of NetCDF
 * dataset with ID "ncid", using separator "separator". */
int
nct_append_attribute(int ncid, int varid, char *attname, char *string, char *separator)
{
  size_t oldlen, newlen = strlen(string);
  int status;

  status = nc_inq_attlen(ncid, varid, attname, &oldlen);

  if (status == NC_ENOTATT) {
    /* New attribute */
    return nc_put_att_text(ncid, varid, attname, newlen, string);
  }
  else if (status == NC_NOERR) {
    /* Existing attribute */
    size_t fulllen = oldlen + newlen + strlen(separator) + 1;
    char *newstring = malloc(fulllen * sizeof(char));
    if (!newstring) {
      return NC_ENOMEM;
    }
    status = nc_get_att_text(ncid, varid, attname, newstring);
    newstring[oldlen] = '\0';
    sprintf(newstring+strlen(newstring), "%s%s", separator, string);

    if (status != NC_NOERR) {
      free(newstring);
      return status;
    }
    status = nc_put_att_text(ncid, varid, attname, strlen(newstring), newstring);
    free(newstring);
    return status;
  }
  else {
    /* An error occurred */
    return status;
  }
}

/* Append information to the "history" global attribute of a NetCDF
 * dataset.  The information is of the form "$action at $time by $user
 * on $host", where action and user are given as arguments (a value of
 * NULL for user causes the username to be used), time is the current
 * time and host is the name of the machine. Histories are separated
 * by semicolons. */
int
nct_add_history(int ncid, char *action, char *user)
{
  struct timeval tv;
  char *name;
  char hostname[MAX_STRING_LENGTH] = "unknown";
  char history[MAX_STRING_LENGTH] = "";

  if (gettimeofday(&tv, NULL)) {
    name = "";
  }
  else {
    name = ctime(&tv.tv_sec);
    name[24] = ' '; /* Get rid of the newline */
  }

  gethostname(hostname, MAX_STRING_LENGTH);

  if (!user) {
    /* user string is NULL; get the username instead */
    if (!(user = getenv("LOGNAME"))) {
      if (!(user = getenv("USER"))) {
	user = "anonymous";
      }    
    }
  }
  
  snprintf(history, MAX_STRING_LENGTH, "%s- %s by %s on %s", name, action, user, hostname);
  history[MAX_STRING_LENGTH-1] = '\0';

  return nct_append_attribute(ncid, NC_GLOBAL, "history", history, "\n");
}

/* Append the full command line to the global attribute "command_line"
 * of a NetCDF file, using semicolons as separators if several
 * programs write to the same file. */
#define COMPARE_ARG_LENGTH 4
#define NUM_SIMILAR_ARGS 5
int
nct_add_command_line(int ncid, int argc, char **argv)
{
  int i, j, len = 0;
  char *cmdline;
  char lastarg[COMPARE_ARG_LENGTH+1];
  int status;
  int numsimilarargs = 1;
  
  if (argv == NULL) {
    return NC_EINVAL;
  }
  for (i = 0; i < argc; i++) {
    char *c = argv[i];
    char need_quotes = 0;
    for (j = 0; c[j] != '\0'; j++) {
      if (c[j] <= ' ') {
	need_quotes = 1;
      }
    }
    len += (1 + j + need_quotes * 2);
  }

  cmdline = malloc(len * sizeof(char));
  if (cmdline == NULL) {
    return NC_ENOMEM;
  }

  len = 0;
  lastarg[0] = '\0';
  for (i = 0; i < argc; i++) {
    char *c = argv[i];
    char need_quotes = 0;
    for (j = 0; c[j] != '\0'; j++) {
      if (c[j] <= ' ') {
	need_quotes = 1;
      }
    }
    /* This next bit of code ensures that if more than
       NUM_SIMILAR_ARGS are found that share the same first
       COMPARE_ARG_LENGTH characters, then an elipsis (...) will be
       output.  This saves on excessively long command lines being
       saved. */
    if (j > COMPARE_ARG_LENGTH) {
      if (c[0] != '-' && strncmp(lastarg, c, COMPARE_ARG_LENGTH) == 0) {
	numsimilarargs++;
	if (numsimilarargs == NUM_SIMILAR_ARGS+1) {
	  c = "...";
	  j = 3;
	}
	else if (numsimilarargs > NUM_SIMILAR_ARGS+1) {
	  continue;
	}
      }
      else {
	snprintf(lastarg, COMPARE_ARG_LENGTH+1, "%s", c);
	lastarg[COMPARE_ARG_LENGTH] = '\0';
	numsimilarargs = 1;
      }
    }
    else {
      numsimilarargs = 1;
    }

    if (need_quotes) {
      sprintf(cmdline + len, "\"%s\"", c);
    }
    else {
      sprintf(cmdline + len, "%s", c);
    }
    len += (1 + j + need_quotes * 2);
    cmdline[len -1] = ' ';
  }
  cmdline[len-1] = '\0';

  status = nct_append_attribute(ncid, NC_GLOBAL, "command_line", cmdline, "\n");
  free(cmdline);

  return status;
}

/* Remove trailing abreviations in parantheses from an attribute */
int
nct_strip_parentheses(int ncid, int varid, char *attname)
{
  size_t len;
  int status = nc_inq_attlen(ncid, varid, attname, &len);

  if (status == NC_NOERR) {
    char *str = malloc(len * sizeof(char));
    char *ch;
    if (!str) {
      return NC_ENOMEM;
    }
    CHECK(nc_get_att_text(ncid, varid, attname, str));
    ch = str + len - 1;
    if (*ch == '\0' && len > 4) {
      --ch;
    }
    if (*ch == ')') {
      while ((--ch) > str) {
	if (*ch == '(') {
	  ch--;
	  if (*ch == ' ') {
	    *ch = '\0';
	  }
	  break;
	}
      }
    }
    if (*ch == '\0') {
      nc_put_att_text(ncid, varid, attname, strlen(str), str);
    }
    free(str);
    return NC_NOERR;
  }
  else {
    /* An error occurred */
    return status;
  }
  
}

int
nct_add_string_attribute(int ncid, int varid, char *attname, char *value)
{
  if (value) {
    return nc_put_att_text(ncid, varid, attname,
			   strlen(value), value);
  }
  else {
    return NC_NOERR;
  }
}

int
nct_add_standard_attributes(int ncid, int varid, char *units,
			    char *long_name, char *standard_name,
			    char *comment)
{
  if (units) {
    CHECK(nct_add_string_attribute(ncid, varid, "units", units));
  }
  if (long_name) {
    CHECK(nct_add_string_attribute(ncid, varid, "long_name", long_name));
  }
  if (standard_name) {
    CHECK(nct_add_string_attribute(ncid, varid,
				   "standard_name", standard_name));
  }
  if (comment) {
    CHECK(nct_add_string_attribute(ncid, varid, "comment", comment));
  }
  return NC_NOERR;
}

int
nct_add_plot_attributes(int ncid, int varid, char *units_html,
			 float plot_min, float plot_max, int is_log_scale)
{
  if (units_html) {
    CHECK(nct_add_string_attribute(ncid, varid, "units_html", units_html));
  }
  if (plot_min < plot_max) {
    float plot_range[2] = {plot_min, plot_max};
    nc_type xtype;
    CHECK(nc_inq_vartype(ncid, varid, &xtype));
    CHECK(nc_put_att_float(ncid, varid, "plot_range", xtype, 2, plot_range));
    if (is_log_scale) {
      CHECK(nct_add_string_attribute(ncid, varid,
				     "plot_scale", "logarithmic"));
    }
    else {
      CHECK(nct_add_string_attribute(ncid, varid, "plot_scale", "linear"));
    }
  }
  return NC_NOERR;
}

int
nct_add_missing_value(int ncid, int varid, float value)
{
  CHECK(nc_put_att_float(ncid, varid, "missing_value",
			 NC_FLOAT, 1, &value));
  return nc_put_att_float(ncid, varid, "_FillValue",
			  NC_FLOAT, 1, &value);
}

int
nct_add_missing_value_type(int ncid, int varid, int type, float value)
{
  CHECK(nc_put_att_float(ncid, varid, "missing_value",
			 type, 1, &value));
  return nc_put_att_float(ncid, varid, "_FillValue",
			  type, 1, &value);
}

int
nct_copy_attributes(int ncid_src, int varid_src,
		    int ncid_dest, int varid_dest)
{
  int iatt, natts;
  char attname[NC_MAX_NAME];
  CHECK(nc_inq_varnatts(ncid_src, varid_src, &natts))

  for (iatt = 0; iatt < natts; iatt++) {
    CHECK(nc_inq_attname(ncid_src, varid_src, iatt, attname));
    CHECK(nc_copy_att(ncid_src, varid_src, attname,
		      ncid_dest, varid_dest));
  }
  return NC_NOERR;
}

int
nct_copy_attributes_by_name(int ncid_src, char *varname_src,
			    int ncid_dest, char *varname_dest)
{
  int varid_src, varid_dest;
  CHECK(nc_inq_varid(ncid_src, varname_src, &varid_src));
  CHECK(nc_inq_varid(ncid_dest, varname_dest, &varid_dest));
  CHECK(nct_copy_attributes(ncid_src, varid_src,
			    ncid_dest, varid_dest));
  return NC_NOERR;
}
