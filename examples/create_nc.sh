#!/bin/bash

if [ "$#" = 0 ]
then
  echo "Usage:"
  echo "  $0 file1.cfg [file2.cfg ...]"
  echo "If you have a lot of time on your hands you could create"
  echo "netcdf files from all the config files in this directory:"
  echo "  $0 *.cfg"
  exit
fi


# Decide which component(s) of the electric field will be simulated
POL=z
#POL=xy
#POL=xyz

# Uncomment this if you want the netcdf files to contain only the
# dielectric constant and the mean Poynting vector
#OPTIONS="nc_skip_time_dependent_fields=1"

# Loop through all command-line arguments, treating each as a config
# file
for CFGFILE in $@
do
  if [ ! -r $CFGFILE ]
  then
    echo "Error: \"$CFGFILE\" is not a readable file"
    exit 1
  fi

  # Extract the prefix from the config file name
  PREFIX=$(echo $CFGFILE | sed 's/\.cfg$//')

  echo "Reading $CFGFILE, writing ${PREFIX}.nc"
  cat default/domain.cfg default/$POL.cfg $CFGFILE \
      | ../src/maxwell2d_nc $OPTIONS nc_file=${PREFIX}.nc -
done
