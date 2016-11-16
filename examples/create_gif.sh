#!/bin/bash

if [ "$#" = 0 ]
then
  echo "Usage:"
  echo "  $0 file1.cfg [file2.cfg ...]"
  echo "If you have a lot of time on your hands you could create"
  echo "animated gifs from all the config files in this directory:"
  echo "  $0 *.cfg"
  exit
fi


# Decide which component(s) of the electric field will be simulated
POL=z
#POL=xy
#POL=xyz

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

  echo "Reading $CFGFILE, writing ${PREFIX}_$POL.gif, ${PREFIX}_epsilon.gif"
  cat default/domain.cfg default/$POL.cfg $CFGFILE \
      | ../src/maxwell2d_gif epsilon_gif_file=${PREFIX}_epsilon.gif \
      - > ${PREFIX}_$POL.gif
  chmod a+r ${PREFIX}_epsilon.gif
done
