#!/bin/bash

#for epsilon in circle1 circle3 circle10 circle20 gradient edge ripple
#for epsilon in circle1 circle3 circle10
for epsilon in  microwave_oven
do
  echo $epsilon
  pol=z
  if [ "" ]
  then
      cat domain.cfg $pol.cfg $epsilon.cfg \
	  | ../src/maxwell_gif epsilon_gif_file=${epsilon}_epsilon.gif \
	  - > ${epsilon}_$pol.gif
  else
      cat domain.cfg $pol.cfg $epsilon.cfg \
	  | ../src/maxwell_nc nc_file=${epsilon}.nc - 
  fi
#  pol=xy
#  cat domain.cfg $pol.cfg $epsilon.cfg \
#      | ./maxwell epsilon_gif_file=${epsilon}_epsilon.gif - > ${epsilon}_$pol.gif
done
 