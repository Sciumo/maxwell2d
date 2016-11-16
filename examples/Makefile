# Simple Makefile to create (and remove) gif and NetCDF files using
# the create_gif.sh and create_nc.sh scripts

all: gifs ncs

gifs:
	create_gif.sh *.cfg

ncs:
	create_nc.sh *.cfg

clean:
	rm -f *.nc *.gif *.png

clean-autosaves:
	rm -f *~
