# Makefile for maxwell2d_nc and maxwell2d_gif

# C compiler
CC = gcc

# Optimization flags
OPTFLAGS = -O3 --fast-math

# Flags to provide to the C compiler
CFLAGS = -I../giflib-4.1.6/lib -Wall -g -s -pipe $(OPTFLAGS) \
	-I/usr/include/netcdf-3 -I/opt/graphics/include

# Object files required by both programs
OBJECTS = mw_alloc.o mw_shape.o mw_step.o mw_print.o mw_start.o \
	mw_frame.o mw_math.o mw_boundaries.o readconfig.o

# Gif-specific object files
GIFOBJECTS = main_gif.o mw_gif.o

# NetCDF-specific object files
NCOBJECTS =  main_nc.o mw_nc.o nctools.o

# Prefix for the two program names
PROGRAM_PREFIX = maxwell2d

# General libraries required
LIBS = -L../giflib-4.1.6/lib/.libs -lm \
	-L../../netcdf-3.6.0-p1/src/libsrc
#-L/opt/graphics/lib -lm \


# The default target will compile both programs
all: $(PROGRAM_PREFIX)_gif $(PROGRAM_PREFIX)_nc

# "make maxwell2d_gif" will compile only the gif version of the program
$(PROGRAM_PREFIX)_gif: $(OBJECTS) $(GIFOBJECTS)
	$(CC) --static -o $(PROGRAM_PREFIX)_gif $(OBJECTS) $(GIFOBJECTS) $(LIBS) -lgif

# "make maxwell2d_nc" will compile only the NetCDF version of the program
$(PROGRAM_PREFIX)_nc: $(OBJECTS) $(NCOBJECTS)
	$(CC) -o $(PROGRAM_PREFIX)_nc $(OBJECTS) $(NCOBJECTS) $(LIBS) -lnetcdf

# Object file dependencies
%.o: %.c *.h
	$(CC) $(CFLAGS) -c $<

# Type "make clean" to remove object files and executables
clean:
	rm -f $(OBJECTS) $(GIFOBJECTS) $(NCOBJECTS) \
		$(PROGRAM_PREFIX)_gif $(PROGRAM_PREFIX)_nc

# Type "make clean-autosaves" to remove Emacs autosave files
clean-autosaves:
	rm -f *~
