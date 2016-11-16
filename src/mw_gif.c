/* mw_gif.c -- Functions to manipulate gif files

   Copyright (C) 2008-2009 Robin Hogan <r.j.hogan@reading.ac.uk> 

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

#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "gif_lib.h"
#include "maxwell.h"
#include "jet_and_black.h"

#define PROGRAM_NAME	"Maxwell2D"
#define GIF_ASM_NAME   "NETSCAPE2.0"
#define COMMENT_GIF_ASM    "(C) Robin Hogan 2009"

#define JET_SIZE 64
#define HALF_JET_SIZE 32
#define MAX_MAG 10
/*static int AsmGifAnimNumIters = 1, /* Loop just once */
  static int AsmGifAnimNumIters = 0, /* Loop indefinitely */
  AsmGifAnimDelay = 10,
  AsmGifAnimUserWait = FALSE;

static GifFileType *gif_file; 
static ColorMapObject *color_map;
static GifByteType *gif_line = NULL;
static int gif_mag = 1;

/* Initialize animated gif file of name "filename" */
int
mw_gif_init(char *filename, mwDomain *domain)
{
  int gif_file_id;
  int width;
  int height;
  unsigned char ExtStr[3];

  gif_mag = domain->mag;
  if (gif_mag > MAX_MAG) {
    gif_mag = MAX_MAG;
  }
  else if (gif_mag <= 0) {
    gif_mag = 1;
  }

  width = gif_mag*(1 + (domain->mode & MW_MODE_VACUUM))*domain->nx;
  height = gif_mag*domain->ny;

  /*
  ExtStr[0] = AsmGifAnimNumIters % 256;
  ExtStr[1] = AsmGifAnimNumIters / 256;
  ExtStr[2] = 0;
  */
  ExtStr[0] = 1;
  ExtStr[1] = AsmGifAnimNumIters % 256;
  ExtStr[2] = AsmGifAnimNumIters / 256;

  if (!(color_map = MakeMapObject(JET_SIZE, JetPalette))) {
    return MW_FAILURE;
  }

  if (filename) {
    gif_file_id = open(filename,O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
    if (gif_file_id == -1) {
      return MW_FAILURE;
    }
  }
  else {
    gif_file_id = 1; /* stdout */
  }

  EGifSetGifVersion("89a");
  if ((gif_file = EGifOpenFileHandle(gif_file_id)) == NULL) {
    return MW_FAILURE;
  }

  if (EGifPutScreenDesc(gif_file, width, height,
 			color_map->BitsPerPixel,
			0, color_map) == GIF_ERROR) {
    return MW_FAILURE;
  }

  if ((gif_line = malloc(width*sizeof(GifByteType))) == NULL) {
    return MW_FAILURE;
  }

  EGifPutExtensionFirst(gif_file, APPLICATION_EXT_FUNC_CODE,
			strlen(GIF_ASM_NAME), GIF_ASM_NAME);
  EGifPutExtensionLast(gif_file, APPLICATION_EXT_FUNC_CODE,
		       3, ExtStr);
  EGifPutExtension(gif_file, COMMENT_EXT_FUNC_CODE,
		   strlen(COMMENT_GIF_ASM), COMMENT_GIF_ASM);

  return MW_SUCCESS;
}

/* Write a frame of the animated gif file */      
int
mw_gif_write_frame(mwDomain *domain)
{
  real **field, **vac;
  real plot_max, scat_max;
  int width = gif_mag*(1 + (domain->mode & MW_MODE_VACUUM))*domain->nx;
  int height = gif_mag*domain->ny;
  int i, j, k;
  unsigned char ExtStr[4] = { 0x04, 0x00, 0x00, 0xff };

  ExtStr[0] = AsmGifAnimUserWait ? 0x06 : 0x04;
  ExtStr[1] = AsmGifAnimDelay % 256;
  ExtStr[2] = AsmGifAnimDelay / 256;
  ExtStr[3] = 0xff;

  EGifPutExtension(gif_file, GRAPHICS_EXT_FUNC_CODE,
		   4, ExtStr);

  if (EGifPutImageDesc(gif_file, 0, 0, width, height,
		       FALSE, NULL) == GIF_ERROR) {
    return MW_FAILURE;
  }

  if (domain->mode & MW_MODE_EZ) {
    field = domain->Ez;
    plot_max = domain->plot_E_max;
    vac = domain->Ez_vacuum;
    scat_max = domain->plot_E_max*domain->plot_scat_ratio;
  }
  else {
    field = domain->Bz;
    plot_max = domain->plot_B_max;
    vac = domain->Bz_vacuum;
    scat_max = domain->plot_B_max*domain->plot_scat_ratio;
  }

  for (j = domain->ny-1; j >= 0; j--) {
    for (i = 0; i < domain->nx; i++) {
      real value = HALF_JET_SIZE*(1.0+field[j][i]/plot_max);
      if (domain->boundaries[j][i] > 0.0) {
	value = JET_SIZE-1;
      }
      else if (value < 0.0) {
	value = 0;
      }
      else if (value >= JET_SIZE-1) {
	value = JET_SIZE-2;
      }
      for (k = 0; k < gif_mag; k++) {
	gif_line[i*gif_mag + k] = (GifByteType) value;
      }
      if (domain->mode & MW_MODE_VACUUM) {
	value = HALF_JET_SIZE*(1.0+(field[j][i]-vac[j][i])/scat_max);
	if (value < 0.0) {
	  value = 0;
	}
	else if (value >= JET_SIZE-1) {
	  value = JET_SIZE-2;
	}
	for (k = 0; k < gif_mag; k++) {
	  gif_line[(i+domain->nx)*gif_mag+k] = (GifByteType) value;
	}
      }
    }
    for (k = 0; k < gif_mag; k++) {
      if (EGifPutLine(gif_file, gif_line, width) == GIF_ERROR) {
	return MW_FAILURE;
      }
    }  
  }
  return MW_SUCCESS;
}

/* Close the gif file */      
int
mw_gif_close()
{
  if (EGifCloseFile(gif_file) == GIF_ERROR) {
    return MW_FAILURE;
  }
  return MW_SUCCESS;
}

/* Write a gif file containing the dielectric constant field */
int
mw_gif_write_epsilon(char *filename, mwDomain *domain)
{
  GifFileType *file; 
  GifByteType *line = NULL;
  int file_id;
  int width;
  int height;
  int i, j, k;
  int mag = domain->mag;
  if (mag > MAX_MAG) {
    mag = MAX_MAG;
  }
  else if (mag <= 0) {
    mag = 1;
  }

  width = mag*domain->nx;
  height = mag*domain->ny;
 
  MW_CHECK(!(color_map = MakeMapObject(JET_SIZE, JetPalette)));

  if (filename) {
    file_id = open(filename, O_WRONLY|O_CREAT|O_TRUNC, S_IRUSR|S_IWUSR);
    if (file_id == -1) {
      return MW_FAILURE;
    }
  }
  else {
    file_id = 1; /* stdout */
  }

  MW_CHECK((file = EGifOpenFileHandle(file_id)) == NULL);

  MW_CHECK(EGifPutScreenDesc(file, width, height,
			     color_map->BitsPerPixel,
			     0, color_map) == GIF_ERROR);

  MW_CHECK((line = malloc(width*sizeof(GifByteType))) == 0);

  MW_CHECK(EGifPutImageDesc(file, 0, 0, width, height,
			    FALSE, NULL) == GIF_ERROR);

  for (j = domain->ny-1; j >= 0; j--) {
    for (i = 0; i < domain->nx; i++) {
      real value = JET_SIZE*domain->epsilon[j][i]/4.0;
      if (value < 0.0) {
	value = 0;
      }
      else if (value >= JET_SIZE-1) {
	value = JET_SIZE-2;
      }
      for (k = 0; k < mag; k++) {
	line[i*mag + k] = (GifByteType) value;
      }
    }
    for (k = 0; k < mag; k++) {
      MW_CHECK(EGifPutLine(file, line, width) == GIF_ERROR);
    }
  }
  MW_CHECK(EGifCloseFile(file) == GIF_ERROR);

  return MW_SUCCESS;
}
