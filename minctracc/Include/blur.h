/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur.h
@DESCRIPTION: File containing declarations of routines for 
              prototypes for blur_image
@COPYRIGHT  :
              Copyright 2008 Charles Yan, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.

@CREATED    : Mon Dec. 1 16:54:02 EST 2008 (CY) 
@MODIFIED   : 
---------------------------------------------------------------------------- */
#include <volume_io.h>
void blur_volume(VIO_Volume data, double fwhmx, double fwhmy, double fwhmz, 
                       int ndim, int kernel_type);