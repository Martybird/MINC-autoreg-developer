/* ----------------------------- MNI Header -----------------------------------
@NAME       : blur.c
@INPUT      : 
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: functions to blur a 2D or 3D image
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
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

@CREATED    : Dec 1, 2008 CY
---------------------------------------------------------------------------- */

#include <math.h>
#include <float.h>
#include <volume_io.h>
#include "blur.h"

// Indices for row, column and slice
#define SLC_IDX 0
#define ROW_IDX 1
#define COL_IDX 2

#define KERN_UNDEF    0
#define KERN_GAUSSIAN 1
#define KERN_RECT     2

#define PI2 6.28318530717959

/* ----------------------------- MNI Header -----------------------------------
@NAME       : muli_vects
@INPUT      : s1 - a zero offset array containing real,imag,real,imag
              s2 - a zero offset array containing real,imag,real,imag
              n  - the number of complex pairs to be multiplied
@OUTPUT     : r  - the result of the multiplication, a zero offset array 
                   containing real,imag,real,imag
@RETURNS    : nothing
@DESCRIPTION: 
              for c = a*b, all real:

                      c.r=a.r*b.r-a.i*b.i;a
                      c.i=a.i*b.r+a.r*b.i;
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void muli_vects(float *r, float *s1, float *s2, int n)
{
  int i;
  float r1,i1,r2,i2;

  r++; s1++; s2++; /* all arrays start at r[1],s1[1] and s2[1], where the real
  part is r[1] and the imag in r[2], and so on... */

  if (r!=s1 && r!=s2) {                /* if separate arrays */
    for (i=0; i< n; ++i) { 
      *r = (*(s1) * *(s2))   - (*(s1+1) * *(s2+1)); 
      r++;
      *r = (*(s1+1) * *(s2)) + (*(s1) * *(s2+1)); 
      r++;
      s1++;s1++;
      s2++;s2++;
    } 
  }
  else {                        /* if one array=result array */
    for (i=0; i< n; ++i) { 
      r1= *(s1); i1= *(s1+1);
      r2= *(s2); i2= *(s2+1);
      *r = (r1 * r2)   - (i1 * i2); 
      r++;
      *r = (i1 * r2) + (r1 * i2);
      r++;
      s1++;s1++;
      s2++;s2++;
    } 
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : next_power_of_two
@INPUT      : x - integer number (positive)
@OUTPUT     : 
@RETURNS    : int - the next higher power of two
@DESCRIPTION: the routine returns the smallest number n > x, such that n = 2^?
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* find the next highest power of 2, greater than or equal  */
/* to x, and return 2^n, when n is the required power.      */
/************************************************************/
int next_power_of_two(int x)
{
  int 
      n, power_of_two;
  
  power_of_two = 1;
  n = 0;
  
  while (power_of_two<x && n<32) {
    power_of_two *= 2;
    n++;
  }
  
  return(power_of_two);
}


void fft1(float *signal, int numpoints, int direction)
{
  int 
      n, m, mmax, 
  i, j, istep;
  double 
        sin_a,
  wp_r,wp_i,
  w_r,w_i,
  angle;
  float 
      temp,
  temp_real,
  temp_imag;
                                /* scramble the entries into
  bit reversed order */
  n = numpoints << 1;
  j = 1;
  for (i = 1; i < n; i += 2) {
    if (j > i) {                /* swap entries i and j */
      temp=signal[j];   signal[j]=signal[i];     signal[i]=temp;
      temp=signal[j+1]; signal[j+1]=signal[i+1]; signal[i+1]=temp;
    }

    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;

  }

                                /* calculate the butterflies, but
  leave normalization to calling
  program if this is an inverse
  xform.  */
  mmax = 2;
  while (n > mmax) {
    angle = PI2/(direction*mmax); /* set up trig recurrance */
    sin_a = sin(0.5*angle);        
    wp_r = -2.0*sin_a*sin_a;        
    wp_i = sin(angle);

    istep = mmax<<1;                /* increment on i      */
    w_r = 1.0;                        /* initial sin and cos */
    w_i = 0.0;

    for (m = 1; m<mmax; m+=2) {        
      for (i = m; i<=n; i+=istep) {
        j = i+mmax;
        temp_real    = w_r*signal[j]   - w_i*signal[j+1]; /* mult */
        temp_imag    = w_r*signal[j+1] + w_i*signal[j]; 
        signal[j]    = signal[i]       - temp_real;       /* sub */
        signal[j+1]  = signal[i+1]     - temp_imag;
        signal[i]   += temp_real;                         /* add */
        signal[i+1] += temp_imag;
      }
      sin_a = w_r;
      w_r = sin_a*wp_r - w_i*wp_i   + w_r; /* trig recurrence  */
      w_i = w_i*wp_r   + sin_a*wp_i + w_i;
    }
    mmax = istep;
  }
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : normal_dist
@INPUT      : c    - height of gaussian
              fwhm - full wifth half max of gaussian
              mu   - center of gaussian
              x    - value of x
@OUTPUT     : 
@RETURNS    : value of gaussian evaluated at x
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
float normal_dist(float c, float fwhm, float mu, float x)
{
  float sigma,t1,t2,t3,f;
  
  sigma = fwhm/2.35482;
  
  if (sigma==0) {
    if (x==mu)
      f = c / (sqrt(2*PI) * sigma);      /* !!!!!!!!!!! */
    else
      f = 0;
  }
  else {
    t1 = c / (sqrt(2*PI) * sigma);
    t2 = (x-mu)*(x-mu)/(2*sigma*sigma);
    t3 = exp(-t2);
    f = t1*t3;
  }
  
  return(f);
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : rect_dist
@INPUT      : c    - height of rect function
              fwhm - width of rect function
              mu   - center of rect function
              x    - value of x
@OUTPUT     : 
@RETURNS    : value of rect function evaluated at x
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
/* return the value of the normal dist at x, given c,sigma  */
/* and mu ----   all in mm                                  */
/************************************************************/
float rect_dist(float c, float fwhm, float mu, float x)
{
  
  float t;
  
  t = x-mu;
  
  /*  t = t<0 ? -t : t ; */
  
  if ( t >= -fwhm/2  && t < fwhm/2 ) {
    return(  (float) (c/fwhm) );
  }
  else
    return ( (float) 0.0 );
  
}


/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_kernel_FT
@INPUT      : kern - a zero offset array containing real,imag,real,imag
                     in which will be stored the kernel for the dirivitive
              size - the number of complex numbers in the kernel array
@OUTPUT     : 
@RETURNS    : 
@DESCRIPTION: 
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
/************************************************************/
void make_kernel_FT(float *kern, int size, float vsize)
{
  
  int 
      kindex,k;
  float
      factor,
  f_sample_size;
  
  (void)memset(kern,(int)0,(size_t)((2*size+1)*sizeof(float)));
  
  f_sample_size = 1.0/(vsize*size);
  factor = -2.0 * PI * f_sample_size;
  
  
  for ( k = -size/2; k<size/2; ++k) {
    kindex = ((k + size) % size)*2 +1; 
    kern[kindex+1] = factor*k;
  }
  
}

/* ----------------------------- MNI Header -----------------------------------
@NAME       : make_kernel
@INPUT      : kern - a zero offset array containing real,imag,real,imag
                     in which will be stored the Gaussian kernel for convolution
              vsize- the size (in mm) of the sample along the kern array
              fwhm - full-width-half-maximum of gaussian (in mm)
              size - the number of complex numbers in the kernel array
@OUTPUT     : kern - the Gaussian kernel used for convolution
@RETURNS    : nothing
@DESCRIPTION: 

   note that kern (the convolution kernel) goes into the array so that the
   peak is at kern[1] (remember unit offset), with the
   positive half of the kernel running from kern[1] to kern[n/2] and the
   negative half running from kern[array_size_pow2 - n/2] to kern[array_size_pow2]:

  -                            size / 2
   \                              V                                 /
  --|-----------------------------+--------------------------------|-
     \_/                                                        \_/ 
  ^                                                                 ^
  0                                                            size ^  

            this is not right? ---^
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : Wed Jun 23 09:04:34 EST 1993 Louis Collins
@MODIFIED   : 
---------------------------------------------------------------------------- */
void make_kernel(float *kern, float vsize, float fwhm, int size, int type)
{

  int kindex,k;
  float c,r,max;

  (void)memset(kern,(int)0,(size_t)((2*size+1)*sizeof(float)));
  
  for ( k = -size/2; k<size/2; ++k) {

    kindex = ((k + size) % size)*2 +1;


    switch (type) {
      case KERN_GAUSSIAN: kern[kindex] = normal_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
      case KERN_RECT:     kern[kindex] = rect_dist(1.0*vsize,fwhm,0.0,(float)(vsize*k)); break;
      default: 
      {
        (void) fprintf (stderr,"Illegal kernel type = %d\n",type);
        (void) fprintf (stderr,"Impossible error in make_kernel(), line %d of %s\n",
        __LINE__,__FILE__);
        k = size/2;
      }
    }
  }

  if (type == KERN_RECT) {
    max = -1e37;
    for(k=1; k<=size; k++)
      if (kern[k]>max) max = kern[k];
    if (max != 0.0) {
      for(k=1; k<=size; k++)
        kern[k] = kern[k] / max;
      c = (int)( fwhm/vsize )==(int)(fwhm/vsize+0.5) ? (int)( fwhm/vsize )+1 : (int)( fwhm/vsize )+2;
      for(k=1; k<=size; k++)
        kern[k] = kern[k] / (fwhm/vsize + 1);
      r = 0.0;
      for(k=1; k<=size; k++)
        r +=kern[k];
      kern[ (int)(c/2 + 0.5) + 1] += (1.0 - r)/2;
      kern[ size - (int)(c/2 + 0.5) + 1] += (1.0 - r)/2;
      
    }
  }

}


/*
For different values of ndim, the image is blurred in different dimensions.
ndim == 1, image is blurred in the slice direction
ndim == 2, image is blurred within the slice
ndim == 3, image is blurred in all dimensions
*/


void blur_volume(VIO_Volume data, double fwhmx, double fwhmy, double fwhmz, 
                       int ndim, int kernel_type)
{ 
 
  float 
    *fdata,                        /* floating point storage for blurred volume */
    *f_ptr,                        /* pointer to fdata */
    tmp,
    *dat_vector,                /* temp storage of original row, col or slice vect. */
    *dat_vecto2,                /* storage of result of dat_vector*kern             */
    *kern;                        /* convolution kernel                               */
                                /*  place it back into data->voxels                 */

  VIO_Real
    lowest_val,
    max_val, 
    min_val;
    
  int                                
    total_voxels,                
    vector_size_data,                /* original size of row, col or slice vector        */
    kernel_size_data,                /* original size of kernel vector                   */
    array_size_pow2,                /* actual size of vector/kernel data used in FFT    */
                                /* routines - needs to be a power of two            */
    array_size;
  int   
    data_offset;                /* offset required to place original data (size n)  */
                                /*  into array (size m=2^n) so that data is centered*/

  
  register int 
    slice_limit,
    row,col,slice,                /* counters to access original data                 */
    vindex;                        /* counter to access vector and vecto2              */

  int 
    slice_size,                        /* size of each data step - in bytes              */
    row_size, col_size;
                
  char
    full_outfilename[256];        /* name of output file */

  VIO_Status 
    status;
  
  int
    sizes[3];                        /* number of rows, cols and slices */
  VIO_Real
    steps[3];                        /* size of voxel step from center to center in x,y,z */

  

  /*---------------------------------------------------------------------------------*/
  /*             start by setting up the raw data.                                   */
  /*---------------------------------------------------------------------------------*/

  get_volume_sizes(data, sizes);          /* rows,cols,slices */
  get_volume_separations(data, steps);
  
  slice_size = sizes[COL_IDX] * sizes[ROW_IDX];    /* sizeof one slice  */
  col_size   = sizes[ROW_IDX];               /* sizeof one column */
  row_size   = sizes[COL_IDX];               /* sizeof one row    */
  
  total_voxels = sizes[COL_IDX]*sizes[ROW_IDX]*sizes[SLC_IDX];

  //printf("VIO_X = %d, VIO_Y = %d, VIO_Z = %d  ", VIO_X, VIO_Y, VIO_Z);
  //printf("sizes[0] = %d, sizes[1] = %d, sizes[2] = %d\n", sizes[0],  sizes[1],  sizes[2]);
  //printf("ndim = %d, kernel type  = %d.\n", ndim, kernel_type);
  
  ALLOC(fdata, total_voxels);

  lowest_val = get_volume_voxel_min(data);
    

  get_volume_real_range(data, &min_val, &max_val);

  max_val = -FLT_MAX;
  min_val = FLT_MAX;

  f_ptr = fdata;
  for(slice=0; slice<sizes[SLC_IDX]; slice++)
    for(row=0; row<sizes[ROW_IDX]; row++)
      for(col=0; col<sizes[COL_IDX]; col++) {
        GET_VOXEL_3D( tmp, data, slice, row, col);
        //GET_VOXEL_3D( tmp, data, col, row, slice);

        if (tmp <= lowest_val)
          tmp = lowest_val;

        *f_ptr = CONVERT_VOXEL_TO_VALUE(data, tmp);
        if (max_val < *f_ptr) max_val = *f_ptr;
        if (min_val > *f_ptr) min_val = *f_ptr;
        f_ptr++;
      }



  /* note data is stored by rows (along x, i.e. increasing col index), then by cols (along y, ie.
     increasing row index) then slices (along z, i.e. increasing slice index) */
  
  /*--------------------------------------------------------------------------------------*/
  /*                get ready to start up the transformation.                             */
  /*--------------------------------------------------------------------------------------*/
    
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  vector_size_data = sizes[COL_IDX]; 
  kernel_size_data = (int)(((4*fwhmx)/ABS(steps[COL_IDX])) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                    /* complex number for FFT, and the plus 1*/
                                    /* is for the zero offset FFT routine    */
  ALLOC(dat_vector, array_size);
  ALLOC(dat_vecto2, array_size);
  ALLOC(kern,       array_size);

  
  /*--------------------------------------------------------------------------------------*/
  /*                start with rows                                                       */
  /*--------------------------------------------------------------------------------------*/
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,(float)(ABS(steps[COL_IDX])),fwhmx,array_size_pow2,kernel_type);
  fft1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[COL_IDX])/2;
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  slice_limit = 0;
  switch (ndim) {
    case 1: slice_limit = 0; break;
    case 2: slice_limit = sizes[SLC_IDX]; break;
    case 3: slice_limit = sizes[SLC_IDX]; break;
  }


  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    for (row = 0; row < sizes[ROW_IDX]; row++) {           /* for each row   */
      
      f_ptr = fdata + slice*slice_size + row*sizes[COL_IDX];
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (col=0; col< sizes[COL_IDX]; col++) {        /* extract the row */
        dat_vector[1 +2*(col+data_offset)  ] = *f_ptr++;
      }
      
      fft1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      fft1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + row*sizes[COL_IDX];
      for (col=0; col< sizes[COL_IDX]; col++) {        /* put the row back */
        
        vindex = 1 + 2*(col+data_offset);
        
        *f_ptr++ = dat_vecto2[vindex]/array_size_pow2;
      }
    }
  }
  
  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);

 
  /*--------------------------------------------------------------------------------------*/
  /*                 now do cols                                                          */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  f_ptr = fdata;
  
  vector_size_data = sizes[ROW_IDX];
  kernel_size_data = (int)(((4*fwhmy)/(ABS(steps[ROW_IDX]))) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                     /* complex number for FFT, and the plus 1*/
                                     /* is for the zero offset FFT routine    */
  
  ALLOC(dat_vector, array_size);
  ALLOC(dat_vecto2, array_size);
  ALLOC(kern,       array_size);
  
  /*    1st calculate kern array for gaussian kernel*/
  
  make_kernel(kern,(float)(ABS(steps[ROW_IDX])),fwhmy,array_size_pow2,kernel_type);
  fft1(kern,array_size_pow2,1);
  
  /*    calculate offset for original data to be placed in vector            */
  
  data_offset = (array_size_pow2-sizes[ROW_IDX])/2;
  
  /*    2nd now convolve this kernel with the rows of the dataset            */
  
  switch (ndim) {
    case 1: slice_limit = 0; break;
    case 2: slice_limit = sizes[SLC_IDX]; break;
    case 3: slice_limit = sizes[SLC_IDX]; break;
  }


  for (slice = 0; slice < slice_limit; slice++) {      /* for each slice */
    
    for (col = 0; col < sizes[COL_IDX]; col++) {           /* for each col   */
      
      f_ptr = fdata + slice*slice_size + col;
      
      memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
      
      for (row=0; row< sizes[ROW_IDX]; row++) {        /* extract the col */
        dat_vector[1 +2*(row+data_offset) ] = *f_ptr;
        f_ptr += row_size;
      }
      
      
      fft1(dat_vector,array_size_pow2,1);
      muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);
      fft1(dat_vecto2,array_size_pow2,-1);
      
      f_ptr = fdata + slice*slice_size + col;
      for (row=0; row< sizes[ROW_IDX]; row++) {        /* put the col back */
        
        vindex = 1 + 2*(row+data_offset);
        
        *f_ptr = dat_vecto2[vindex]/array_size_pow2;
        
        f_ptr += row_size;
      }
    }
  }
  
  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);
  
  
  
  /*--------------------------------------------------------------------------------------*/
  /*                 now do slices                                                        */
  /*--------------------------------------------------------------------------------------*/
  
  /*-----------------------------------------------------------------------------*/
  /*             determine   size of data structures needed                      */
  
  
  f_ptr = fdata;
  
  vector_size_data = sizes[SLC_IDX];
  kernel_size_data = (int)(((4*fwhmz)/(ABS(steps[SLC_IDX]))) + 0.5);
  
  if (kernel_size_data > MAX(vector_size_data,256))
    kernel_size_data =  MAX(vector_size_data,256);
  
  /*             array_size_pow2 will hold the size of the arrays for FFT convolution,
                 remember that ffts require arrays 2^n in length                          */
  
  array_size_pow2  = next_power_of_two(vector_size_data+kernel_size_data+1);
  array_size = 2*array_size_pow2+1;  /* allocate 2*, since each point is a    */
                                     /* complex number for FFT, and the plus 1*/
                                     /* is for the zero offset FFT routine    */

  ALLOC(dat_vector, array_size); 
  ALLOC(dat_vecto2, array_size); 
  ALLOC(kern,       array_size); 
  
  max_val = -FLT_MAX;
  min_val = FLT_MAX;
    
  if (ndim==1 || ndim==3) {
    
    /*    1st calculate kern array for gaussian kernel*/
    
    make_kernel(kern,(float)(ABS(steps[SLC_IDX])),fwhmz,array_size_pow2,kernel_type);
    fft1(kern,array_size_pow2,1);
    
    /*    calculate offset for original data to be placed in vector            */
    
    data_offset = (array_size_pow2-sizes[SLC_IDX])/2;
    
    /*    2nd now convolve this kernel with the slices of the dataset            */
    
    for (col = 0; col < sizes[COL_IDX]; col++) {      /* for each column */
      
      for (row = 0; row < sizes[ROW_IDX]; row++) {           /* for each row   */
        
        f_ptr = fdata + col*col_size + row;
        
        memset(dat_vector,0,(2*array_size_pow2+1)*sizeof(float));
        
        for (slice=0; slice< sizes[SLC_IDX]; slice++) {        /* extract the slice vector */
          dat_vector[1 +2*(slice+data_offset) ] = *f_ptr;
          f_ptr += slice_size;
        }
        
        fft1(dat_vector,array_size_pow2,1);
        muli_vects(dat_vecto2,dat_vector,kern,array_size_pow2);

        fft1(dat_vecto2,array_size_pow2,-1);
        
        f_ptr = fdata + col*col_size + row;
        
        for (slice=0; slice< sizes[SLC_IDX]; slice++) {        /* put the vector back */
          
          vindex = 1 + 2*(slice+data_offset);
          
          *f_ptr = dat_vecto2[vindex]/array_size_pow2;
          
          if (max_val<*f_ptr) max_val = *f_ptr;
          if (min_val>*f_ptr) min_val = *f_ptr;
          
          f_ptr += slice_size;
        }
      }
    }
  }  /* if ndim */
  else {

    for (slice = 0; slice < sizes[SLC_IDX]; slice++) {      /* for each slice */
      for (col = 0; col < sizes[COL_IDX]; col++) {             /* for each column */
        for (row = 0; row < sizes[ROW_IDX]; row++) {           /* for each row   */
          if (max_val<*f_ptr) max_val = *f_ptr;
          if (min_val>*f_ptr) min_val = *f_ptr;
          f_ptr++;
        }
      }
    }
  }

  FREE(dat_vector);
  FREE(dat_vecto2);
  FREE(kern);

  
/* set up the correct info to copy the data back out in mnc */

  f_ptr = fdata;
  
  set_volume_real_range(data, min_val, max_val);

  for(slice=0; slice<sizes[SLC_IDX]; slice++)
    for(row=0; row<sizes[ROW_IDX]; row++)
      for(col=0; col<sizes[COL_IDX]; col++) {
        tmp = CONVERT_VALUE_TO_VOXEL(data, *f_ptr);
        SET_VOXEL_3D( data, slice, row, col, tmp);
        //SET_VOXEL_3D( data, col, row, slice, tmp);
        f_ptr++;
      }

  FREE(fdata);

  return;
}

