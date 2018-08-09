/*
  Periodicity
 */
#include <stdint.h>
#include <stdlib.h>
#include <math.h>

int pad_box(const double inv_boxsize, const double r_pad, const uint32_t num_pos, 
	    const double *restrict pos, double *restrict pad_pos, 
	    int64_t *restrict pad_idx, const int max_images)
{
  /*
    Scale every position into the unit cube and find the image in [0,1) and
    append it into [1, 1+r_pad*inv_boxsize)
    
   */
  const double pad = r_pad * inv_boxsize; // Scaled r_pad
  int image = 0; // image to add

  double *new_pos = &pad_pos[(size_t)num_pos*3];
  const double *im_start = new_pos;

  for (size_t i=0; i<num_pos;++i, im_start = new_pos)
    {

      // Scale to unit cube
      const double ux = pos[i*3]*inv_boxsize,
	uy = pos[i*3+1]*inv_boxsize,
	uz = pos[i*3+2]*inv_boxsize;
      
      const double px = ux - floor(ux),
	py = uy - floor(uy),
	pz = uz - floor(uz);	

      pad_pos[i*3]   = px;
      pad_pos[i*3+1] = py;
      pad_pos[i*3+2] = pz;
      
      int prev_ims = 0;
      if (px<pad)
	{
	  prev_ims =1;
	  if (image==max_images)
	    return -1; // Not enough space for x images
	  pad_idx[image++] = i;
	  *new_pos++ = px + 1.0;
	  *new_pos++ = py;
	  *new_pos++ = pz;
	}

      if (py<pad)
	{
	  const double *my_ims = im_start;
	  for (int j=prev_ims;j;--j)
	    {
	      if (image==max_images)
		return -1; // Not enough space for y images
	      pad_idx[image++] = i;
	      *new_pos++ = *my_ims++;
	      *new_pos++ = 1.0 + *my_ims++;
	      *new_pos++ = *my_ims++;
	    }

	  if (image==max_images)
	    return -1;

	  pad_idx[image++] = i;
	  *new_pos++ = px;
	  *new_pos++ = py+1.0;
	  *new_pos++ = pz;
	  
	  prev_ims = prev_ims*2+1;
	}
      if (pz<pad)
	{
	  const double *my_ims = im_start;
	  for (int j=prev_ims;j;--j)
	    {
	      if (image==max_images)
		return -1; // Not enough space for z images

	      pad_idx[image++] = i;
	      *new_pos++ = *my_ims++;
	      *new_pos++ = *my_ims++;
	      *new_pos++ = 1 + *my_ims++;
	    }
	  
	  if (image==max_images)
	    return -1;

	  pad_idx[image++] = i;
	  *new_pos++ = px;
	  *new_pos++ = py;
	  *new_pos++ = pz+1;
	}
    }

  return image;
}
