/*
  Utilities for testing (e.g. find morton index)
 */
#include <stdint.h>
#include <stdlib.h>

void get_morton_idx(const double *pos, const uint32_t num_pos, const double inv_cell_width, int64_t *restrict out)
{
  /*
    Morton idx for pos in [0,1)^3
   */

  for (size_t i=0;i<num_pos;i++)
    {
      // Integer indices in x,y,z
      int64_t x = pos[i*3]*inv_cell_width,
	y = pos[i*3+1]*inv_cell_width,
	z = pos[i*3+2]*inv_cell_width;

      // convert to morton (works up to 21 bits)
      x = (x | x << 32) & 0x1f00000000ffff;
      x = (x | x << 16) & 0x1f0000ff0000ff;
      x = (x | x << 8) & 0x100f00f00f00f00f;
      x = (x | x << 4) & 0x10c30c30c30c30c3;
      x = (x | x << 2) & 0x1249249249249249;

      y = (y | y << 32) & 0xf00000000ffff;
      y = (y | y << 16) & 0x1f0000ff0000ff;
      y = (y | y << 8) & 0x100f00f00f00f00f;
      y = (y | y << 4) & 0x10c30c30c30c30c3;
      y = (y | y << 2) & 0x1249249249249249;

      z = (z | z << 32) & 0xf00000000ffff;
      z = (z | z << 16) & 0x1f0000ff0000ff;
      z = (z | z << 8) & 0x100f00f00f00f00f;
      z = (z | z << 4) & 0x10c30c30c30c30c3;
      z = (z | z << 2) & 0x1249249249249249;

      out[i] = x | (y << 1) | (z << 2); // Z-order
    }
}
