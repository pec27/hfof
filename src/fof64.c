/*
  Like fof.c but with cells grouped into 4x4x4 (=64) blocks
 */
#include "fof.h"
#include <math.h>
#include <stdlib.h>

#ifdef DEBUG
  #include <stdio.h>
static int max_stack_usage, pt_cmp, fill_collisions, search_collisions, ngb_found;
#endif

// A walk over the 13 cells within sqrt(3) cell widths, and having idx<me, 
// where idx(i,j,k)=Mi+Nj+k
#define WALK_NGB(M,N) {\
    N, M, 1, M+N, N+1, M-N, N-1, M+1, M-1, M+N-1, M-N-1, M-N+1, M+N+1}

// hashed position of the neighbour (i.e. (ngb*prime) mod N)
#define HASH_NGB(A,P,M) {\
    A[0]*P&M, A[1]*P&M, A[2]*P&M, A[3]*P&M, A[4]*P&M, A[5]*P&M, A[6]*P&M, \
      A[7]*P&M, A[8]*P&M, A[9]*P&M, A[10]*P&M, A[11]*P&M, A[12]*P&M}

static const unsigned int quad_masks[64] = {
  0xdc743210, 0xdc742310, 0xd986310, 0xd986310, 0xdc430721, 0xdc473201, 0xd968301, 0xd963081, 0xdb5721, 0xdb7521, 0xda851, 0xda581, 0xdb7521, 0xdb7251, 0xda851, 0xda851, 0xdc731420, 0xdc743210, 0xd986310, 0xd983160, 0xdc374102, 0xd743210, 0xd86310, 0xd938610, 0xdb5712, 0xd7521, 0xd851, 0xda581, 0xdb7512, 0xdb7521, 0xda851, 0xda851, 0xd420, 0xd420, 0xd60, 0xd60, 0xd402, 0xd420, 0xd60, 0xd60, 0xd2, 0xd2, 0xd, 0xd, 0xd2, 0xd2, 0xd, 0xd, 0xd420, 0xd420, 0xd60, 0xd60, 0xd402, 0xd420, 0xd60, 0xd60, 0xd2, 0xd2, 0xd, 0xd, 0xd2, 0xd2, 0xd, 0xd};

static const unsigned char h_ngb[65] = {
  0, 7, 14, 20, 26, 33, 40, 46, 52, 57, 62, 66, 70, 75, 80, 84, 88, 95, 102, 108, 114, 121, 127, 132, 138, 143, 147, 150, 154, 159, 164, 168, 172, 175, 178, 180, 182, 185, 188, 190, 192, 193, 194, 194, 194, 195, 196, 196, 196, 199, 202, 204, 206, 209, 212, 214, 216, 217, 218, 218, 218, 219, 220, 220, 220};

static const uint64_t ngb_bits[284] = {
  0x37707770777, 0x730077007700, 0x777037700000000, 0x8cc0ccc0ccc, 
  0x7700730000000000, 0xc800cc00cc00, 0xccc08cc00000000, 0xcc00c80000000000L, 
  0x7ff0fff0fff, 0xf700ff00ff00, 0xfff07ff00000000, 0xff00f70000000000L, 
  0x8808880888, 0x800088008800, 0x888008800000000, 0x8800800000000000L, 
  0xeff0fff0fff, 0xfe00ff00ff00, 0xfff0eff00000000, 0xff00fe0000000000L, 
  0x100011001100, 0x111001100000000, 0x1100100000000000, 0xcee0eee0eee, 
  0xec00ee00ee00, 0xeee0cee00000000, 0xee00ec0000000000L, 0x310033003300, 
  0x333013300000000, 0x3300310000000000, 0x377777777777, 0x7777377700000000, 
  0x8ccccccccccc, 0xcccc8ccc00000000L, 0x300070007000, 0x7000300000000000, 
  0x8000c000c000, 0xc000800000000000L, 0x7fffffffffff, 0xffff7fff00000000L, 
  0x7000f000f000, 0x88888888888, 0xf000700000000000L, 0x8888088800000000L, 
  0x80008000, 0x8000000000000000L, 0xefffffffffff, 0xffffefff00000000L, 
  0xe000f000f000, 0xf000e00000000000L, 0x1111011100000000, 0x10001000, 
  0x1000000000000000, 0xceeeeeeeeeee, 0xeeeeceee00000000L, 0x3333133300000000, 
  0xc000e000e000, 0xe000c00000000000L, 0x100030003000, 0x3000100000000000, 
  0x777377777777, 0x7777777300000000, 0xccc8cccccccc, 0xccccccc800000000L, 
  0x7000300000000, 0xc000800000000, 0xfff7ffffffff, 0xfffffff700000000L, 
  0x888088888888, 0xf000700000000, 0x8888888000000000L, 0x8000000000000, 
  0xfffeffffffff, 0xfffffffe00000000L, 0xf000e00000000, 0x1111111000000000, 
  0x1000000000000, 0xeeeceeeeeeee, 0xeeeeeeec00000000L, 0x3333333100000000, 
  0xe000c00000000, 0x3000100000000, 0x773077707770, 0x7770773000000000, 
  0xcc80ccc0ccc0, 0x77003700000000, 0xccc0cc8000000000L, 0xcc008c00000000, 
  0xff70fff0fff0, 0xfff0ff7000000000L, 0xff007f00000000, 0x880088808880, 
  0x8880880000000000L, 0x88000800000000, 0xffe0fff0fff0, 0xfff0ffe000000000L, 
  0xff00ef00000000, 0x1110110000000000, 0x11000100000000, 0xeec0eee0eee0, 
  0xeee0eec000000000L, 0xee00ce00000000, 0x3330331000000000, 0x33001300000000, 
  0x377077707770777, 0x7300770077007700, 0x8cc0ccc0ccc0ccc, 0xc800cc00cc00cc00L, 
  0x377000000000000, 0x7300000000000000, 0x8cc000000000000, 0xc800000000000000L, 
  0x7ff0fff0fff0fff, 0xf700ff00ff00ff00L, 0x7ff000000000000, 0x88088808880888, 
  0xf700000000000000L, 0x8000880088008800L, 0x88000000000000, 0x8000000000000000L, 
  0xeff0fff0fff0fff, 0xfe00ff00ff00ff00L, 0xeff000000000000, 0xfe00000000000000L, 
  0x1000110011001100, 0x11000000000000, 0x1000000000000000, 0xcee0eee0eee0eee, 
  0xec00ee00ee00ee00L, 0x3100330033003300, 0xcee000000000000, 0xec00000000000000L, 
  0x133000000000000, 0x3100000000000000, 0x3777777777777777, 0x8cccccccccccccccL, 
  0x3000700070007000, 0x3777000000000000, 0x8000c000c000c000L, 0x8ccc000000000000L, 
  0x3000000000000000, 0x8000000000000000L, 0x7fffffffffffffff, 0x7000f000f000f000, 
  0x7fff000000000000, 0x888888888888888, 0x7000000000000000, 0x800080008000, 
  0x888000000000000, 0xefffffffffffffffL, 0xe000f000f000f000L, 0xefff000000000000L, 
  0xe000000000000000L, 0x100010001000, 0x111000000000000, 0xceeeeeeeeeeeeeeeL, 
  0xc000e000e000e000L, 0xceee000000000000L, 0x1000300030003000, 0x1333000000000000, 
  0xc000000000000000L, 0x1000000000000000, 0x7773777777777777, 0xccc8ccccccccccccL, 
  0x7773000000000000, 0xccc8000000000000L, 0x3000000000000, 0x8000000000000, 
  0xfff7ffffffffffffL, 0xfff7000000000000L, 0x8880888888888888L, 0x7000000000000, 
  0x8880000000000000L, 0xfffeffffffffffffL, 0xfffe000000000000L, 0xe000000000000, 
  0x1110000000000000, 0xeeeceeeeeeeeeeeeL, 0xeeec000000000000L, 0x3331000000000000, 
  0xc000000000000, 0x1000000000000, 0x7730777077707770, 0xcc80ccc0ccc0ccc0L, 
  0x7730000000000000, 0x37000000000000, 0xcc80000000000000L, 0x8c000000000000, 
  0xff70fff0fff0fff0L, 0xff70000000000000L, 0x8800888088808880L, 0x7f000000000000, 
  0x8800000000000000L, 0x8000000000000, 0xffe0fff0fff0fff0L, 0xffe0000000000000L, 
  0xef000000000000, 0x1100000000000000, 0x1000000000000, 0xeec0eee0eee0eee0L, 
  0xeec0000000000000L, 0xce000000000000, 0x3310000000000000, 0x13000000000000, 
  0x777077707770377, 0x7700770077007300, 0xccc0ccc0ccc08cc, 0xcc00cc00cc00c800L, 
  0xfff0fff0fff07ff, 0xff00ff00ff00f700L, 0x888088808880088, 0x8800880088008000L, 
  0xfff0fff0fff0eff, 0xff00ff00ff00fe00L, 0x1100110011001000, 0xeee0eee0eee0cee, 
  0xee00ee00ee00ec00L, 0x3300330033003100, 0x7777777777773777, 0xcccccccccccc8cccL, 
  0x7000700070003000, 0xc000c000c0008000L, 0xffffffffffff7fffL, 0xf000f000f0007000L, 
  0x8888888888880888L, 0x8000800080000000L, 0xffffffffffffefffL, 0xf000f000f000e000L, 
  0x1000100010000000, 0xeeeeeeeeeeeeceeeL, 0xe000e000e000c000L, 0x3000300030001000, 
  0x7777777777777773, 0xccccccccccccccc8L, 0xfffffffffffffff7L, 0x8888888888888880L, 
  0xfffffffffffffffeL, 0xeeeeeeeeeeeeeeecL, 0x7770777077707730, 0xccc0ccc0ccc0cc80L, 
  0xfff0fff0fff0ff70L, 0x8880888088808800L, 0xfff0fff0fff0ffe0L, 0xeee0eee0eee0eec0L, 
  0x777077703770000, 0x7700770073000000, 0xccc0ccc08cc0000, 0xcc00cc00c8000000L, 
  0xfff0fff07ff0000, 0xff00ff00f7000000L, 0x888088800880000, 0x8800880080000000L, 
  0xfff0fff0eff0000, 0xff00ff00fe000000L, 0x1100110010000000, 0xeee0eee0cee0000, 
  0xee00ee00ec000000L, 0x3300330031000000, 0x7777777737770000, 0xcccccccc8ccc0000L, 
  0x7000700030000000, 0xc000c00080000000L, 0xffffffff7fff0000L, 0xf000f00070000000L, 
  0x8888888808880000L, 0x8000800000000000L, 0xffffffffefff0000L, 0xf000f000e0000000L, 
  0x1000100000000000, 0xeeeeeeeeceee0000L, 0xe000e000c0000000L, 0x3000300010000000, 
  0x7777777777730000, 0xccccccccccc80000L, 0xfffffffffff70000L, 0x8888888888800000L, 
  0xfffffffffffe0000L, 0xeeeeeeeeeeec0000L, 0x7770777077300000, 0xccc0ccc0cc800000L, 
  0xfff0fff0ff700000L, 0x8880888088000000L, 0xfff0fff0ffe00000L, 0xeee0eee0eec00000L};


/* SubCells */
typedef struct {
  uint32_t parent, start; // Parent (disjoint set) and start point per cell
  unsigned int octant:8, rank:8; // 1-8
} SubCell;


void blocks_cells(const double *pos, const int num_pos, const double inv_cell_width, const int ny, const int nx, int64_t *out)
{
  /*
    Find the bucket index (int)(p[0]*inv_cell_width)*nx + (int)(p[1]*inv_cell_width)*ny + (int)(p[2]*inv_cell_width)
    for every point
   */

  const int64_t NY=ny, NX = nx;
  for (size_t i=0;i<num_pos;i++)
    {
      int64_t ix = (int64_t)(pos[i*3]*inv_cell_width),
	iy = (int64_t)(pos[i*3+1]*inv_cell_width),
	iz = (int64_t)(pos[i*3+2]*inv_cell_width);

      out[i] = ((ix>>2)*NX+(iy>>2)*NY + (iz>>2))<<6 | (ix&3)<< 4 | (iy&3)<<2 | (iz&3);

    }

}
static inline unsigned int find_path_compress(const unsigned int idx, SubCell *restrict cells)
{
  /*
    Find root (domain) with path compression (i.e. set all parents to final root 
  */

  unsigned int stack_used=0, stack[16], domain = cells[idx].parent;

  for(unsigned int child=idx;child!=domain;domain = cells[child=domain].parent)
    stack[stack_used++] = child;
  
#ifdef DEBUG
  if (stack_used>max_stack_usage)
    printf("Increasing stack used to %d\n",(max_stack_usage = stack_used));
#endif

  // Set all those paths to have final parent
  while (stack_used)
    cells[stack[--stack_used]].parent = domain;

  return domain;
}

static int connected_pwise(const int cur_size, const int adj_size,
			   const double b2,
			   const int64_t* restrict cur_idx,
			   const int64_t *restrict adj_idx,
			   const double *restrict xyz)
{
  /*
    Check pairwise for any connections, i.e. any point p1 in my
    cell within b of any point p2 in adj cell

    Returns 1 if connected, 0 otherwise
  */

  for (int my_k=cur_size;my_k;--my_k, cur_idx++)
    {
      const double *xyz1 = &xyz[(*cur_idx)*3];
      const int64_t* p2 = adj_idx;      
      for (int adj_k=adj_size;
	   adj_k; --adj_k, p2++)
	{
	  const double *xyz2 = &xyz[(*p2)*3];
	  const double dx = xyz2[0]-xyz1[0],
	    dy = xyz2[1]-xyz1[1],
	    dz = xyz2[2]-xyz1[2];
#ifdef DEBUG		      
	  pt_cmp++; // Count comparisons
#endif
	  if ((dx*dx + dy*dy + dz*dz) > b2)
	    continue;
	  
	  return 1;
	}
    }
  return 0;
}

typedef struct {
  int64_t high_id; // ID of cell
  int32_t idx; // index of cell (when ordered)
  unsigned int num_subcells;
} HashCell;

int fof(const int num_pos, const int N, const int M, const double b, 
	const double *restrict xyz, const int64_t *restrict cell_ids, 
	const int64_t *restrict sort_idx, int32_t *restrict domains,
	const double desired_load)
{
  /*
    Friends of Friends by binning into cells.

    num_pos  - The number of points
    N        - Multiplicative factor such that cell(i,j,k) = M i + N j + k. 
               Usually 3 mod 4 for maximum bitwise independence
    M        - prime bigger than N^2
    b        - The linking length
    xyz      - (num_pos * 3) array of (unsorted) points
    cell_ids - The cell values for each point (s.t. i=x/cell_width etc.)
    sort_idx - An index such that cell_ids[sort_idx] is ascending (ordered points
               by cell)
    domains  - (num_pos,) integer output array for the domains
    desired_load - Load for the hash table (0.6 probably best)
    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */
  
  int num_cells=0, num_big_cells=1;

  const unsigned int walk_ngbs[13] = WALK_NGB(M,N);

  const double b2 = b*b;

  /* Count the number of cells */
  for (int64_t i=1,last_id = cell_ids[sort_idx[num_cells++]];i<num_pos;++i)
    if (cell_ids[sort_idx[i]]!=last_id)
      {
	num_cells++;
	// Check if high bits are different
	if ((cell_ids[sort_idx[i]]^last_id)>>6!=0) 
	  num_big_cells++;

	last_id = cell_ids[sort_idx[i]];
      }


  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*desired_load<num_big_cells) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  const int64_t hsize = TAB_START<<tab_size;
  const unsigned int hprime = HASH_PRIMES[tab_size],
    hmask = hsize-1;

  const unsigned int hash_ngb[13] = HASH_NGB(walk_ngbs, hprime, hmask);

#ifdef DEBUG
  float actual_load = (float)num_big_cells / hsize, expected_collisions= (-1.0 - log(1-actual_load)/actual_load);
  printf("Number of cells %d, %.2f%% of number of positions (%u)\n", num_big_cells, num_big_cells*100.0/num_pos, (unsigned int)num_pos);
  printf("Number of octancts %d, %.2f%% of number of positions (%u)\n", num_cells, num_cells*100.0/num_pos, (unsigned int)num_pos);
  printf("Table size 1024<<%d, prime %u\n", tab_size, hprime);
  printf("Desired load %.1f%%\nActual load %.1f%%\n", 100.0*desired_load, 100.0*actual_load);

  printf("Platform int size %d\nHash cell size %d bytes\n",
	 (int)sizeof(int), (int)sizeof(HashCell));
  max_stack_usage = search_collisions = ngb_found = fill_collisions = pt_cmp = 0;
#endif

  // Hashtable of indices
  HashCell *htable = malloc(hsize * sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  SubCell* cells = malloc(num_cells*sizeof(SubCell));

  if (!cells)
    return -1;

  // Hash table fill bits
  uint32_t *hfill = calloc(hsize>>5, sizeof(uint32_t));
  if (!hfill)
    return -5;

  for (unsigned int cur_start=0, big_cell=0, my_idx=0;
       cur_start<num_pos;big_cell=my_idx)
    {
      const int64_t my_high_id = cell_ids[sort_idx[cur_start]]>>6;
      const unsigned int my_hash = (unsigned int)my_high_id*hprime&hmask;
      //      int my_oct_fill = 0;

      // Loop over octants
      for (int cur_end=cur_start+1;
	   cur_start<num_pos && cell_ids[sort_idx[cur_start]]>>6==my_high_id; // TODO try masks not shifts?
	   ++my_idx, cur_start=cur_end)
	{
	  const int64_t my_id = cell_ids[sort_idx[cur_start]];
	  // Find all points in octant
	  while (cur_end<num_pos && cell_ids[sort_idx[cur_end]]==my_id)
	    cur_end++;

	  // Put in own domain
	  cells[my_idx].parent = my_idx; 
	  cells[my_idx].start = cur_start;
	  unsigned int my64 = cells[my_idx].octant = (unsigned int)my_id & 63;

	  cells[my_idx].rank = 0;

	  const uint64_t *nbits = ngb_bits + (h_ngb[my64]+my64);
	  // Octant done, connect to previous octants (if any)
	  for (unsigned int adj_idx=my_idx;adj_idx!=big_cell;) 
	    {
	      if (!(*nbits>>cells[--adj_idx].octant&1))
		continue;
	      const unsigned int my_root = cells[my_idx].parent,
		adj_root = find_path_compress(adj_idx, cells);
	      
	      if (adj_root!=my_root && 
		  connected_pwise(cur_end-cur_start,
				  cells[adj_idx+1].start - cells[adj_idx].start, 
				  b2, sort_idx+cur_start, sort_idx + cells[adj_idx].start, xyz))
		{
		  // Connect the domains - Disjoint sets union algorithm
		  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
		    cells[my_idx].parent = cells[my_root].parent = adj_root;
		  else if (cells[cells[adj_root].parent = my_root].rank == cells[adj_root].rank) // Add adj to me
		    cells[my_root].rank++; // Arbitrarily choose, in tests add adj to me slightly faster
		}
	    }

	  // Link to adjacent cells
	  for (unsigned int mu=quad_masks[my64], adj;(adj=mu&0xF)!=13;mu>>=4, nbits++)
	    for (unsigned int adj_hash=(my_hash-hash_ngb[adj])&hmask;
		 hfill[adj_hash>>5]>>(adj_hash&31)&1; adj_hash=(adj_hash+1)&hmask)
	      if (htable[adj_hash].high_id==my_high_id-walk_ngbs[adj]) // My neighbour or collision?
		{
#ifdef DEBUG
		  ngb_found++;
#endif
		  // First draft just try to connect every octant pairwise (backwards)
		  for (unsigned int adj_idx = htable[adj_hash].idx+htable[adj_hash].num_subcells;
		       adj_idx!=htable[adj_hash].idx; )
		    {
		      if (!(nbits[1]>>cells[--adj_idx].octant&1))
			continue;

		      const unsigned int adj_root = find_path_compress(adj_idx, cells),
			my_root = cells[my_idx].parent;

		      if (adj_root!=my_root &&
			  connected_pwise(cur_end-cur_start,
					  cells[adj_idx+1].start - cells[adj_idx].start, 
					  b2, sort_idx+cur_start, sort_idx + cells[adj_idx].start, xyz))
			{

			  // Connect the domains - Disjoint sets union algorithm
			  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
			    cells[my_idx].parent = cells[my_root].parent = adj_root;
			  else if (cells[cells[adj_root].parent = my_root].rank == cells[adj_root].rank) // Add adj to me
			    cells[my_root].rank++; // Arbitrarily choose, in tests add adj to me slightly faster
			}
		    }
		  break; // Found neighbour so done
		}
#ifdef DEBUG
	      else
		search_collisions++; // Wrong cell, move to next
#endif
	}
       
      // Add to hash table      
      unsigned int cur = my_hash;
#ifdef DEBUG
      while (hfill[cur>>5]&1U<<(cur&31))
	{
	  fill_collisions++;
	  cur = (cur+1)&hmask;
	}
     
#else
      while (hfill[cur>>5]&(1U<<(cur&31)))
	cur = (cur+1)&hmask;
#endif
      
      hfill[cur>>5] |= 1U<<(cur&31);
      
      htable[cur].idx = big_cell;
      htable[cur].high_id = my_high_id;
      htable[cur].num_subcells = my_idx - big_cell;
    }
  
      
  free(hfill);
  free(htable);

  // Renumber the domains from 0...num_doms-1
  // A bit hack-y, but use the cell starts to hold the domains
  int num_doms = 0;
  for (size_t i=0;i<num_cells;++i)
    if (cells[i].parent==i)
      cells[i].start = num_doms++; 

  unsigned int domain = find_path_compress(0, cells);
  domains[sort_idx[0]] = cells[domain].start;
  int64_t cur_cell_id = cell_ids[sort_idx[0]];
  for (size_t i=0,j=1;j<num_pos; ++j)
    {
      /* Find root of this domain (path compression of disjoint sets) */
      if (cell_ids[sort_idx[j]]!=cur_cell_id)
	{
	  cur_cell_id = cell_ids[sort_idx[j]];
	  domain = find_path_compress(++i, cells);
	}
      domains[sort_idx[j]] = cells[domain].start; // see hack above
    }

#ifdef DEBUG
  int max_rank = cells[0].rank;
  for (size_t i=1;i<num_cells;++i)
    if (cells[i].rank>max_rank)
	max_rank = cells[i].rank;

  float frac_found = ngb_found/ (13.0*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f average distance comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, search_collisions, num_doms, frac_found, ngb_found);
  printf("%.2f%% fill collisions (c.f. %.2f%% for perfect hashing)\n", (float)fill_collisions*100.0/num_big_cells, 100.0*expected_collisions);

  printf("Max stack usage %d, max rank %d, point comparisons %d\n",max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif

  //  free(ranks);
  free(cells);

  return num_doms;
}
