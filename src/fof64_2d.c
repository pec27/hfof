/* Like fof64 but in 2-d */
#include "fof.h"
#include <math.h>
#include <stdlib.h>

#ifdef DEBUG
  #include <stdio.h>
static int max_stack_usage, pt_cmp, fill_collisions, search_collisions, ngb_found;
#endif

// Index for starting point of walk in ngb_bits
static const unsigned char ngb_start[65] = {
0, 3, 6, 7, 8, 9, 10, 12, 14, 17, 19, 20, 21, 22, 23, 24, 26, 27, 28, 28, 28, 28, 28, 28, 28, 29, 30, 30, 30, 30, 30, 30, 30, 31, 32, 32, 32, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34, 34, 34, 35, 36, 36, 36, 36, 36, 36, 36, 37, 38, 38, 38, 38, 38, 38, 38};

/* quad_walks - for each (64) subcells, integer with each hex digit (0-3) 
                denotes one of the neighbouring blocks. 13 (0xD) is the stop signal */
// TODO store with ngb_start ?
const uint16_t quad_walks[64] = {
0xd310, 0xd310, 0xd0, 0xd0, 0xd0, 0xd0, 0xd20, 0xd20, 0xd301, 0xd10, 0xd0, 0xd0, 0xd0, 0xd0, 0xd0, 0xd20, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd1, 0xd1, 0xd, 0xd, 0xd, 0xd, 0xd, 0xd};

/*
  Ragged array of
  ngb_mask[0,...N],

  where ngb_mask[..] - each 8x8 cell in ngb block I am within sqrt(2) of.
*/

// TODO this could be reduced, since there are only 8x2+6x2 = 28 cells on a block that any other block is adjacent to (when going forward in raster order) 
static const uint64_t ngb_bits[38] = {
0x703000000000000, 0x80c0c0, 0xc080000000000000, 0xf07000000000000, 0x8080, 0x8000000000000000, 0x1f0e000000000000, 0x3e1c000000000000, 0x7c38000000000000, 0xf870000000000000, 0xf0e0000000000000, 0x100000000000000, 0xe0c0000000000000, 0x301000000000000, 0x80c0c0c0, 0x300000000000000, 0x8000000000000000, 0x700000000000000, 0x808080, 0xe00000000000000, 0x1c00000000000000, 0x3800000000000000, 0x7000000000000000, 0xe000000000000000, 0xc000000000000000, 0x100000000000000, 0x80c0c0c080, 0x80808000, 0x80c0c0c08000, 0x8080800000, 0x80c0c0c0800000, 0x808080000000, 0x80c0c0c080000000, 0x80808000000000, 0xc0c0c08000000000, 0x8080800000000000, 0xc0c0800000000000, 0x8080000000000000};

// All the 8x8 cells in my block I am within sqrt(2) of
static const uint64_t self_bits[64] = {
0x30707, 0x70f0f, 0xe1f1f, 0x1c3e3e, 0x387c7c, 0x70f8f8, 0xe0f0f0, 0xc0e0e0, 0x3070707, 0x70f0f0f, 0xe1f1f1f, 0x1c3e3e3e, 0x387c7c7c, 0x70f8f8f8, 0xe0f0f0f0, 0xc0e0e0e0, 0x307070703, 0x70f0f0f07, 0xe1f1f1f0e, 0x1c3e3e3e1c, 0x387c7c7c38, 0x70f8f8f870, 0xe0f0f0f0e0, 0xc0e0e0e0c0, 0x30707070300, 0x70f0f0f0700, 0xe1f1f1f0e00, 0x1c3e3e3e1c00, 0x387c7c7c3800, 0x70f8f8f87000, 0xe0f0f0f0e000, 0xc0e0e0e0c000, 0x3070707030000, 0x70f0f0f070000, 0xe1f1f1f0e0000, 0x1c3e3e3e1c0000, 0x387c7c7c380000, 0x70f8f8f8700000, 0xe0f0f0f0e00000, 0xc0e0e0e0c00000, 0x307070703000000, 0x70f0f0f07000000, 0xe1f1f1f0e000000, 0x1c3e3e3e1c000000, 0x387c7c7c38000000, 0x70f8f8f870000000, 0xe0f0f0f0e0000000, 0xc0e0e0e0c0000000, 0x707070300000000, 0xf0f0f0700000000, 0x1f1f1f0e00000000, 0x3e3e3e1c00000000, 0x7c7c7c3800000000, 0xf8f8f87000000000, 0xf0f0f0e000000000, 0xe0e0e0c000000000, 0x707030000000000, 0xf0f070000000000, 0x1f1f0e0000000000, 0x3e3e1c0000000000, 0x7c7c380000000000, 0xf8f8700000000000, 0xf0f0e00000000000, 0xe0e0c00000000000};

/* SubCells */
typedef struct {
  uint32_t parent, start; // Parent (disjoint set) and start point per cell
  unsigned int subcell:8, rank:8; // 1-8
} SubCell;


void get_min_max_2d(const double *pos, const uint32_t num_pos, double *restrict out)
{
  /*
    Find the minimum and maximum (faster than numpy min & max calls)

    pos     - (num_pos,2) doubles of x,y
    num_pos - number of positions
    out     - 4 doubles for result of (min x,y max x,y)
   */

  if (!num_pos)
    return;

  double min_x = pos[0], max_x = min_x,
    min_y = pos[1], max_y = min_y;
  
  for (size_t i=1;i<num_pos;i++)
    {
      if (pos[i*2]<min_x)
	min_x = pos[i*2];
      if (pos[i*2]>max_x)
	max_x = pos[i*2];

      if (pos[i*2+1]<min_y)
	min_y = pos[i*2+1];
      if (pos[i*2+1]>max_y)
	max_y = pos[i*2+1];
    }

  out[0] = min_x;
  out[1] = min_y;
  out[2] = max_x;
  out[3] = max_y;
}

void blocks_cells_2d(const double min_x, const double min_y, 
		  const double *restrict pos, const uint32_t N, 
		  const double inv_cell_width, const int64_t Px, 
		  int64_t *restrict out)
{
  /*
    Convert (x,y) -> to int64_t key := block<<6 | subcell 
    for N positions

    i.e. blocks of 8x8 cells

    where
      block := (ix>>3)*Px + (iy>>3)
      cell := (ix&7)<<3 | (iy&7) # (i.e. 0-63)
    
      ix := (int64_t)((x - min_x)*inv_cell_width)
      iy := ...


    Args:  
    min_x  - minimum in x-direction
    min_y  - minimum in y-direction
    pos    - (N,2) double array of xyz
    N      - number of positions
    inv_cell_width - inverse of subcell width
    nx     - prime 
    out    - (N,) int64_t array for output

    returns void
   */


  // Pre-multiplied top_left corner of box
  const double x_begin = -min_x*inv_cell_width, 
    y_begin = -min_y*inv_cell_width; 

  for (size_t i=0;i<N;i++)
    {
      int64_t ix = (int64_t)(pos[i*2]*inv_cell_width + x_begin),
	iy = (int64_t)(pos[i*2+1]*inv_cell_width + y_begin);

      out[i] = ((ix>>3)*Px + (iy>>3))<<6 | 
	(ix&7)<< 3 | (iy&7);

    }

}

static inline unsigned int find_path_compress(const unsigned int idx, SubCell *restrict cells)
{
  /*
    Find root (domain) with path compression (i.e. set all parents to final root)
  */

  unsigned int stack_used=0, stack[16], domain = cells[idx].parent;

  for(unsigned int child=idx;
      child != domain;
      domain = cells[child=domain].parent)
    {
      stack[stack_used++] = child;
    }
  
#ifdef DEBUG
  if (stack_used>max_stack_usage)
    printf("Increasing stack used to %d\n",(max_stack_usage = stack_used));
#endif

  // Set all those paths to have final parent
  while (stack_used)
    cells[stack[--stack_used]].parent = domain;

  return domain;
}

static inline int connected_pwise(const int64_t *restrict cur_stop, const double b2,
				  const int64_t* restrict cur_idx,
				  const int64_t *restrict adj_idx,
				  const double *restrict xy,
				  const int64_t *restrict cells)
{
  /*
    Check pairwise for any connections, i.e. any point p1 in my
    cell within b of any point p2 in adj cell

    Returns 1 if connected, 0 otherwise
  */
  const int64_t adj_cell = cells[*adj_idx];
  do
    {
      const double *xy2 = &xy[(*adj_idx++)*2];
      
      for (const int64_t* p1 = cur_idx; p1!=cur_stop; p1++)
	{
	  const double *xy1 = &xy[(*p1)*2];
	  
	  const double dx = xy2[0]-xy1[0],
	    dy = xy2[1]-xy1[1];

#ifdef DEBUG		      
	  pt_cmp++; // Count comparisons
#endif
	  if ((dx*dx + dy*dy) > b2)
	    continue;
	  
	  return 1;
	}
    } while (cells[*adj_idx]==adj_cell);
  return 0;
}

typedef struct {
  int64_t block; // ID of cell
  uint32_t idx; // index of cell (when ordered)
  unsigned int free_subcells:6, fill:1;
} HashBlock;

static void loop_cells(const unsigned int num_pos, const int64_t delta_ngb[], 
		       const unsigned int delta_hash[], 
		       const int64_t *restrict cell_ids, 
		       const int64_t *restrict sort_idx, 
		       const unsigned int hprime, const unsigned int hmask, 
		       const double b2, const double *restrict xy, 
		       SubCell *restrict cells, HashBlock *restrict htable)
{
  /*
    Sequentially insert every point (+cell & block) into a hash-table a link to 
    neighbours.
    
    num_pos   - number of xyz coords
    delta_ngb - Delta to index of neighbouring blocks (4 with idx<me)
    delta_hash- Delta to hash of neighboues (i.e. delta_ngb[i]*P mod N)
    cell_ids  - cell idx for every point
    sort_idx  - indices to sort points s.t. cell_ids[sort_idx[...]] are ordered
    hprime    - hash prime
    hmask     - hash mask (2<<L - 1)
    b2        - squared linking length
    xy        - positions (unsorted)
    cells     - SubCell data (to fill)
    htable    - zerod hash table of blocks
   */

  for (unsigned int i=0, my_cell=0; i<num_pos;)
    {
      const int64_t my_block = cell_ids[sort_idx[i]]>>6;
      const unsigned int my_hash = (unsigned int)my_block*hprime&hmask,
	first_cell_in_block = my_cell;

      // Loop over subcells
      do {
	const int64_t my_id = cell_ids[sort_idx[i]];
	// Put in own domain
	const unsigned int my_start = cells[cells[my_cell].parent = my_cell].start = i,
	  my64 = cells[my_cell].subcell = (unsigned int)my_id & 63;
	
	cells[my_cell].rank = 0;
	
	// Find all points in subcell
	while (++i<num_pos && cell_ids[sort_idx[i]]==my_id);
	
	const uint64_t *connection_mask = ngb_bits + ngb_start[my64],
	  self_mask = self_bits[my64];
	
	// Connect to subcells in same block (if any)
	for (unsigned int adj_idx=my_cell;(adj_idx--)!=first_cell_in_block;) 
	  if (self_mask>>cells[adj_idx].subcell&1) // |p_me - p_adj|<b possible
	    {
	      const unsigned int my_root = cells[my_cell].parent,
		adj_root = find_path_compress(adj_idx, cells);
	      
	      if (adj_root!=my_root && 
		  connected_pwise(sort_idx+i, b2, sort_idx + my_start, 
				  sort_idx + cells[adj_idx].start, xy, cell_ids))
		{
		  // Connect the domains - Disjoint sets union algorithm
		  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
		    cells[my_cell].parent = cells[my_root].parent = adj_root;
		  else if (cells[cells[adj_root].parent = my_root].rank == cells[adj_root].rank) // Add adj to me
		    cells[my_root].rank++; // Arbitrarily choose, in tests add adj to me slightly faster
		}
	    }
	// Finished connecting cells in same block
	
	// Connect cells in adjacent blocks
	for (unsigned int walk_blocks=quad_walks[my64], adj=walk_blocks&0xF;
	     adj!=13;adj=(walk_blocks>>=4)&0xF, connection_mask++)
	  // Find block (or none) in hashtable
	  for (unsigned int adj_hash=(my_hash-delta_hash[adj])&hmask;
	       htable[adj_hash].fill; adj_hash=(adj_hash+1)&hmask)
	    if (htable[adj_hash].block==my_block-delta_ngb[adj]) // My neighbour or collision?
	      {
#ifdef DEBUG
		ngb_found++;
#endif
		// Go over every subcell in adj
		for (unsigned int adj_idx = htable[adj_hash].idx, ctr=64-htable[adj_hash].free_subcells;
		     ctr--; ++adj_idx)
		  if (connection_mask[0]>>cells[adj_idx].subcell&1) // |p_me - p_adj|<b possible
		    {
		      
		      const unsigned int adj_root = find_path_compress(adj_idx, cells),
			my_root = cells[my_cell].parent;
		    
		      if (adj_root!=my_root &&
			  connected_pwise(sort_idx+i, b2, sort_idx + my_start, 
					  sort_idx + cells[adj_idx].start, xy, cell_ids))
			{
			  // Connect the domains - Disjoint sets union algorithm
			  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
			    cells[my_cell].parent = cells[my_root].parent = adj_root;
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
	++my_cell;
      } while (i<num_pos && cell_ids[sort_idx[i]]>>6==my_block); // TODO try masks not shifts?
      
      // Add to hash table      
      unsigned int cur = my_hash;
#ifdef DEBUG
      while (htable[cur].fill)
	{
	  fill_collisions++;
	  cur = (cur+1)&hmask;
	}
     
#else
      while (htable[cur].fill)
	cur = (cur+1)&hmask;
#endif
      
      htable[cur].fill=1;
      
      htable[cur].idx = first_cell_in_block;
      htable[cur].block = my_block;
      htable[cur].free_subcells = 64-(my_cell - first_cell_in_block);
    }
}

static void loop_cells_periodic(const unsigned int num_pos, const int64_t delta_ngb[], 
				const unsigned int delta_hash[], 
				const int64_t *restrict cell_ids, const int64_t *restrict sort_idx, 
				const unsigned int hprime, const unsigned int hmask, const double b2, const double *restrict xy, const int64_t* pad_idx, const int num_orig, 
				SubCell *restrict cells, HashBlock *restrict htable)
{
  /*
    As for loop_cells except for each point check it is not an image

    num_orig - First num_orig points in the xyz array are real, subsequent are images
    pad_idx  - (num_pos - num_orig) array of indices of real point per image.
   */

  for (unsigned int i=0, my_cell=0; i<num_pos;)
    {
      const int64_t my_block = cell_ids[sort_idx[i]]>>6;
      const unsigned int my_hash = (unsigned int)my_block*hprime&hmask,
	first_cell_in_block = my_cell;

      // Loop over subcells
      do {
	const int64_t my_id = cell_ids[sort_idx[i]];
	// Put in own domain
	const unsigned int my_start = cells[cells[my_cell].parent = my_cell].start = i,
	  my64 = cells[my_cell].subcell = (unsigned int)my_id & 63;
	
	cells[my_cell].rank = 0;

	// Check for any false images and link to cell of original
	do {
	  const int64_t p1 = sort_idx[i];
	  
	  if (cell_ids[p1]!=my_id)
	    break; // diff cell
	  
	  if (p1<num_orig) // Original (not an image)
	    continue; 
	  
	  const int64_t orig_cell = cell_ids[pad_idx[p1-num_orig]];
	  const int64_t orig_block = orig_cell>>6;
	  const unsigned int orig_subcell = orig_cell&63;
	  
	  // Find cell in hash table (must be there, since earlier insertion)
	  unsigned int adj_hash=(unsigned int)orig_block*hprime&hmask;
	  while (htable[adj_hash].block!=orig_block) // Found or collision?
	    adj_hash=(adj_hash+1) & hmask;
	  
	  unsigned int adj_cell = htable[adj_hash].idx;
	  while (cells[adj_cell].subcell != orig_subcell)
	    adj_cell++;
	  
	  // Link me (if not already)
	  // Find root of this domain (path compression of disjoint sets)
	  const unsigned int adj_root = find_path_compress(adj_cell, cells),
	    my_root = cells[my_cell].parent;	
	  
	  if (adj_root==my_root) // Already linked
	    continue;
	  
	  // Connected (since my image in both cells)
	  // Connect the domains - Disjoint sets union algorithm
	  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
	    cells[my_cell].parent = cells[my_root].parent = adj_root;
	  else if (cells[cells[adj_root].parent = my_root].rank == cells[adj_root].rank) // Add adj to me
	    cells[my_root].rank++; // Arbitrarily choose, in tests add adj to me slightly faster
	} while ((++i)<num_pos);

	const uint64_t *connection_mask = ngb_bits + ngb_start[my64],
	  self_mask = self_bits[my64];
	
	// Connect to subcells in same block (if any)
	for (unsigned int adj_idx=my_cell;(adj_idx--)!=first_cell_in_block;) 
	  if (self_mask>>cells[adj_idx].subcell&1) // |p_me - p_adj|<b possible
	    {
	      const unsigned int my_root = cells[my_cell].parent,
		adj_root = find_path_compress(adj_idx, cells);
	      
	      if (adj_root!=my_root && 
		  connected_pwise(sort_idx+i, b2, sort_idx + my_start, 
				  sort_idx + cells[adj_idx].start, xy, cell_ids))
		{
		  // Connect the domains - Disjoint sets union algorithm
		  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
		    cells[my_cell].parent = cells[my_root].parent = adj_root;
		  else if (cells[cells[adj_root].parent = my_root].rank == cells[adj_root].rank) // Add adj to me
		    cells[my_root].rank++; // Arbitrarily choose, in tests add adj to me slightly faster
		}
	    }
	// Finished connecting cells in same block
	
	// Connect cells in adjacent blocks
	for (unsigned int walk_blocks=quad_walks[my64], adj=walk_blocks&0xF;
	     adj!=13;adj=(walk_blocks>>=4)&0xF, connection_mask++)
	  // Find block (or none) in hashtable
	  for (unsigned int adj_hash=(my_hash-delta_hash[adj])&hmask;
	       htable[adj_hash].fill; adj_hash=(adj_hash+1)&hmask)
	    if (htable[adj_hash].block==my_block-delta_ngb[adj]) // My neighbour or collision?
	      {
#ifdef DEBUG
		ngb_found++;
#endif
		// Go over every subcell in adj
		for (unsigned int adj_idx = htable[adj_hash].idx, ctr=64-htable[adj_hash].free_subcells;
		     ctr--; ++adj_idx)
		  if (connection_mask[0]>>cells[adj_idx].subcell&1) // |p_me - p_adj|<b possible
		    {
		      
		      const unsigned int adj_root = find_path_compress(adj_idx, cells),
			my_root = cells[my_cell].parent;
		    
		      if (adj_root!=my_root &&
			  connected_pwise(sort_idx+i, b2, sort_idx + my_start, 
					  sort_idx + cells[adj_idx].start, xy, cell_ids))
			{
			  
			  // Connect the domains - Disjoint sets union algorithm
			  if (cells[my_root].rank < cells[adj_root].rank) // Add me to adj
			    cells[my_cell].parent = cells[my_root].parent = adj_root;
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
	++my_cell;
      } while (i<num_pos && cell_ids[sort_idx[i]]>>6==my_block); // TODO try masks not shifts?
      
      // Add to hash table      
      unsigned int cur = my_hash;
#ifdef DEBUG
      while (htable[cur].fill)
	{
	  fill_collisions++;
	  cur = (cur+1)&hmask;
	}
     
#else
      while (htable[cur].fill)
	cur = (cur+1)&hmask;
#endif
      
      htable[cur].fill=1;
      
      htable[cur].idx = first_cell_in_block;
      htable[cur].block = my_block;
      htable[cur].free_subcells = 64 - (my_cell - first_cell_in_block);
    }
}

int fof64_2d(const uint32_t num_pos, const int N, const uint32_t num_orig, const double b, 
	     const double *restrict xy, const int64_t *restrict cell_ids, 
	     const int64_t *restrict sort_idx, const int64_t* restrict pad_idx,
	     int32_t *restrict domains, const double desired_load)
{
  /*
    Friends of Friends by binning into cells.

    num_pos  - The number of points
    N        - Multiplicative factor such that cell(i,j) = N i + j
               Usually prime for maximum bitwise independence
    num_orig - Number of points that are not 'false' images
    b        - The linking length
    xy2      - (num_pos * 2) array of (unsorted) points
    cell_ids - The cell values for each point (s.t. i=x/cell_width etc.)
    sort_idx - An index such that cell_ids[sort_idx] is ascending (ordered points
               by cell)
    [pad_idx]- Optional indices of the images in the original array (for linking)
    domains  - (num_pos,) integer output array for the domains
    desired_load - Load for the hash table (0.6 probably best)
    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */
  
  unsigned int num_cells=0, num_blocks=1;

  const double b2 = b*b;

  /* Count the number of cells */
  for (int64_t i=1,last_id = cell_ids[sort_idx[num_cells++]];i<num_pos;++i)
    if (cell_ids[sort_idx[i]]!=last_id)
      {
	num_cells++;
	// Check if high bits (i.e. block) is different
	if ((cell_ids[sort_idx[i]]^last_id)>>6!=0) 
	  num_blocks++;

	last_id = cell_ids[sort_idx[i]];
      }

  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*desired_load<num_blocks) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  const int64_t hsize = TAB_START<<tab_size;
  const unsigned int hprime = HASH_PRIMES[tab_size],
    hmask = hsize-1;

#ifdef DEBUG
  float actual_load = (float)num_blocks / hsize, expected_collisions= (-1.0 - log(1-actual_load)/actual_load);
  printf("Number of blocks %d, %.2f%% of number of positions (%u)\n", num_blocks, num_blocks*100.0/num_pos, (unsigned int)num_pos);
  printf("Number of cells %d, %.2f%% of number of positions (%u)\n", num_cells, num_cells*100.0/num_pos, (unsigned int)num_pos);
  printf("Table size 1024<<%d, prime %u\n", tab_size, hprime);
  printf("Desired load %.1f%%\nActual load %.1f%%\n", 100.0*desired_load, 100.0*actual_load);

  printf("Platform int size %d\nHash cell size %d bytes\n",
	 (int)sizeof(int), (int)sizeof(HashBlock));
  max_stack_usage = search_collisions = ngb_found = fill_collisions = pt_cmp = 0;
#endif

  SubCell* cells = malloc(num_cells * sizeof(SubCell));
  if (!cells)
    return -1; // Not enough mem.

  // A walk over the <=4 blocks within sqrt(2) cell widths, and having idx<me, 
  // where idx(i,j)=Ni+j. Store (hash(idx),idx) for each neighbour,
  // where hash(idx)= (ngb * prime) mod N

  const int64_t delta_ngb[4] = {N, 1, N-1, N+1};

  const unsigned int delta_hash[4] = {N*hprime&hmask, hprime&hmask, (N-1)*hprime&hmask, (N+1)*hprime&hmask};  

  // Hashtable of indices
  HashBlock *htable = calloc(hsize, sizeof(HashBlock));

  if (!htable)
    {
      free(cells);
      return -2; // not enough mem.
    }

  if (num_orig==num_pos) // Non-periodic case
    loop_cells(num_pos, delta_ngb, delta_hash, cell_ids, sort_idx, hprime, 
	       hmask, b2, xy, cells, htable);
  else
    loop_cells_periodic(num_pos, delta_ngb, delta_hash, cell_ids, sort_idx, 
			hprime, hmask, b2, xy, pad_idx, num_orig, cells, htable);

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
      const unsigned int p1 = sort_idx[j];
      if (cell_ids[p1]!=cur_cell_id)
	{
	  cur_cell_id = cell_ids[sort_idx[j]];
	  domain = find_path_compress(++i, cells);
	}
      if (p1<num_orig) // Ignore images
	domains[p1] = cells[domain].start; // see hack above
    }

#ifdef DEBUG
  int max_rank = cells[0].rank;
  for (size_t i=1;i<num_cells;++i)
    if (cells[i].rank>max_rank)
	max_rank = cells[i].rank;

  float frac_found = ngb_found/ (4.0*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f average distance comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, search_collisions, num_doms, frac_found, ngb_found);
  printf("%.2f%% fill collisions (c.f. %.2f%% for perfect hashing)\n", (float)fill_collisions*100.0/num_blocks, 100.0*expected_collisions);

  printf("Max stack usage %d, max rank %d, point comparisons %d\n",max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif

  free(cells);

  return num_doms;
}

