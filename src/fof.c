/*
  Friends of Friends via hash grid
 */

#include "fof.h"
#include <math.h>
#include <stdlib.h>


#ifdef DEBUG
  #include <stdio.h>
static int max_stack_usage, pt_cmp, fill_collisions, search_collisions, ngb_found;
#endif

// A walk over the 58 cells within sqrt(3) cell widths, and having idx<me, 
// where idx(i,j,k)=Mi+Nj+k
#define WALK_NGB_SIZE 58
#define WALK_NGB(M,N) {\
    M, 1, N, M+N, M-N, M-1, M+1, N+1, N-1, M-N+1, M+N+1, M-N-1, M+N-1, 2*M, 2*N,\
      2, 2*M+N, M+2*N, 2*M-N, M-2*N, N+2, N-2, M+2, 2*M+1, M-2, 2*N-1, 2*N+1, \
      2*M-1, 2*M-N-1, 2*M+N-1, M+2*N-1, M-2*N+1, M+N-2, 2*M-N+1, 2*M+N+1, \
      M+2*N+1, M-N+2, M+N+2, M-N-2, M-2*N-1, 2*M+2*N, 2*N-2, 2*M-2*N, 2*M-2, \
      2*M+2, 2*N+2, 2*M+N-2, 2*M-N-2, M+2*N+2, 2*M-2*N-1, 2*M+2*N-1, 2*M-2*N+1, \
      2*M+2*N+1, M-2*N+2, 2*M-N+2, 2*M+N+2, M+2*N-2, M-2*N-2}

// hashed position of the neighbour (i.e. (ngb*prime) mod N)
#define HASH_NGB(A,P,M) {\
    A[0]*P&M, A[1]*P&M, A[2]*P&M, A[3]*P&M, A[4]*P&M, A[5]*P&M, A[6]*P&M, \
      A[7]*P&M, A[8]*P&M, A[9]*P&M, A[10]*P&M, A[11]*P&M, A[12]*P&M, A[13]*P&M,\
      A[14]*P&M, A[15]*P&M, A[16]*P&M, A[17]*P&M, A[18]*P&M, A[19]*P&M, \
      A[20]*P&M, A[21]*P&M, A[22]*P&M, A[23]*P&M, A[24]*P&M, A[25]*P&M, \
      A[26]*P&M, A[27]*P&M, A[28]*P&M, A[29]*P&M, A[30]*P&M, A[31]*P&M, \
      A[32]*P&M, A[33]*P&M, A[34]*P&M, A[35]*P&M, A[36]*P&M, A[37]*P&M, \
      A[38]*P&M, A[39]*P&M, A[40]*P&M, A[41]*P&M, A[42]*P&M, A[43]*P&M, \
      A[44]*P&M, A[45]*P&M, A[46]*P&M, A[47]*P&M, A[48]*P&M, A[49]*P&M, \
      A[50]*P&M, A[51]*P&M, A[52]*P&M, A[53]*P&M, A[54]*P&M, A[55]*P&M, \
      A[56]*P&M, A[57]*P&M}

// Magic number that specifies how full the hash-table should be to be most efficient
// Hash tables perform very badly beyond about 0.7
// We use a lower load since spatial hashing uses lots of speculative searches
#define DESIRED_LOAD 0.3 


typedef struct {
  uint32_t parent, start; // Parent (disjoint set) and start point per cell
} CellParentStart;

typedef struct {
  int64_t cell; // ID of cell
  uint32_t ancestor, // Hint for root
    idx; // index of cell (when ordered)
} HashCell;

void find_lattice(const double *pos, const int num_pos, const double inv_cell_width, const int ny, const int nx, int64_t *out)
{
  /*
    Find the bucket index (int)(p[0]*inv_cell_width)*nx + (int)(p[1]*inv_cell_width)*ny + (int)(p[2]*inv_cell_width)
    for every point
   */

  const int64_t NY=ny, NX = nx;
  for (size_t i=0;i<num_pos;i++)
    out[i] = (int64_t)(pos[i*3]*inv_cell_width)*NX + (int64_t)(pos[i*3+1]*inv_cell_width)*NY + (int64_t)(pos[i*3+2]*inv_cell_width);

}
static inline unsigned int find_path_compress(const unsigned int idx, CellParentStart *restrict cells)
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
int fof_link_cells(const int num_pos, const int N, const int M, const double b, 
		   const double *restrict xyz, const int64_t *restrict cell_ids, 
		   const int64_t *restrict sort_idx, int32_t *restrict domains)
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

    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */

  int num_cells=0;

  const int walk_ngbs[WALK_NGB_SIZE] = WALK_NGB(M,N);

  const double b2 = b*b;

  if (!num_pos)
    return 0; // No points => no domains!

  /* Count the number of cells */
  for (int64_t i=1,last_id = cell_ids[sort_idx[num_cells++]];i<num_pos;++i)
    if (cell_ids[sort_idx[i]]!=last_id)
      {
	num_cells++;
	last_id = cell_ids[sort_idx[i]];
      }


  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  const int64_t hsize = TAB_START<<tab_size;
  const unsigned int hprime = HASH_PRIMES[tab_size],
    hmask = hsize-1;

  const unsigned int hash_ngb[WALK_NGB_SIZE] = HASH_NGB(walk_ngbs, hprime, hmask);

#ifdef DEBUG
  float actual_load = (float)num_cells / hsize, expected_collisions= (-1.0 - log(1-actual_load)/actual_load);
  printf("Number of cells %d, %.2f%% of number of positions (%u)\n", num_cells, num_cells*100.0/num_pos, (unsigned int)num_pos);
  printf("Table size 1024<<%d, prime %u\n", tab_size, hprime);
  printf("Desired load %.1f%%\nActual load %.1f%%\n", 100.0*DESIRED_LOAD, 100.0*actual_load);

  printf("Platform int size %d\nHash cell size %d bytes\n",
	 (int)sizeof(int), (int)sizeof(HashCell));
  max_stack_usage = search_collisions = ngb_found = fill_collisions = pt_cmp = 0;
#endif

  // Hashtable of indices
  HashCell *htable = malloc(hsize * sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  CellParentStart* cells = malloc(num_cells*sizeof(CellParentStart));

  if (!cells)
    return -1;

  uint8_t *ranks = malloc(num_cells*sizeof(uint8_t));
  if (!ranks)
    return -4; // not enough mem.

  // Hash table fill bits
  uint32_t *hfill = calloc(hsize>>5, sizeof(uint32_t));
  if (!hfill)
    return -5;

  // Loop over every (filled) cell
  for (size_t i=0, cur_end=0;cur_end<num_pos; ++i)
    {
      // Put this in a new domain
      ranks[cells[i].parent=i] = 0; // Own domain

      const size_t cur_start = cells[i].start = cur_end;
      const int64_t cur_cell_id = cell_ids[sort_idx[cur_end++]]; // ID for this cell

      unsigned int cur = (unsigned int)cur_cell_id*hprime&hmask;

      /* Start and end of point data for this cell */
      while (cur_end<num_pos && cell_ids[sort_idx[cur_end]]==cur_cell_id)
	cur_end++;

      // Loop over adjacent cells (within sqrt(3)*cell_width)
      for (int adj=0;adj<WALK_NGB_SIZE;++adj)
	for (unsigned int j=(cur - hash_ngb[adj])&hmask; // Look up in hash table (search for cell or 0)
	     hfill[j>>5]>>(j&31)&1;j=(j+1)&hmask)
	  if (htable[j].cell==cur_cell_id - walk_ngbs[adj]) 
	    {
	      // Found my cell, Find root of this domain (path compression of disjoint sets)
#ifdef DEBUG
	      ngb_found++;
#endif
	      const unsigned int my_root = cells[i].parent;
	      
	      if (my_root==htable[j].ancestor) // Some false negatives (i.e. just an ancestor), never false positives.
		break;
	      
	      const int adj_idx = htable[j].idx;
	      const unsigned int adj_root = htable[j].ancestor = find_path_compress(adj_idx, cells);
	      
	      /* Other domain to check for connection? */
	      if (adj_root!=my_root && 
		  connected_pwise(cur_end-cur_start,
				  cells[adj_idx+1].start - cells[adj_idx].start, 
				  b2, sort_idx+cur_start, sort_idx + cells[adj_idx].start, xyz))
		{
		  
		  // Connect the domains
		  // Disjoint sets union algorithm
		  if (ranks[my_root] < ranks[adj_root]) // Add me to adj
		    cells[i].parent = cells[my_root].parent = adj_root;
		  else if (ranks[htable[j].ancestor = cells[adj_root].parent = my_root] == ranks[adj_root]) // Add adj to me
		    ranks[my_root]++; // Arbitrarily choose, in tests add adj to me slightly faster
		}
	      break;
	    }
#ifdef DEBUG
	  else
	    search_collisions++; // Wrong cell, move to next
#endif
  
      // Insert my cell into hash-table at the next free spot
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
      
      htable[cur].idx = i;
      htable[cur].cell = cur_cell_id;
      htable[cur].ancestor = cells[i].parent;
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
  int max_rank = ranks[0];
  for (size_t i=1;i<num_cells;++i)
    if (ranks[i]>max_rank)
	max_rank = ranks[i];

  float frac_found = ngb_found/ (WALK_NGB_SIZE*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f average distance comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, search_collisions, num_doms, frac_found, ngb_found);
  printf("%.2f%% fill collisions (c.f. %.2f%% for perfect hashing)\n", (float)fill_collisions*100.0/num_cells, 100.0*expected_collisions);

  printf("Max stack usage %d, max rank %d, point comparisons %d\n",max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif

  free(ranks);
  free(cells);

  return num_doms;
}

int fof_periodic(const int num_pos, const int N, const int M, const int num_orig, const double b, 
		 const double *restrict xyz, const int64_t *restrict cell_ids, 
		 const int64_t *restrict sort_idx, const int64_t* restrict pad_idx,
		 int32_t *restrict domains)
{
  /*
    Friends of Friends by binning into cells.

    num_pos  - The number of points
    N        - Multiplicative factors such that cell(i,j,k) = M i + N j + k. 
               Usually 3 mod 4 for maximum bitwise independence
    M        - prime bigger than N^2
    num_orig - Number of points that are not 'false' images
    b        - The linking length
    xyz      - (num_pos * 3) array of (unsorted) points
    cell_ids - The cell values for each point (s.t. i=x/cell_width etc.)
    sort_idx - An index such that cell_ids[sort_idx] is ascending (ordered points
               by cell)
    pad_idx  - the indices of the images in the original array (for linking)
    domains  - (num_pos,) integer output array for the domains

    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */

  int num_cells=1;

  const int walk_ngbs[WALK_NGB_SIZE] = WALK_NGB(M,N);

  const double b2 = b*b;

  if (!num_pos)
    return 0; // No points no domains;

  /* Count the number of cells */
  for (int64_t i=1,last_id = cell_ids[sort_idx[num_cells++]];i<num_pos;++i)
    if (cell_ids[sort_idx[i]]!=last_id)
      {
	num_cells++;
	last_id = cell_ids[sort_idx[i]];
      }

  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  const int64_t hsize = TAB_START<<tab_size;
  const unsigned int hprime = HASH_PRIMES[tab_size], hmask = hsize-1;

  const unsigned int hash_ngb[WALK_NGB_SIZE] = HASH_NGB(walk_ngbs, hprime, hmask);

#ifdef DEBUG
  float actual_load = (float)num_cells / hsize, expected_collisions= (-1.0 - log(1-actual_load)/actual_load);
  printf("Number of cells %d, %.2f%% of number of positions (%u)\n", num_cells, num_cells*100.0/num_pos, (unsigned int)num_pos);
  printf("Table size 1024<<%d, prime %u\n", tab_size, hprime);
  printf("Desired load %.1f%%\nActual load %.1f%%\n", 100.0*DESIRED_LOAD, 100.0*actual_load);

  printf("Platform int size %d\nHash cell size %d bytes\n",
	 (int)sizeof(int), (int)sizeof(HashCell));
  max_stack_usage = search_collisions = ngb_found = fill_collisions = pt_cmp = 0;
#endif

  // Hashtable of indices
  HashCell *htable = malloc(hsize * sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  CellParentStart* cells = malloc(num_cells*sizeof(CellParentStart));
  if (!cells)
    return -1;

  uint8_t *ranks = malloc(num_cells*sizeof(uint8_t));
  if (!ranks)
    return -4; // not enough mem.

  // Hash table fill bits
  uint32_t *hfill= calloc(hsize>>5, sizeof(uint32_t));
  if (!hfill)
    return -5;



  // Loop over every (filled) cell
  for (size_t i=0, cur_end=0;i<num_cells;++i)
    {
      // Put this in a new domain
      ranks[cells[i].parent = i] = 0; // Own domain

      const size_t cur_start =  cells[i].start = cur_end;
      const int64_t cur_cell_id = cell_ids[sort_idx[cur_start]]; // ID for this cell

      unsigned int cur = (unsigned int)cur_cell_id*hprime&hmask; // desired spot in hash table

      // Check for any false images and link to cell of original
      do {
	const size_t p1 = sort_idx[cur_end];
	
	if (cell_ids[p1]!=cur_cell_id)
	  break; // diff cell
	
	if (p1<num_orig) // Original (not an image)
	  continue; 

	const int64_t orig_cell = cell_ids[pad_idx[p1-num_orig]];
	
	// Find cell in hash table (must be there)
	unsigned int j = (unsigned int)orig_cell*hprime;
	while (htable[j&=hmask].cell!=orig_cell)
	  j++;
	
	// Link me (if not already)
	// Find root of this domain (path compression of disjoint sets)
	const unsigned int adj_root = find_path_compress(htable[j].idx, cells),
	  my_root = cells[i].parent;	

	if (adj_root==my_root) // Already linked
	  continue;
	
	// Connected (since my image in both cells)
	// Disjoint sets union algorithm
	if (ranks[my_root] < ranks[adj_root]) // Add me to adj
	  cells[i].parent = cells[my_root].parent = adj_root;
	else if (ranks[htable[j].ancestor = cells[adj_root].parent = my_root] == ranks[adj_root]) // Add adj to me
	  ranks[my_root]++; // Arbitrarily choose, in tests add adj to me slightly faster

      } while ((++cur_end)<num_pos);

      // Loop over adjacent cells (within sqrt(3)*cell_width)
      for (int adj=0;adj<WALK_NGB_SIZE;++adj)
	for (unsigned int j=(cur - hash_ngb[adj])&hmask; // Look up in hash table (search for cell or 0)
	     hfill[j>>5]>>(j&31)&1;j=(j+1)&hmask)
	  if (htable[j].cell==cur_cell_id - walk_ngbs[adj]) 
	    {
	      // Found my cell, Find root of this domain (path compression of disjoint sets)
#ifdef DEBUG
	      ngb_found++;
#endif
	      const unsigned int my_root = cells[i].parent;
	      
	      if (my_root==htable[j].ancestor) // Some false negatives (i.e. just an ancestor), never false positives.
		break;
	      
	      const int adj_idx = htable[j].idx;
	      const unsigned int adj_root = htable[j].ancestor = find_path_compress(adj_idx, cells);
	      
	      /* Other domain to check for connection? */
	      if (adj_root!=my_root && 
		  connected_pwise(cur_end-cur_start,
				  cells[adj_idx+1].start - cells[adj_idx].start, 
				  b2, sort_idx+cur_start, sort_idx + cells[adj_idx].start, xyz))
		{
		  
		  // Connect the domains
		  // Disjoint sets union algorithm
		  if (ranks[my_root] < ranks[adj_root]) // Add me to adj
		    cells[i].parent = cells[my_root].parent = adj_root;
		  else if (ranks[htable[j].ancestor = cells[adj_root].parent = my_root] == ranks[adj_root]) // Add adj to me
		    ranks[my_root]++; // Arbitrarily choose, in tests add adj to me slightly faster
		}
	      break;
	    }
#ifdef DEBUG
	  else
	    search_collisions++; // Wrong cell, move to next
#endif

      // Insert my cell into hash-table at the next free spot
#ifdef DEBUG
      while (hfill[cur>>5]&1U<<(cur&31))
	{
	  fill_collisions++;
	  cur = (cur+1)&hmask;
	}

#else
      while (hfill[cur>>5]&1U<<(cur&31))
	cur = (cur+1)&hmask;

#endif
      hfill[cur>>5] |= 1U<<(cur&31);
      
      htable[cur].idx = i;
      htable[cur].cell = cur_cell_id;
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
  if (sort_idx[0]<num_orig)
    domains[sort_idx[0]] = cells[domain].start;
  int64_t cur_cell_id = cell_ids[sort_idx[0]];
  for (size_t i=0,j=1;j<num_pos; ++j)
    {
      const int64_t p1=sort_idx[j];
      /* Find root of this domain (path compression of disjoint sets) */
      if (cell_ids[p1]!=cur_cell_id)
	{
	  cur_cell_id = cell_ids[p1];
	  domain = find_path_compress(++i, cells);
	}
      if (p1<num_orig) // Ignore images
	domains[p1] = cells[domain].start; // see hack above
    }
#ifdef DEBUG
  int max_rank = ranks[0];
  for (size_t i=1;i<num_cells;++i)
    if (ranks[i]>max_rank)
	max_rank = ranks[i];

  float frac_found = ngb_found/ (WALK_NGB_SIZE*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f average distance comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, search_collisions, num_doms, frac_found, ngb_found);
  printf("%.2f%% fill collisions (c.f. %.2f%% for perfect hashing)\n", 
	 fill_collisions*100.0/num_cells, 100.0*expected_collisions);
  printf("Max stack usage %d, max rank %d, point comparisons %d\n",
	 max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif

  free(ranks);
  free(cells);
  
  return num_doms;
}
