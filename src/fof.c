/*
  Friends of Friends via hash grid
 */
//#define DEBUG  // Enable debug

#include <math.h>

#include <stdlib.h>
#include <stdint.h>
#ifdef DEBUG
  #include <stdio.h>
#endif

// Magic number that specifies how full the hash-table should be to be most efficient
// Hash tables perform very badly beyond about 0.7
#define DESIRED_LOAD 0.6

// Hash table primes from planetmath.org
#define MAX_TAB_NO 22 // number of table sizes
static const unsigned int TAB_START = 1024; // Corresponds to smallest hash table (& hash prime)
static const unsigned int HASH_PRIMES[MAX_TAB_NO] = {769,1543,3079,6151,12289, 24593,49157, 98317, 196613, 393241, 786433, 1572869, 3145739,6291469,12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};


typedef struct {
    unsigned int parent, rank;
} DisjointSet;

typedef struct {
  int64_t icell; // Bit twiddle cell to make sure nonzero
  unsigned int idx; // index of cell (when ordered)
} HashCell;



void find_lattice(const double *pos, const int num_pos, const double inv_cell_width, const int nx, int64_t *out)
{
  /*
    Find the bucket index ((int)(p[0]*nx)*nx + (int)(p[1]*nx))*nx + (int)p[2]*nx 
    for every point
   */

  const int64_t NX = nx;
  for (unsigned int i=0;i<num_pos;i++)
    out[i] = ((int64_t)(pos[i*3]*inv_cell_width)*NX + (int64_t)(pos[i*3+1]*inv_cell_width))*NX + (int64_t)(pos[i*3+2]*inv_cell_width);

}

int fof_link_cells(const int num_pos, const int N, const double b, 
		   const double *restrict xyz, const int64_t *restrict cells, 
		   const int64_t *restrict sort_idx, int32_t *restrict domains)
{
  /*
    Friends of Friends by binning into cells.

    num_pos  - The number of points
    N        - Multiplicative factor such that cell(i,j,k) = N^2 i + N j + k. 
               Usually 3 mod 4 for maximum bitwise independence
    b        - The linking length
    xyz      - (num_pos * 3) array of (unsorted) points
    cells    - The cell values for each point (s.t. i=x/cell_width etc.)
    sort_idx - An index such that cells[sort_idx] is ascending (ordered points
               by cell)
    domains  - (num_pos,) integer output array for the domains

    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */

  int stack_used=0, num_doms=0, stack[16], num_cells=0;

  const int M = N*N;
  const int walk_ngbs[58] = {M, 1, N, M+N, M-N, M-1, M+1, N+1, N-1, M-N+1, 
			     M+N+1, M-N-1, M+N-1, 2*M, 2*N, 2, 2*M+N, M+2*N, 
			     2*M-N, M-2*N, N+2, N-2, M+2, 2*M+1, M-2, 2*N-1, 
			     2*N+1, 2*M-1, 2*M-N-1, 2*M+N-1, M+2*N-1, M-2*N+1, 
			     M+N-2, 2*M-N+1, 2*M+N+1, M+2*N+1, M-N+2, M+N+2, 
			     M-N-2, M-2*N-1, 2*M+2*N, 2*N-2, 2*M-2*N, 2*M-2, 
			     2*M+2, 2*N+2, 2*M+N-2, 2*M-N-2, M+2*N+2, 2*M-2*N-1, 
			     2*M+2*N-1, 2*M-2*N+1, 2*M+2*N+1, M-2*N+2, 2*M-N+2, 
			     2*M+N+2, M+2*N-2, M-2*N-2};

  const double b2 = b*b;


  if (!num_pos)
    return 0; // No points => no domains!

  /* Count the number of cells */
  for (int64_t i=1,last_id = cells[sort_idx[num_cells++]];i<num_pos;++i)
    if (cells[sort_idx[i]]!=last_id)
      {
	num_cells++;
	last_id = cells[sort_idx[i]];
      }

#ifdef DEBUG
  printf("Number of cells %d\n", num_cells);
#endif

  uint32_t* cell_start_end = malloc((num_cells+1)*sizeof(uint32_t));
  if (!cell_start_end)
    return -1;

  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  // Unsigned ints have guaranteed 'wrapping'
  const unsigned int hsize = TAB_START<<tab_size, hprime = HASH_PRIMES[tab_size];
  const unsigned int hmask = hsize-1;

  // Hashtable of indices
  HashCell *htable = calloc(hsize, sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  DisjointSet *ds = malloc(num_cells*sizeof(DisjointSet));
  if (!ds)
    return -4; // not enough mem.

#ifdef DEBUG
  printf("Platform int size %d\nHash cell size %d bytes, DisjointSet cell size %d bytes\n", 
	 (int)sizeof(int), (int)sizeof(HashCell), (int)sizeof(DisjointSet));
  int pt_cmp = 0, collisions=0, max_stack_usage=0;
#endif


  cell_start_end[0] = 0; // Start point of first cell

  // Loop over every (filled) cell
  for (unsigned int i=0;i<num_cells;++i)
    {
      // Put this in a new domain
      num_doms++;

      const uint32_t cur_start = cell_start_end[i];
      const int64_t cur_cell_id = cells[sort_idx[cur_start]]; // ID for this cell

      /* Start and end of point data for this cell */
      uint32_t cur_end = cur_start+1;
      if ((i+1)==num_cells)
	cur_end = num_pos;
      else
	while (cells[sort_idx[cur_end]]==cur_cell_id)
	  cur_end++;

      cell_start_end[i+1] = cur_end;

      // Insert my cell into hash-table at the next free spot
      unsigned int cur=((unsigned int)cur_cell_id*hprime)&hmask; 
      while (htable[cur].icell) 
	cur=(cur+1)&hmask;

      htable[cur].icell = ~cur_cell_id; // Bit twiddle to make sure never 0
      htable[cur].idx = i;
      unsigned int my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;
      
      // Loop over adjacent cells (within sqrt(3)*cell_width)
      int64_t wanted_cell = cur_cell_id - walk_ngbs[0];
      for (unsigned int adj=0, j=((unsigned int)wanted_cell*hprime)&hmask;1;)
	{
	  // Look up in hash table (search for cell or 0)
	  if (!htable[j].icell)
	    goto connected; // No cell (absent)

	  if ((~htable[j].icell)!=wanted_cell)
	    {
	      // Wrong cell, move to next
#ifdef DEBUG
	      collisions++;
#endif
	      j = (j+1)&hmask;
	      continue;
	    }
	  // Found my cell
	  const unsigned int adj_idx = htable[j].idx;
	  /*************************************/
	  // Find root of this domain (path compression of disjoint sets)

	  unsigned int adj_root = ds[adj_idx].parent;
	  
	  for(unsigned int child=adj_idx;child!=adj_root;adj_root = ds[child=adj_root].parent)
	    stack[stack_used++] = child;
	  
#ifdef DEBUG
	  if (stack_used>max_stack_usage)
	    printf("Increasing stack used to %d\n",(max_stack_usage = stack_used));
#endif
	  // Set all those paths to have final parent
	  while (stack_used)
	    ds[stack[--stack_used]].parent = adj_root;

	  /* Done with (path-compressing) find */
	  /*************************************/
	  
	  /* Other domain to check for connection? */
	  if (adj_root!=my_root) 
	    {
	      /*
		Check pairwise for any connections, i.e. any point p1 in my
		cell within b of any point p2 in adj cell
	      */
	      for (unsigned int my_k=cur_start;my_k<cur_end;++my_k)
		{
		  const unsigned int p1 = sort_idx[my_k];
		  const double x1 = xyz[p1*3], y1 = xyz[p1*3+1], z1 = xyz[p1*3+2];
		  for (unsigned int adj_k=cell_start_end[adj_idx],p2=sort_idx[adj_k];
		       cells[p2]==wanted_cell; 
		       p2=sort_idx[++adj_k])
		    {
		      const double dx = xyz[p2*3]-x1,
			dy = xyz[p2*3+1]-y1,
			dz = xyz[p2*3+2]-z1;
#ifdef DEBUG		      
		      pt_cmp++; // Count comparisons
#endif
		      if ((dx*dx + dy*dy + dz*dz) > b2)
			continue;
		      
		      /* Disjoint sets union algorithm */
		      num_doms--;
		      if (ds[my_root].rank < ds[adj_root].rank) // Add me to adj
			my_root = ds[my_root].parent = adj_root;
		      else if (ds[my_root].rank > ds[adj_root].rank) // Add adj to me
			ds[adj_root].parent = my_root;
		      else // Arbitrarily choose, in tests add adj to me slightly faster
			ds[ds[adj_root].parent = my_root].rank++;
		      /* Union done */
		      
		      goto connected;
		    }
		}
	    }
		
	connected:
	  if ((++adj)==58)
	    break; // All done
	  
	  j=((unsigned int)(wanted_cell = cur_cell_id - walk_ngbs[adj])*hprime)&hmask;
	}
    }
#ifdef DEBUG
  printf("Connections complete, pt comparisons %d, hash collisions %d, domains created %d\n", pt_cmp, collisions, num_doms);
#endif

  for (int i=0;i<num_cells; ++i)
    {

      /* Find root of this domain (path compression of disjoint sets) */

      unsigned int domain = ds[i].parent;
      for(unsigned int child=i;child!=domain;domain = ds[child=domain].parent)
	stack[stack_used++] = child;

#ifdef DEBUG
      max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
      
      while (stack_used) // Set all those paths to have final parent
	ds[stack[--stack_used]].parent = domain;

      /* Done with (path-compressing) find */

      // Set all particles in this cell to have this domain
      for (int j=cell_start_end[i];j<cell_start_end[i+1];++j)
	domains[sort_idx[j]] = domain;
    }
#ifdef DEBUG
  printf("Max stack usage %d, point comparisons %d\n",max_stack_usage, pt_cmp);
#endif

  free(ds);
  free(htable);
  free(cell_start_end);

  return num_doms;
}

int fof_periodic(const int num_pos, const int N, const int num_orig, const double b, 
		 const double *restrict xyz, const int64_t *restrict cells, 
		 const int64_t *restrict sort_idx, const int64_t* restrict pad_idx,
		 int32_t *restrict domains)
{
  /*
    Friends of Friends by binning into cells.

    num_pos  - The number of points
    N        - Multiplicative factor such that cell(i,j,k) = N^2 i + N j + k. 
               Usually 3 mod 4 for maximum bitwise independence
    num_orig - Number of points that are not 'false' images
    b        - The linking length
    xyz      - (num_pos * 3) array of (unsorted) points
    cells    - The cell values for each point (s.t. i=x/cell_width etc.)
    sort_idx - An index such that cells[sort_idx] is ascending (ordered points
               by cell)
    pad_idx  - the indices of the images in the original array (for linking)
    domains  - (num_pos,) integer output array for the domains

    Returns filled domains, integers s.t.  (same integer)=>connected

    Implementation uses a hash-table and the disjoint sets algorithm to link
    cells.
  */

  int stack_used=0, num_doms=0, stack[16];

  const int M = N*N;
  const int walk_ngbs[58] = {M, 1, N, M+N, M-N, M-1, M+1, N+1, N-1, M-N+1, 
			     M+N+1, M-N-1, M+N-1, 2*M, 2*N, 2, 2*M+N, M+2*N, 
			     2*M-N, M-2*N, N+2, N-2, M+2, 2*M+1, M-2, 2*N-1, 
			     2*N+1, 2*M-1, 2*M-N-1, 2*M+N-1, M+2*N-1, M-2*N+1, 
			     M+N-2, 2*M-N+1, 2*M+N+1, M+2*N+1, M-N+2, M+N+2, 
			     M-N-2, M-2*N-1, 2*M+2*N, 2*N-2, 2*M-2*N, 2*M-2, 
			     2*M+2, 2*N+2, 2*M+N-2, 2*M-N-2, M+2*N+2, 2*M-2*N-1, 
			     2*M+2*N-1, 2*M-2*N+1, 2*M+2*N+1, M-2*N+2, 2*M-N+2, 
			     2*M+N+2, M+2*N-2, M-2*N-2};

  const double b2 = b*b;

  /* Count the number of cells */
  int num_cells = 1;
  int64_t last_id = cells[sort_idx[0]];
  for (int i=1;i<num_pos;++i)
    if (cells[sort_idx[i]]!=last_id)
      {
	num_cells++;
	last_id = cells[sort_idx[i]];
      }

#ifdef DEBUG
  printf("Number of cells %d\n", num_cells);
#endif

  uint32_t* cell_start_end = malloc((num_cells+1)*sizeof(uint32_t));
  if (!cell_start_end)
    return -1;

  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) 
    tab_size++;

  if (tab_size==MAX_TAB_NO)
    return -2; // Table too big

  // Unsigned ints have guaranteed 'wrapping'
  const unsigned int hsize = TAB_START<<tab_size, hprime = HASH_PRIMES[tab_size];
  const unsigned int hmask = hsize-1;

  // Hashtable of indices
  HashCell *htable = calloc(hsize, sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  DisjointSet *ds = malloc(num_cells*sizeof(DisjointSet));
  if (!ds)
    return -4; // not enough mem.

#ifdef DEBUG
  printf("Platform int size %d\nHash cell size %d bytes, DisjointSet cell size %d bytes\n", 
	 (int)sizeof(int), (int)sizeof(HashCell), (int)sizeof(DisjointSet));
  int pt_cmp = 0, collisions=0, max_stack_usage=0;
#endif

  cell_start_end[0] = 0; // Start point of first cell

  // Loop over every (filled) cell
  for (unsigned int i=0;i<num_cells;++i)
    {
      // Put this in a new domain
      num_doms++;

      const uint32_t cur_start = cell_start_end[i];
      const int64_t cur_cell_id = cells[sort_idx[cur_start]]; // ID for this cell

      /* Start and end of point data for this cell */
      uint32_t cur_end = cur_start+1;
      if ((i+1)==num_cells)
	cur_end = num_pos;
      else
	while (cells[sort_idx[cur_end]]==cur_cell_id)
	  cur_end++;

      cell_start_end[i+1] = cur_end;

      // Insert my cell into hash-table at the next free spot
      unsigned int cur=((unsigned int)cur_cell_id*hprime)&hmask; 
      while (htable[cur].icell) 
	cur=(cur+1)&hmask;

      htable[cur].icell = ~cur_cell_id; // Bit twiddle to make sure never 0
      htable[cur].idx = i;
      unsigned int my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;

      // Check for any false images and link to cell of original
      
      for (int k=cur_start;k<cur_end;++k)
	{
	  const int64_t p1 = sort_idx[k];
	  if (p1<num_orig) // Original (not an image)
	    continue; 
	  const int64_t orig_cell=cells[pad_idx[p1-num_orig]];

	  // Find cell in hash table (must be there)
	  unsigned int j;
	  for (j=((unsigned int)orig_cell*hprime)&hmask;
	       (~htable[j].icell)!=orig_cell; j&=hmask)
	    j++;

	  const int adj_idx = htable[j].idx;

	  // Link me (if not already)

	  // *************************************
	  // Find root of this domain (path compression of disjoint sets)

	  unsigned int adj_root = ds[adj_idx].parent;

	  for(unsigned int child=adj_idx;child!=adj_root;adj_root = ds[child=adj_root].parent)
	    stack[stack_used++] = child;

#ifdef DEBUG
	  max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
	  // Set all those paths to have final parent
	  while (stack_used)
	    ds[stack[--stack_used]].parent = adj_root;

	  // Done with (path-compressing) find
	  // *************************************
	  
	  if (adj_root==my_root) // Already linked
	    continue;

	  // Connected (since my image in both cells)
	  // Disjoint sets union algorithm
	  num_doms--;
	  if (ds[my_root].rank < ds[adj_root].rank) // Add me to adj
	    my_root = ds[my_root].parent = adj_root;
	  else if (ds[my_root].rank > ds[adj_root].rank) // Add adj to me
	    ds[adj_root].parent = my_root;
	  else // Arbitrarily choose, in tests add adj to me slightly faster
	    ds[ds[adj_root].parent = my_root].rank++;
	  // Union done 
	}
      
      // Loop over adjacent cells (within sqrt(3)*cell_width)
      int64_t wanted_cell = cur_cell_id - walk_ngbs[0];
      for (unsigned int adj=0, j=((unsigned int)wanted_cell*hprime)&hmask;1;)
	{
	  // Look up in hash table (search for cell or 0)
	  if (!htable[j].icell)
	    goto connected; // No cell (absent)

	  if ((~htable[j].icell)!=wanted_cell)
	    {
	      // Wrong cell, move to next
#ifdef DEBUG
	      collisions++;
#endif
	      j = (j+1)&hmask;
	      continue;
	    }
	  // Found my cell
	  const int adj_idx = htable[j].idx;
	  /*************************************/
	  // Find root of this domain (path compression of disjoint sets)

	  unsigned int adj_root = ds[adj_idx].parent;
	  
	  for(unsigned int child=adj_idx;child!=adj_root;adj_root = ds[child=adj_root].parent)
	    stack[stack_used++] = child;
	  
#ifdef DEBUG
	  max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
	  // Set all those paths to have final parent
	  while (stack_used)
	    ds[stack[--stack_used]].parent = adj_root;

	  /* Done with (path-compressing) find */
	  /*************************************/
	  
	  /* Other domain to check for connection? */
	  if (adj_root!=my_root) 
	    {
	      /*
		Check pairwise for any connections, i.e. any point p1 in my
		cell within b of any point p2 in adj cell
	      */
	      for (int my_k=cur_start;my_k<cur_end;++my_k)
		{
		  const unsigned int p1 = sort_idx[my_k];
		  const double x1 = xyz[p1*3], y1 = xyz[p1*3+1], z1 = xyz[p1*3+2];
		  for (unsigned int adj_k=cell_start_end[adj_idx],p2=sort_idx[adj_k];
		       cells[p2]==wanted_cell; 
		       p2=sort_idx[++adj_k])
		    {
		      const double dx = xyz[p2*3]-x1,
			dy = xyz[p2*3+1]-y1,
			dz = xyz[p2*3+2]-z1;
#ifdef DEBUG		      
		      pt_cmp++; // Count comparisons
#endif
		      if ((dx*dx + dy*dy + dz*dz) > b2)
			continue;
		      
		      /* Disjoint sets union algorithm */
		      num_doms--;
		      if (ds[my_root].rank < ds[adj_root].rank) // Add me to adj
			my_root = ds[my_root].parent = adj_root;
		      else if (ds[my_root].rank > ds[adj_root].rank) // Add adj to me
			ds[adj_root].parent = my_root;
		      else // Arbitrarily choose, in tests add adj to me slightly faster
			ds[ds[adj_root].parent = my_root].rank++;
		      /* Union done */
		      
		      goto connected;
		    }
		}
	    }
		
	connected:
	  if ((++adj)==58)
	    break; // All done
	  
	  j=((unsigned int)(wanted_cell = cur_cell_id - walk_ngbs[adj])*hprime)&hmask;
	}
    }
#ifdef DEBUG
  printf("Connections complete, pt comparisons %d, hash collisions %d, domains created %d\n", pt_cmp, collisions, num_doms);
#endif

  for (int i=0;i<num_cells; ++i)
    {

      /* Find root of this domain (path compression of disjoint sets) */

      unsigned int domain = ds[i].parent;
      for(unsigned int child=i;child!=domain;domain = ds[child=domain].parent)
	stack[stack_used++] = child;

#ifdef DEBUG
      max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
      
      while (stack_used) // Set all those paths to have final parent
	ds[stack[--stack_used]].parent = domain;

      /* Done with (path-compressing) find */

      // Set all particles in this cell to have this domain
      for (int j=cell_start_end[i];j<cell_start_end[i+1];++j)
	{
	  const int p1=sort_idx[j];
	  if (p1<num_orig)
	    domains[p1] = domain;
	}
    }
#ifdef DEBUG
  printf("Max stack usage %d, point comparisons %d\n",max_stack_usage, pt_cmp);
#endif

  free(ds);
  free(htable);
  free(cell_start_end);

  return num_doms;
}

