/*
  Friends of Friends via hash grid
 */
//#define DEBUG  // Enable debug

#include <math.h>

#include <stdlib.h>
#include <stdint.h>
#ifdef DEBUG
  #include <stdio.h>
static int max_stack_usage=0;
#endif

// Magic number that specifies how full the hash-table should be to be most efficient
// Hash tables perform very badly beyond about 0.7
#define DESIRED_LOAD 0.6

// Hash table primes from planetmath.org
#define MAX_TAB_NO 22 // number of table sizes
static const unsigned int TAB_START = 1024; // Corresponds to smallest hash table (& hash prime)
static const int64_t HASH_PRIMES[MAX_TAB_NO] = {769,1543,3079,6151,12289, 24593,49157, 98317, 196613, 393241, 786433, 1572869, 3145739,6291469,12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};


typedef struct {
  unsigned int parent, rank;
} DisjointSet;

typedef struct {
  int64_t cell, // ID of cell
    idx; // index of cell (when ordered)
} HashCell;



void find_lattice(const double *pos, const int num_pos, const double inv_cell_width, const int nx, int64_t *out)
{
  /*
    Find the bucket index ((int)(p[0]*nx)*nx + (int)(p[1]*nx))*nx + (int)p[2]*nx 
    for every point
   */

  const int64_t NX = nx;
  for (size_t i=0;i<num_pos;i++)
    out[i] = ((int64_t)(pos[i*3]*inv_cell_width)*NX + (int64_t)(pos[i*3+1]*inv_cell_width))*NX + (int64_t)(pos[i*3+2]*inv_cell_width);

}
static inline unsigned int unite(const unsigned int my_root, const unsigned int adj_root, DisjointSet *restrict ds)
{
  // Disjoint sets union algorithm
  if (ds[my_root].rank < ds[adj_root].rank) // Add me to adj
    return (ds[my_root].parent = adj_root);
  else if (ds[my_root].rank > ds[adj_root].rank) // Add adj to me
    ds[adj_root].parent = my_root;
  else // Arbitrarily choose, in tests add adj to me slightly faster
    ds[ds[adj_root].parent = my_root].rank++;
  return my_root;
}
static inline unsigned int find_path_compress(const unsigned int idx, DisjointSet *restrict ds)
{
  /*
    Find root (domain) with path compression (i.e. set all parents to final root 
  */

  unsigned int stack_used=0, stack[16], domain = ds[idx].parent;

  for(unsigned int child=idx;child!=domain;domain = ds[child=domain].parent)
    stack[stack_used++] = child;
  
#ifdef DEBUG
  if (stack_used>max_stack_usage)
    printf("Increasing stack used to %d\n",(max_stack_usage = stack_used));
#endif

  // Set all those paths to have final parent
  while (stack_used)
    ds[stack[--stack_used]].parent = domain;

  return domain;
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

  int num_doms=0, num_cells=0;

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

  const int64_t hsize = TAB_START<<tab_size, hprime = HASH_PRIMES[tab_size],
    hmask = hsize-1;

  // Hashtable of indices
  HashCell *htable = malloc(hsize * sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  // Hash table fill bits
  uint32_t *hfill= calloc(hsize>>5, sizeof(uint32_t));
  if (!hfill)
    return -5;

  DisjointSet *ds = malloc(num_cells*sizeof(DisjointSet));
  if (!ds)
    return -4; // not enough mem.

#ifdef DEBUG
  printf("Platform int size %d\nHash cell size %d bytes, DisjointSet cell size %d bytes\n", 
	 (int)sizeof(int), (int)sizeof(HashCell), (int)sizeof(DisjointSet));
  int pt_cmp = 0, collisions=0, found=0;
#endif

  cell_start_end[0] = 0;
  // Loop over every (filled) cell
  for (size_t i=0, cur_end=0;i<num_cells; cell_start_end[++i] = cur_end)
    {
      // Put this in a new domain
      num_doms++;
      int64_t my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;

      const size_t cur_start = cur_end++,
	cur_cell_id = cells[sort_idx[cur_start]]; // ID for this cell

      /* Start and end of point data for this cell */
      if ((i+1)==num_cells)
	cur_end = num_pos;
      else
	while (cells[sort_idx[cur_end]]==cur_cell_id)
	  cur_end++;

      // Loop over adjacent cells (within sqrt(3)*cell_width)
      int64_t wanted_cell = cur_cell_id - walk_ngbs[0];
      for (int64_t adj=0, j=(wanted_cell*hprime)&hmask;1;)
	{
	  // Look up in hash table (search for cell or 0)
	  if (!(hfill[j>>5]>>(j&31)&1))
	    goto next_ngb; // No cell (absent)

	  if (htable[j].cell!=wanted_cell)
	    {
	      // Wrong cell, move to next
#ifdef DEBUG
	      collisions++;
#endif
	      j = (j+1)&hmask;
	      continue;
	    }
#ifdef DEBUG
	  found++;
#endif

	  // Found my cell, Find root of this domain (path compression of disjoint sets)
	  const int adj_idx = htable[j].idx;
	  const int64_t adj_root = find_path_compress(adj_idx, ds);
	  
	  /* Other domain to check for connection? */
	  if (adj_root!=my_root) 
	    {
	      /*
		Check pairwise for any connections, i.e. any point p1 in my
		cell within b of any point p2 in adj cell
	      */
	      for (size_t my_k=cur_start;my_k<cur_end;++my_k)
		{
		  const double *xyz1 = xyz + sort_idx[my_k]*3;

		  for (size_t adj_k=cell_start_end[adj_idx];
		       adj_k<cell_start_end[adj_idx+1]; adj_k++)
		    {
		      const double *xyz2 = xyz + sort_idx[adj_k]*3;
		      const double dx = xyz2[0]-xyz1[0],
			dy = xyz2[1]-xyz1[1],
			dz = xyz2[2]-xyz1[2];
#ifdef DEBUG		      
		      pt_cmp++; // Count comparisons
#endif
		      if ((dx*dx + dy*dy + dz*dz) > b2)
			continue;
		      
		      // Connect the domains
		      num_doms--;
		      my_root = unite(my_root, adj_root, ds);
		      
		      goto next_ngb;
		    }
		}
	    }
		
	next_ngb:
	  if ((++adj)==58)
	    break; // All done
	  
	  j=((wanted_cell = cur_cell_id - walk_ngbs[adj])*hprime)&hmask;
	}

      // Insert my cell into hash-table at the next free spot
      int64_t cur = (cur_cell_id*hprime)&hmask;
      while (hfill[cur>>5]&(1U<<(cur&31)))
	cur = (cur+1)&hmask;
      hfill[cur>>5] |= 1U<<(cur&31);
      
      htable[cur].idx = i;
      htable[cur].cell = cur_cell_id;
    }
#ifdef DEBUG
  float frac_found = found/ (58.0*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f Comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, collisions, num_doms, frac_found, found);
  int max_rank=0;
#endif

  for (size_t i=0;i<num_cells; ++i)
    {
#ifdef DEBUG
      if (ds[i].rank>max_rank)
	max_rank = ds[i].rank;
#endif
      /* Find root of this domain (path compression of disjoint sets) */
      const int64_t domain = find_path_compress(i, ds);

      // Set all particles in this cell to have this domain
      for (size_t j=cell_start_end[i];j<cell_start_end[i+1];++j)
	domains[sort_idx[j]] = domain;
    }
#ifdef DEBUG
  printf("Max stack usage %d, max rank %d, point comparisons %d\n",max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif
  free(hfill);
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

  int num_doms=0, num_cells=1;

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
    return 0; // No points no domains;

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

  const int64_t hsize = TAB_START<<tab_size, hprime = HASH_PRIMES[tab_size],
    hmask = hsize-1;

  // Hashtable of indices
  HashCell *htable = malloc(hsize * sizeof(HashCell));
  if (!htable)
    return -3; // not enough mem.

  // Hash table fill bits
  uint32_t *hfill= calloc(hsize>>5, sizeof(uint32_t));
  if (!hfill)
    return -5;

  DisjointSet *ds = malloc(num_cells*sizeof(DisjointSet));
  if (!ds)
    return -4; // not enough mem.

#ifdef DEBUG
  printf("Platform int size %d\nHash cell size %d bytes, DisjointSet cell size %d bytes\n", 
	 (int)sizeof(int), (int)sizeof(HashCell), (int)sizeof(DisjointSet));
  int pt_cmp = 0, collisions=0, found=0;
#endif

  cell_start_end[0] = 0; // Start point of first cell

  // Loop over every (filled) cell
  for (size_t i=0, cur_end=0;i<num_cells;cell_start_end[++i]=cur_end)
    {
      // Put this in a new domain
      num_doms++;
      int64_t my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;

      const size_t cur_start = cur_end++,
	cur_cell_id = cells[sort_idx[cur_start]]; // ID for this cell

      /* Start and end of point data for this cell */
      if ((i+1)==num_cells)
	cur_end = num_pos;
      else
	while (cells[sort_idx[cur_end]]==cur_cell_id)
	  cur_end++;


      // Check for any false images and link to cell of original
      
      for (size_t k=cur_start;k<cur_end;++k)
	{
	  const size_t p1 = sort_idx[k];
	  if (p1<num_orig) // Original (not an image)
	    continue; 
	  const int64_t orig_cell=cells[pad_idx[p1-num_orig]];

	  // Find cell in hash table (must be there)
	  int64_t cur = orig_cell*hprime;
	  while (htable[cur&=hmask].cell!=orig_cell)
	    cur++;

	  // Link me (if not already)
	  // Find root of this domain (path compression of disjoint sets)
	  const int64_t adj_root = find_path_compress(htable[cur].idx, ds);
	  
	  if (adj_root==my_root) // Already linked
	    continue;

	  // Connected (since my image in both cells)
	  num_doms--;
	  my_root = unite(my_root, adj_root, ds);

	}

      // Loop over adjacent cells (within sqrt(3)*cell_width)
      int adj=0;
      for (int64_t wanted_cell = cur_cell_id - walk_ngbs[0],j = (wanted_cell*hprime)&hmask; 1;) 
	{
	  // Look up in hash table (search for cell or 0)
	  if (!(hfill[j>>5]>>(j&31)&1))
	    goto next_ngb; // No cell (absent)
	  
	  if (htable[j].cell!=wanted_cell)
	    {
	      // Wrong cell, move to next
#ifdef DEBUG
	      collisions++;
#endif
	      j = (j+1)&hmask;
	      continue;
	    }
	  // Found my cell, Find root of this domain (path compression of disjoint sets)
#ifdef DEBUG
	  found++;
#endif

	  const int adj_idx = htable[j].idx;
	  const int64_t adj_root = find_path_compress(adj_idx, ds);
	  
	  /* Other domain to check for connection? */
	  if (adj_root!=my_root) 
	    {
	      /*
		Check pairwise for any connections, i.e. any point p1 in my
		cell within b of any point p2 in adj cell
	      */
	      for (size_t my_k=cur_start;my_k<cur_end;++my_k)
		{
		  const double *xyz1 = xyz + sort_idx[my_k]*3;
		  
		  for (size_t adj_k=cell_start_end[adj_idx];
		       adj_k<cell_start_end[adj_idx+1]; adj_k++)
		    {
		      const double *xyz2 = xyz + sort_idx[adj_k]*3;
		      const double dx = xyz2[0]-xyz1[0],
			dy = xyz2[1]-xyz1[1],
			dz = xyz2[2]-xyz1[2];
#ifdef DEBUG		      
		      pt_cmp++; // Count comparisons
#endif
		      if ((dx*dx + dy*dy + dz*dz) > b2)
			continue;
		      
		      // Connect the domains
		      num_doms--;
		      my_root = unite(my_root, adj_root, ds);
		      
		      goto next_ngb;
		    }
		}
	    }
	  
	next_ngb:
	  if ((++adj)==58)
	    break; // All done
	  
	  j=((wanted_cell = cur_cell_id - walk_ngbs[adj])*hprime)&hmask;
	}

      // Insert my cell into hash-table at the next free spot
      int64_t cur = (cur_cell_id*hprime)&hmask;
      while (hfill[cur>>5]&(1U<<(cur&31)))
	cur = (cur+1)&hmask;
      hfill[cur>>5] |= 1U<<(cur&31);
      
      htable[cur].idx = i;
      htable[cur].cell = cur_cell_id;

    }
#ifdef DEBUG
  float frac_found = found/ (58.0*num_cells),
    cmp_per_pt = (float)pt_cmp / (float)num_pos;
  printf("%.3f Comparisons per point (%d)\n%d hash collisions, %d domains created\n%.4f cells found/search (%d total)\n", cmp_per_pt, pt_cmp, collisions, num_doms, frac_found, found);
  int max_rank=0;
#endif

  for (size_t i=0;i<num_cells; ++i)
    {
#ifdef DEBUG
      if (ds[i].rank>max_rank)
	max_rank = ds[i].rank;
#endif

      /* Find root of this domain (path compression of disjoint sets) */
      const int64_t domain = find_path_compress(i, ds);

      // Set all particles in this cell to have this domain
      for (size_t j=cell_start_end[i];j<cell_start_end[i+1];++j)
	{
	  const size_t p1=sort_idx[j];
	  if (p1<num_orig) // Ignore images
	    domains[p1] = domain;
	}
    }
#ifdef DEBUG
  printf("Max stack usage %d, max rank %d, point comparisons %d\n",max_stack_usage, max_rank, pt_cmp);
  max_stack_usage=0;
#endif

  free(ds);
  free(hfill);
  free(htable);
  free(cell_start_end);

  return num_doms;
}
