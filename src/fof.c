/*
  Friends of Friends via hash grid
 */
// #define DEBUG  // Enable debug

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

void find_lattice(const double *pos, const int num_pos, const double inv_cell_width, const int nx, int64_t *out)
{
  /*
    Find the bucket index ((int)(p[0]*nx)*nx + (int)(p[1]*nx))*nx + (int)p[2]*nx 
    for every point
   */

  const int64_t NX = nx;
  for (int i=0;i<num_pos;i++)
    out[i] = ((int64_t)(pos[i*3]*inv_cell_width)*NX + (int64_t)(pos[i*3+1]*inv_cell_width))*NX + (int64_t)(pos[i*3+2]*inv_cell_width);

}


int fof_link_cells(const int num_pos, const int N,const double rcut, const int64_t *restrict cells, const int64_t *restrict sort_idx, int32_t *restrict domains, const double *restrict xyzw)
{
  // This version uses disjoint sets algo, sep list
  /*
     cells - (num_cells x 3) array of cell lattice coord, start (in pos) and end (in pos)
  */

  int stack_used=0, num_doms=0;
  int stack[16];
  const int N2 = N*N;
  const int64_t cell_shift[58] = {N2, +1, +N, N2+N, N2-N, N2-1, N2+1, +N+1, +N-1, N2-N+1, N2+N+1, N2-N-1, N2+N-1, 2*N2, +2*N, +2, 2*N2+N, N2+2*N, 2*N2-N, N2-2*N, +N+2, +N-2, N2+2, 2*N2+1, N2-2, +2*N-1, +2*N+1, 2*N2-1, 2*N2-N-1, 2*N2+N-1, N2+2*N-1, N2-2*N+1, N2+N-2, 2*N2-N+1, 2*N2+N+1, N2+2*N+1, N2-N+2, N2+N+2, N2-N-2, N2-2*N-1, 2*N2+2*N, +2*N-2, 2*N2-2*N, 2*N2-2, 2*N2+2, +2*N+2, 2*N2+N-2, 2*N2-N-2, N2+2*N+2, 2*N2-2*N-1, 2*N2+2*N-1, 2*N2-2*N+1, 2*N2+2*N+1, N2-2*N+2, 2*N2-N+2, 2*N2+N+2, N2+2*N-2, N2-2*N-2};

  const double rcut2 = rcut*rcut;


  /* Count the number of cells */
  int num_cells = 1;
  int64_t last_id = cells[sort_idx[0]];
  for (int i=1;i<num_pos;++i)
    {
      if (cells[sort_idx[i]]!=last_id)
	{
	  num_cells++;
	  last_id = cells[sort_idx[i]];
	}
    }
#ifdef DEBUG
  printf("Number of cells %d\n", num_cells);
#endif

  uint32_t* cell_start_end = malloc((num_cells+1)*sizeof(uint32_t));
  if (!cell_start_end)
    return -1;

  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) tab_size++;
  if (tab_size==MAX_TAB_NO)
    return -1; // Table too big

  // Unsigned ints have guaranteed 'wrapping'
  const unsigned int hsize = TAB_START<<tab_size, hprime = HASH_PRIMES[tab_size];
  const unsigned int hmask = hsize-1;

  // Zero-ed hashtable (index + start + end for each cell)
  struct HashCell{
    int64_t icell;
    unsigned int idx;
  };

  struct HashCell *htable = calloc(hsize, sizeof(struct HashCell));
  if (!htable)
    return -2; // not enough mem.

  struct DisjointSet{
    unsigned int parent, rank;
  };

  struct DisjointSet *ds = malloc(num_cells*sizeof(struct DisjointSet));
  if (!ds)
    return -3; // not enough mem.

#ifdef DEBUG
  printf("integer size %d hash cell size %d DisjointSet cell size %d\n", (int)sizeof(int), (int)sizeof(struct HashCell), (int)sizeof(struct DisjointSet));
  int pt_cmp = 0, collisions=0, insertions=0, max_stack_usage=0;
#endif



  cell_start_end[0] = 0; // Start point of first cell

  // Loop over every (filled) cell
  for (unsigned int i=0;i<num_cells;++i)
    {
#ifdef DEBUG
      if (i%100000==0)
	{
	  printf("Working on cell %d, pt comparisons %d, hash collisions %d list insertions %d domains created %d\n",i, pt_cmp, collisions, insertions, num_doms);
	  pt_cmp = 0;
	  collisions=0;
	  insertions=0;
	  //	  num_doms = 0;
	}
#endif

      // Put this in a new domain
      num_doms++;

      const uint32_t cur_start = cell_start_end[i];
      const int64_t cur_cell_id = cells[sort_idx[cur_start]]; // ID for this cell

      /* Start and end of point data for this cell */
      uint32_t cur_end = cur_start+1;
      if ((i+1)<num_cells)
	{
	  while (cells[sort_idx[cur_end]]==cur_cell_id)
	    cur_end++;
	}
      else
	cur_end = num_pos;
      cell_start_end[i+1] = cur_end;

      // Insert my cell into hash-table at the next free spot
      unsigned int cur; // Cursor for linked list
      for (cur=((unsigned int)cur_cell_id*hprime)&hmask; htable[cur].icell; cur=(cur+1)&hmask);

      htable[cur].icell = ~cur_cell_id; // Bit twiddle to make sure never 0
      htable[cur].idx = i;
      unsigned int my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;
      
      // Loop over adjacent cells (within sqrt(3)*cell_width)
      for (int adj=0;adj<58;++adj)
	{
	  const int64_t wanted_cell = cur_cell_id - cell_shift[adj];

	  // Look up in hash table (search for cell or 0)

	  for (unsigned int j=((unsigned int)wanted_cell*hprime)&hmask;htable[j].icell;j=(j+1)&hmask)
	    {
	      if ((~htable[j].icell)!=wanted_cell)
		{
#ifdef DEBUG
		  collisions++;
#endif
		  continue;
		}

	      // Found
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

	      if (adj_root!=my_root) // Diff domain, might be connected
		{
		  // Check pairwise for any connections
		  for (int k1=cur_start;k1<cur_end;++k1)
		    {
		      const int v1 = sort_idx[k1];
		      const double x1 = xyzw[v1*3], y1 = xyzw[v1*3+1], z1 = xyzw[v1*3+2];
		      for (int k2=cell_start_end[adj_idx],v2=sort_idx[k2];cells[v2]==wanted_cell; v2=sort_idx[++k2])
			{
			  const double dx = xyzw[v2*3]-x1,
			    dy = xyzw[v2*3+1]-y1,
			    dz = xyzw[v2*3+2]-z1;
#ifdef DEBUG		      
			  pt_cmp++;
#endif
			  if ((dx*dx + dy*dy + dz*dz) > rcut2)
			    continue;

			  num_doms--;
			  
			  // Disjoint sets union
			  if (ds[my_root].rank < ds[adj_root].rank)
			    my_root = ds[my_root].parent = adj_root;
			  else if (ds[my_root].rank > ds[adj_root].rank)
			    ds[adj_root].parent = my_root;
			  else
			    {
			      // Arbitrarily choose (in tests this slightly faster)
			      ds[my_root].rank++;
			      ds[adj_root].parent = my_root;
			    }

			  goto connected;
			}
		    }
		}

	      break;
	    }
	    connected:
	  continue;
	}
 
    }
#ifdef DEBUG
  printf("Number of domains %d\n", num_doms);
#endif

  for (int i=0;i<num_cells; ++i)
    {
      // Find root of this domain (path compression of disjoint sets)

      unsigned int domain = ds[i].parent;

      for(unsigned int child=i;child!=domain;domain = ds[child=domain].parent)
	stack[stack_used++] = child;

#ifdef DEBUG
      max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
      // Set all those paths to have final parent
      while (stack_used)
	ds[stack[--stack_used]].parent = domain;
      /* Done with (path-compressing) find */

      // Set all particles to have this domain
      for (int j=cell_start_end[i];j<cell_start_end[i+1];++j)
	domains[sort_idx[j]] = domain;
    }
#ifdef DEBUG
  printf("Max stack usage %d, point comparisons %d\n",max_stack_usage, pt_cmp);
#endif
  free(cell_start_end);
  free(htable);
  free(ds);
  return num_doms;
}

int fof_link_cells_nosort(const int num_cells, const int N,const double rcut, const int64_t *restrict cell_ids, const int32_t *restrict cell_start_end, int32_t *restrict domains, const double *restrict xyzw)
{
  // This version uses disjoint sets algo, sep list
  /*
     cells - (num_cells x 3) array of cell lattice coord, start (in pos) and end (in pos)
  */

  int stack_used=0, num_doms=0;
  int stack[16];
  const int N2 = N*N;
  const int64_t cell_shift[58] = {N2, +1, +N, N2+N, N2-N, N2-1, N2+1, +N+1, +N-1, N2-N+1, N2+N+1, N2-N-1, N2+N-1, 2*N2, +2*N, +2, 2*N2+N, N2+2*N, 2*N2-N, N2-2*N, +N+2, +N-2, N2+2, 2*N2+1, N2-2, +2*N-1, +2*N+1, 2*N2-1, 2*N2-N-1, 2*N2+N-1, N2+2*N-1, N2-2*N+1, N2+N-2, 2*N2-N+1, 2*N2+N+1, N2+2*N+1, N2-N+2, N2+N+2, N2-N-2, N2-2*N-1, 2*N2+2*N, +2*N-2, 2*N2-2*N, 2*N2-2, 2*N2+2, +2*N+2, 2*N2+N-2, 2*N2-N-2, N2+2*N+2, 2*N2-2*N-1, 2*N2+2*N-1, 2*N2-2*N+1, 2*N2+2*N+1, N2-2*N+2, 2*N2-N+2, 2*N2+N+2, N2+2*N-2, N2-2*N-2};

  const double rcut2 = rcut*rcut;
  // Find the next power of 2 large enough to hold table at desired load
  int tab_size=0;
  while (tab_size<MAX_TAB_NO && (TAB_START<<tab_size)*DESIRED_LOAD<num_cells) tab_size++;
  if (tab_size==MAX_TAB_NO)
    return -1; // Table too big

  // Unsigned ints have guaranteed 'wrapping'
  const unsigned int hsize = (TAB_START<<tab_size), hprime = HASH_PRIMES[tab_size];
  const unsigned int hmask = hsize-1;

  // Zero-ed hashtable (index + start + end for each cell)
  struct HashCell{
    int64_t icell;
    unsigned int idx;
  };

  struct HashCell *htable = calloc(hsize, sizeof(struct HashCell));
  if (!htable)
    return -2; // not enough mem.

  struct DisjointSet{
    unsigned int parent, rank;
  };

  struct DisjointSet *ds = malloc(num_cells*sizeof(struct DisjointSet));
  if (!ds)
    return -3; // not enough mem.

#ifdef DEBUG
  printf("integer size %d hash cell size %d DisjointSet cell size %d\n", (int)sizeof(int), (int)sizeof(struct HashCell), (int)sizeof(struct DisjointSet));
  int pt_cmp = 0, collisions=0, insertions=0, max_stack_usage=0;
#endif

  // Loop over every (filled) cell
  for (unsigned int i=0;i<num_cells;++i)
    {
#ifdef DEBUG
      if (i%100000==0)
	{
	  printf("Working on cell %d, pt comparisons %d, hash collisions %d list insertions %d domains created %d\n",i, pt_cmp, collisions, insertions, num_doms);
	  pt_cmp = 0;
	  collisions=0;
	  insertions=0;
	  num_doms = 0;
	}
#endif

      num_doms++;

      // Put this in a new domain
      unsigned int cur; // Cursor for linked list      
      // Insert my cell into hash-table at the next free spot
      for (cur=((unsigned int)cell_ids[i]*hprime)&hmask; htable[cur].icell; cur=(cur+1)&hmask);

      htable[cur].icell = ~cell_ids[i]; // Bit twiddle to make sure never 0
      htable[cur].idx = i;
      unsigned int my_root = ds[i].parent = i; // Own domain
      ds[i].rank = 0;
      
      // Loop over adjacent cells (within sqrt(3)*cell_width)
      for (int adj=0;adj<58;++adj)
	{
	  const int64_t wanted_cell = cell_ids[i] - cell_shift[adj];

	  // Look up in hash table (search for cell or 0)

	  for (unsigned int j=((unsigned int)wanted_cell*hprime)&hmask;htable[j].icell;j=(j+1)&hmask)
	    {
	      if ((~htable[j].icell)!=wanted_cell)
		{
#ifdef DEBUG
		  collisions++;
#endif
		  continue;
		}

	      // Found
	      const int adj_idx = htable[j].idx;
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


	      if (adj_root!=my_root) // Diff domain, might be connected
		{
		  // Check pairwise for any connections
		  const int a1 = cell_start_end[adj_idx*2], a2 = cell_start_end[adj_idx*2+1];
		  for (int k1=cell_start_end[i*2];k1<cell_start_end[i*2+1];++k1)
		    {
		      const double x1 = xyzw[k1*3], y1 = xyzw[k1*3+1], z1 = xyzw[k1*3+2];
		      for (int k2=a1;k2<a2;++k2)
			{
			  const double dx = xyzw[k2*3]-x1,
			    dy = xyzw[k2*3+1]-y1,
			    dz = xyzw[k2*3+2]-z1;
#ifdef DEBUG		      
			  pt_cmp++;
#endif
			  if ((dx*dx + dy*dy + dz*dz) > rcut2)
			    continue;

			  num_doms--;
			  
			  // Disjoint sets union
			  if (ds[my_root].rank < ds[adj_root].rank)
			    my_root = ds[my_root].parent = adj_root;
			  else if (ds[my_root].rank > ds[adj_root].rank)
			    ds[adj_root].parent = my_root;
			  else
			    {
			      // Arbitrarily choose (in tests this slightly faster)
			      ds[my_root].rank++;
			      ds[adj_root].parent = my_root;
			    }

			  goto connected;
			}
		    }
		}

	      break;
	    }
	    connected:
	  continue;
	}
 
    }
  // TODO put domain number into rank of root.

  for (int i=0;i<num_cells; ++i)
    {
      // Find root of this domain (path compression of disjoint sets)

      unsigned int domain = ds[i].parent;

      for(unsigned int child=i;child!=domain;domain = ds[child=domain].parent)
	stack[stack_used++] = child;

#ifdef DEBUG
      max_stack_usage = stack_used>max_stack_usage ? stack_used : max_stack_usage;
#endif
      // Set all those paths to have final parent
      while (stack_used)
	ds[stack[--stack_used]].parent = domain;
      /* Done with (path-compressing) find */
      
      domains[i] = domain;
    }
#ifdef DEBUG
  printf("Max stack usage %d, point comparisons %d\n",max_stack_usage, pt_cmp);
#endif

  free(htable);
  free(ds);
  return num_doms;
}


