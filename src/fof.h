#ifndef FOF_H_
#define FOF_H_

//#define DEBUG  // Enable debug

/* Hash tables (primes per size */
#include <stdint.h>

// Hash table primes = nearest_prime(N/pi) (where N is the hash table size)
#define MAX_TAB_NO 23 // number of table sizes
static const int64_t TAB_START = 1024; // Corresponds to smallest hash table (& hash prime)
static const uint32_t HASH_PRIMES[MAX_TAB_NO] = 
  {331, 653, 1307, 2609, 5209, 10427, 20849, 41719, 83443, 166867, 333769, 
   667547, 1335043, 2670181, 5340359, 10680707, 21361421, 42722831, 85445653, 
   170891291, 341782633, 683565271, 1367130559};


#endif /* FOF_H_ included */
