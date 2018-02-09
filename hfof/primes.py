"""
For good hash tables we need large primes. Apply the Miller-Rabin test for 
small (<2^31) primes
"""
from __future__ import absolute_import, print_function

def _miller_rabin(a, n, u, t):
    """ 
    
    Miller-Rabin that a number is *definitely* composite (if not then it is
    *probably* prime, though you should try several to be sure

    a - test prime 
    n - number for testing
    u,t - decomposition of n,s.t. n=u 2^t +1 
    
    is a^t mod n != 1 *and*  (2^(ajt) mod n != -1 for all j in 0...r-1 
    """
    if pow(a, u, n) == 1:
        return False
    for j in range(t):
        if pow(a, u*(1<<j), n) == n-1:
            return False
    return True # n  is definitely composite
 
def is_prime(n):
    """ 
    Technically probababilistic prime testing method (Miller-Rabin), but 
    known to be exact for n<2^31
    """
    if n==1:
        return False
    elif n in (2,3,5,7):
        return True
    elif any((n % p)==0 for p in (2,3,5,7)):
        return False

    u, t = n - 1, 0
    while not u % 2:
        u, t = u >> 1, t + 1
    # Returns exact for n < 2^31, according 
    # C. Pomerance, J. L. Selfridge and Wagstaff, Jr., S. S., Math. Comp., 35:151 (1980) 1003-1026.
    return not any(_miller_rabin(a, n, u, t) for a in (2, 3, 5, 7))

def smallest_prime_atleast(m):
    """
    Find the first prime >= m
    """
    n=m
    while not is_prime(n):
        n+=1
    return n

if __name__=='__main__':
    assert(is_prime(17)==True)
    assert(is_prime(18)==False)
    print('Finding large hashtable prime for mod 2<32')
    n = int(1.5*2**31)
    best = smallest_prime_atleast(n)
    print('Searched', best-n, 'composites')
    print('Use', best)

    # Better primes
    tab0 = 1024
    from math import pi, e
    hash_primes = []
    for i in range(23):
        hsize = tab0<<i
        guess = hsize/pi
        nprime = smallest_prime_atleast(int(guess))
        guess = 2*guess -nprime # look lower as well as higher
        nprime2 = smallest_prime_atleast(int(guess))
        if guess-nprime2<nprime-guess:
            nprime = nprime2

        print(hsize, nprime)
        hash_primes.append(nprime)
    print('{'+', '.join(str(N) for N in hash_primes)+'};')
    print('{'+','.join(hex(N) for N in hash_primes)+'};')

