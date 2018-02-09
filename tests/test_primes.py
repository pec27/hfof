from hfof.primes import is_prime, smallest_prime_atleast
# Primes up to 200
small_primes = (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,
                83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,
                167,173,179,181,191,193,197,199)
def test_200():
    """
    test known primes
    """
    for i in range(2,200):
        if is_prime(i) != (i in small_primes):
            raise Exception('Failed at %d'%i)

def test_smallest_prime_atleast():
    """ Finding next prime """
    for i in range(len(small_primes)-1):
        assert(smallest_prime_atleast(small_primes[i])==small_primes[i])
        for j in range(small_primes[i]+1, small_primes[i+1]):
            assert(smallest_prime_atleast(j)==small_primes[i+1])            

if __name__=='__main__':    
    test_200()
