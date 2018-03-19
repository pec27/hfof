from __future__ import absolute_import, print_function
from hfof import fof
from hfof.lib import get_cells
from numpy import unique
from utils import get_sim_pos
import numpy as np
        
def time_kdcount():
    """ Clustered n^3 points from real simulations """
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()

    print('Loading lib', file=log)
    cells = get_cells([(0.5, 0.2, 0.1)], 0.333, 3, 11, log)

#    for i, ngrid in [(1024, 8867), (128,1107), (256, 2217), (512, 4433)]:
    for i in [128, 128,256,512]:#,1024]:
        pos = get_sim_pos(i)-0.5

        rcut = 0.2/i

        domains = fof(pos, rcut, log=log)

        print('Number of unique domains {:,}'.format(len(unique(domains))),file=log)
        continue
        
        if False:
            boxsize = 1.0
            domains = fof_periodic(pos, boxsize, rcut, log=log)

            print('Number of unique domains {:,}'.format(len(unique(domains))),file=log)
        
            continue

        print('Trying with kdcount',file=log)

        d = dataset(pos)#, boxsize=1.0)
        f = kdfof(d, rcut)
        print('kdcount size {:,}'.format(len(f.length)), file=log)
        continue
        print('Trying with periodic kdcount',file=log)
        d = dataset(pos, boxsize=1.0)
        f = kdfof(d, rcut)
        print('Periodic kdcount Size {:,}'.format(len(f.length)), file=log)
    
def test_scaling():
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()

    print('Loading lib', file=log)
    cells = get_cells([(0.5, 0.2, 0.1)], 0.333, 3, log)

    i, ngrid = (128,1107)
    pos = get_sim_pos(i)

    rcut = (3**0.5)/ngrid
    
    boxsize = 1.0
    scale = 2.5
    domains = fof_periodic(pos, boxsize, rcut, log=log)
    domains_s = fof_periodic(pos*scale, boxsize*scale, rcut*scale, log=log)
    import numpy as np
    assert(np.all(np.equal(domains, domains_s)))

def test_1024():
    """ Why does the 1024 simn break for periodic? """
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()

    print('Loading lib', file=log)
    cells = get_cells([(0.5, 0.2, 0.1)], 0.333, 3, log)

#    for i, ngrid in [(1024, 8867), (128,1107), (256, 2217), (512, 4433)]:
    for i in [1024]:
        pos = get_sim_pos(i)

        rcut = 0.2/i

        boxsize = 1.0
        domains = fof_periodic(pos, boxsize, rcut, log=log)

        print('Number of unique domains {:,}'.format(len(unique(domains))),file=log)
        domains = fof(pos, rcut, log=log)

        print('Number of unique domains {:,}'.format(len(unique(domains))),file=log)

#        continue

        print('Trying with kdcount',file=log)

        d = dataset(pos)#, boxsize=1.0)
        f = kdfof(d, rcut)
        print('kdcount size {:,}'.format(len(f.length)), file=log)
        print('Trying with periodic kdcount',file=log)
        d = dataset(pos, boxsize=1.0)
        f = kdfof(d, rcut)
        print('Periodic kdcount Size {:,}'.format(len(f.length)), file=log)

def time_rcut():
    """ Timings with rcut """
    from lizard.log import VerboseTimingLog
    from kdcount.cluster import fof as kdfof, dataset
    log = VerboseTimingLog()

    print('Loading lib', file=log)
    cells = get_cells([(0.5, 0.2, 0.1)], 0.333, 3, 11, log)

    b = np.linspace(-2, 0.0, 5)
    i = 128
    pos = get_sim_pos(i)#  - 0.5

    rcut = 10.0**b / i
    from time import time
    from nbodyshop import nbodyshop_fof
    time_hfof = []
    time_kd = []
    time_kdcount = []

    for bi,r in zip(b,rcut):
        t0 = time()
        domains = fof(pos, r)
        hfof_dt = time()- t0
        time_hfof.append('%.3f'%hfof_dt)
        n_hfof = len(np.unique(domains)) 
        kdt, foft, ngrps = nbodyshop_fof(pos, r)
        assert(ngrps==n_hfof)

        time_kd.append(('(%.3f, %.3f)'%(kdt, foft)))


        t0 = time()
        kdcount_pos = dataset(pos)
        f = kdfof(kdcount_pos, r)
        kdc_dt = time() - t0

        n_kdcount = len(f.length)
        assert(ngrps==n_kdcount)
        time_kdcount.append('%.3f'%kdc_dt)

        print('Link length', bi, 'Time', hfof_dt, kdt+foft, kdc_dt,'Groups', n_hfof, ngrps, file=log)


    print('Time kd', ', '.join(time_kd), file=log)
    print('Time hfof', ', '.join(time_hfof), file=log)
    print('Time kdcount', ', '.join(time_kdcount), file=log)


if __name__=='__main__':
#    test_scaling()
#    time_kdcount()
    time_rcut()
#    test_1024()
