from __future__ import absolute_import, print_function
from test import get_sim_pos
from lizard.log import VerboseTimingLog
from nbodyshop import nbodyshop_fof

def time(drt):
    """ Clustered n^3 points from real simulations """
    log = VerboseTimingLog()

    for i in [128, 128,256,512,1024]:
        pos = get_sim_pos(i, drt)

        rcut = 0.2/i


        print('Trying with nbodyshop fof',file=log)

        kd_setup, search, ngroups = nbodyshop_fof(pos, rcut)
        print('Times (%.3f, %.3f)'%(kd_setup, search), 'Groups {:,}'.format(ngroups), file=log)


if __name__=='__main__':
    
    from sys import argv
    drt = argv[1]
    time(drt)
