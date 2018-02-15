from __future__ import absolute_import, print_function
from test import get_sim_pos
from lizard.log import VerboseTimingLog
from hfof import fof, fof_periodic

def time(drt):
    """ Clustered n^3 points from real simulations """
    
    log = VerboseTimingLog()

    # Do smallest twice to average
    for i in [128, 128,256,512,1024]:
        pos = get_sim_pos(i, drt)

        rcut = 0.2/i

        domains = fof(pos, rcut, log=log)


if __name__=='__main__':
    
    from sys import argv
    drt = argv[1]
    time(drt)
