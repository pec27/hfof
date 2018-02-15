from __future__ import absolute_import, print_function
from test import get_sim_pos
from lizard.log import VerboseTimingLog
from kdcount.cluster import fof as kdfof, dataset

def time(drt):
    """ Clustered n^3 points from real simulations """
    
    
    log = VerboseTimingLog()

    for i in [128, 128,256,512,1024]:
        pos = get_sim_pos(i, drt)

        rcut = 0.2/i


        print('Trying with kdcount',file=log)

        d = dataset(pos)#, boxsize=1.0)
        f = kdfof(d, rcut)
        print('kdcount size {:,}'.format(len(f.length)), file=log)


if __name__=='__main__':
    
    from sys import argv
    drt = argv[1]
    time(drt)
