from __future__ import print_function
from hfof import fof, example_data
import numpy as np

def test_fof():
    """ Test friends-of-friends """

    pos = example_data.get_pos() # (32768,3) positions in [0,1]
    pos[:,2] *= 0.0 # project (squash) z-dimension
    r_cut = 0.004 # linking length

    fof_labels = fof(pos, r_cut, boxsize=1.0) # integer labels from 0...
    print('Number of labels', len(np.unique(fof_labels)), 'minimum', min(fof_labels), 'maximum', max(fof_labels))

    # color by number of particles in group
    clrs = np.bincount(fof_labels)[fof_labels]
    import pylab as pl
    pl.scatter(pos[:,0],pos[:,1], c=np.power(clrs,0.3), s=1.0, edgecolors='none')
#    pl.xlim(0,1)
#    pl.ylim(0,1)
#    pl.axis('off')
#    pl.savefig('figs/fof.png', dpi=50, bbox_inches='tight')
    pl.show()
if __name__=='__main__':
    test_fof()
