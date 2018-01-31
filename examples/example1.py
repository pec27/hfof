from hfof import fof
import numpy as np

def test_fof2d(n=32):
    """ Test friends-of-friends """
    from numpy.random import RandomState
    rs = RandomState(seed=123)

    ndim = 3

    pos = rs.rand(ndim*(n**ndim)).reshape((n**ndim,ndim))
    pos[:,2] = 0.0 # flatten to 2d
    r_cut = 1.05*(float(n)**(-1.5)) # linking length

    fof_labels = fof(pos, r_cut)
    print 'Number of labels', len(np.unique(fof_labels)), 'minimum', min(fof_labels), 'maximum', max(fof_labels)

    # color by number of particles in group
    clrs = np.bincount(fof_labels)[fof_labels]
    import pylab as pl
    pl.scatter(pos[:,0],pos[:,1], c=clrs, s=4.1, edgecolors='none')
    pl.show()

if __name__=='__main__':
    test_fof2d()
