hfof
======

### Friends of Friends group finding via spatial hashing

hfof is an open-source friends-of-friends (FoF) group finder in 3d (periodic and non-periodic), based on the following paper (ZZZ arxiv link when public).

Once you have downloaded hfof you will probably want to do the following:

## Install

Build the C-functions with

```bash
make
```

Run the tests, e.g. with

```bash
nosetests tests/
```
(use `-v` verbose or `-vs` for really verbose)

Add hfof to your PYTHONPATH, perhaps in your .profile e.g.
```bash
PYTHONPATH="${PYTHONPATH}:${HOME}/codes/hfof"
export PYTHONPATH
```

## Examples
These probably need at least python 2.7 and a (non-ancient) numpy and scipy. 

### Labal a random set of particles with by their connectivity
    
In a python environment try
    
```python
from hfof import fof
from numpy.random import RandomState
rs = RandomState(seed=123)

ndim = 3 # 3-d
n = 32 # 32^3 points

pos = rs.rand(ndim*(n**ndim)).reshape((n**ndim,ndim))
pos[:,2] = 0.0 # Make 2d
r_cut = float(n)**(-1.5) # linking length

fof_labels = fof(pos, r_cut)

# color by number of particles in group
clrs = np.bincount(fof_labels)[fof_labels]
import pylab as pl
pl.scatter(pos[:,0],pos[:,1], c=clrs, s=4.1, edgecolors='none')
pl.show()
```

All of these are equivalent to running

```
python examples/example1.py
```
