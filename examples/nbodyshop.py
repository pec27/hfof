"""
Stuff for wrapping nbodyshop
"""
from os import system
import subprocess
import numpy as np

# Path for Nbodyshop fof in double precision
fof_dbl = '/home/peter.creasey/nbodyshop/fof/fof_dbl/fof'

def nbodyshop_fof(pos, rcut, period=None):
    if period is None:
        run_str = ' '.join((fof_dbl,'-e %.12f -m 1 -d -v -o /dev/null'%rcut))
    else:
        run_str = ' '.join((fof_dbl, '-p %.12f'%period, '-e %.12f -m 1 -d -v -o /dev/null'%rcut))
        
    p = subprocess.Popen(run_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, stdin=subprocess.PIPE)

    pos_to_tipsy(pos, p.stdin)


    kd_setup, search,ngroups = 0.0, 0.0, 0
    for line in p.stdout.readlines():
        if 'FOF tree build' in line:
            kd_setup = float(line.split(' ')[3])
        elif 'FOF CPU TIME' in line:
            search = float(line.split(' ')[3])
        elif 'Number of groups:' in line:
            ngroups = int(line.split(':')[1])            
#        print line,

    retval = p.wait()
    return kd_setup, search, ngroups

def pos_to_tipsy(pos, f):

    n = pos.shape[0]

    time, nbodies, ndim, nsph, ndark, nstar, pad = (0.0, n, 3, 0,n,0, 0)
    
    np.array(time, dtype=np.float64).tofile(f)
    np.array((nbodies, ndim, nsph, ndark, nstar, pad), dtype=np.int32).tofile(f)
    
    

    """	
	if (!xdr_double(pxdrs,&ph->time)) return 0;
	if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
	if (!xdr_int(pxdrs,&ph->ndim)) return 0;
	if (!xdr_int(pxdrs,&ph->nsph)) return 0;
	if (!xdr_int(pxdrs,&ph->ndark)) return 0;
	if (!xdr_int(pxdrs,&ph->nstar)) return 0;
	if (!xdr_int(pxdrs,&pad)) return 0;
        """
    """
    struct dark_particle {
    Real mass;
    Real pos[MAXDIM];
    Real vel[MAXDIM];
    Real eps;
    Real phi ;
    } ;
    """

    data = np.zeros((n,9), dtype=np.float64)
    data[:,1:4] = pos
    
    data.tofile(f)



if __name__=='__main__':
    print nbodyshop_fof(128, 0.2/128, period=1.0)
