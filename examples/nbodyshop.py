"""
Stuff for wrapping nbodyshop
"""
from os import system
import subprocess

# Path for Nbodyshop fof in double precision
fof_dbl = '/home/peter.creasey/nbodyshop/fof/fof_dbl/fof'
# where the tipsy data is TODO make tipsy files into pipe
data_drt = '/mainvol/peter.creasey/lizbench/data/tipsy/pos_%d.tipsy'

def nbodyshop_fof(size, rcut, period=None):
    if period is None:
        run_str = ' '.join((fof_dbl,'-e %.12f -m 1 -d -v -o /dev/null'%rcut, '< '+data_drt%size))
    else:
        run_str = ' '.join((fof_dbl, '-p %.12f'%period, '-e %.12f -m 1 -d -v -o /dev/null'%rcut, '< '+data_drt%size))        
    p = subprocess.Popen(run_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

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


if __name__=='__main__':
    print nbodyshop_fof(128, 0.2/128, period=1.0)
