#!/bin/bash 

#SBATCH --partition=mainq
#SBATCH -n 1 # tasks
#SBATCH --cpus-per-task=1 # raptor likes to add tasks per core, avoid
##SBATCH --mem-per-cpu=250000 # MB
#SBATCH --mem=250000 # MB
#SBATCH -t 48:00:00  # time to run
#SBATCH -J time_all # job name
#SBATCH -o out/time_all.o%j # append job number
#SBATCH --mail-user=p.e.creasey.00@gmail.com
#SBATCH --mail-type=end    # email me when the job finishes

# Load modules
module load openmpi-x86_64

# kd count
python examples/time_kdcount.py /mainvol/peter.creasey/lizbench/data

# hash based fof
python examples/time_hfof.py /mainvol/peter.creasey/lizbench/data

# The nbodyshop bits
/home/peter.creasey/nbodyshop/fof/fof_dbl/fof -e 0.0015625 -m 1 -d -v -o out/test128.fof < /mainvol/peter.creasey/lizbench/data/tipsy/pos_128.tipsy
/home/peter.creasey/nbodyshop/fof/fof_dbl/fof -e 0.00078125 -m 1 -d -v -o out/test256.fof < /mainvol/peter.creasey/lizbench/data/tipsy/pos_256.tipsy
/home/peter.creasey/nbodyshop/fof/fof_dbl/fof -e 0.000390625 -m 1 -d -v -o out/test512.fof < /mainvol/peter.creasey/lizbench/data/tipsy/pos_512.tipsy
/home/peter.creasey/nbodyshop/fof/fof_dbl/fof -e 0.0001953125 -m 1 -d -v -o out/test1024.fof < /mainvol/peter.creasey/lizbench/data/tipsy/pos_1024.tipsy





