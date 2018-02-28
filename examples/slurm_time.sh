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

#nbodyshop
python examples/time_tipsy.py /mainvol/peter.creasey/lizbench/data

# hash based fof
python examples/time_hfof.py /mainvol/peter.creasey/lizbench/data

# kd count
python examples/time_kdcount.py /mainvol/peter.creasey/lizbench/data







