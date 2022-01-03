#!/bin/bash

#SBATCH --job-name=t-hurdles
#SBATCH --time=07-00:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --partition=general-bionic
#SBATCH --mail-type=ALL

######################
# Comments
######################

# max time as of 2021 12 is 7 days 7-00:00:00


######################
# Begin work section #
######################



matlab-2018a -nodesktop -nodisplay -nosplash -r \
    "runstartup; disp(ls); many_estimates('$1',$2,$3); quit"    


