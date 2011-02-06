#!/bin/sh
ulimit -s unlimited
cat  | uniq > /home/users/cbart/par_lab_2011/nodes
export BOOST_LIBS=/home/users/cbart/lib
export LD_LIBRARY_PATH=:
export LD_RUN_PATH=:
time mpiexec --mca btl self,openib -n 64 --machinefile /home/users/cbart/par_lab_2011/nodes /home/users/cbart/par_lab_2011/cc/cannon
