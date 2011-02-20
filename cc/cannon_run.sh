#!/bin/sh
ulimit -s unlimited
time mpiexec -n $1 cannon
