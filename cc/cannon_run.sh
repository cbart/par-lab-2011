#!/bin/sh
ulimit -s unlimited
time mpiexec -n 4 cannon
