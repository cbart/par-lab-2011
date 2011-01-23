#!/bin/sh
ulimit -s unlimited
time mpiexec -n 16 canon
