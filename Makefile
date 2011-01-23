all:
	mpicc -Wall -pedantic -O3 -o canon_c canon.c
	g++ -Wall -pedantic -O0 -g -o canon canon.cc
