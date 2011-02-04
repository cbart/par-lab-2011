CXX=mpiCC -Wall -Wpointer-arith -pedantic
ARCH=-march=native
CXXFLAGS=-O3 -fno-exceptions
LDFLAGS=-Wl,--hash-style=gnu -Wl,--as-needed -Wl,--rpath
INCLUDES=-I/usr/include/mpi -I/usr/include/boost/mpi -I/usr/lib/openmpi/include/openmpi/ompi/mpi
LIBRARIES=-L/usr/lib
LIBS=-lboost_mpi-mt -lboost_mpi -lboost_serialization-mt -lboost_serialization

all:
	${CXX} cannon.cc ${ARCH} ${CXXFLAGS} ${LDFLAGS} ${INCLUDES} ${LIBS} -o cannon

clean:
	@rm -f cannon

.PHONY: all clean
