CXX=g++ -Wall -pedantic
ARCH=-march=native
CXXFLAGS=-O3 -fno-exceptions
LDFLAGS=-Wl,--hash-style=gnu -Wl,--as-needed
all:
	${CXX} ${ARCH} ${CXXFLAGS} ${LDFLAGS} -o canon canon.cc
