
CC=g++
#CXXFLAGS= -O3 -Wextra -std=c++11 -lm 
CXXFLAGS= -fopenmp -O3 -Wextra -lm -lgomp 


all: 
	$(CC) -o laseracv10_omp_basic_V0  laseracv10_omp_basic_V0.cpp   pcg_basic.h pcg_basic.c $(CXXFLAGS)

clean:
	rm -rf laseracv10_omp_basic_V0   *.o

