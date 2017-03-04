CXX = g++
CXXFLAGS = -O3 -fomit-frame-pointer -fopenmp -mavx2 -march=core-avx2 -lfftw3 -lfftw3_threads -fopenmp -Wall

all: fftw

fftw: fftw.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@

clean:
	rm -f *.out *.o
