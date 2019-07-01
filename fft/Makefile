CXX = g++
CXXFLAGS = -O3 -fomit-frame-pointer -fopenmp -mavx2 -march=core-avx2 -Wall
LDFLAGS = -lfftw3 -lfftw3_threads

all: fftw

fftw: fftw.cpp
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *.out *.o fftw
