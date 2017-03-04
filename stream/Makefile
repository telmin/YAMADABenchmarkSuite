CXX = g++
CFLAGS = -O2 -std=c++11 -Wall -fopenmp

all: stream_cxx.out

stream_cxx.out: stream.cpp
	$(CXX) $(CFLAGS) $^ -o $@

clean:
	rm -f *.out *.o
