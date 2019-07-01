CXX=g++
OPENBLAS_DIR = $(HOME)/libs/OpenBLAS
INCLUDE=-I $(OPENBLAS_DIR)/include
CXXFLAGS=-std=c++11 -O3 $(INCLUDE) -fopenmp
SRCS=main.cpp dgemm.cpp sgemm.cpp
OBJS=$(SRCS:.cpp=.o)
LDFLAGS=-L $(OPENBLAS_DIR)/lib -lopenblas -fopenmp
TARGET=gemms

all:$(TARGET) $(OBJS)

.SUFFIXES: .cpp

$(TARGET):$(OBJS)
	$(CXX) -o $@ $(OBJS)  $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -rf $(TARGET) *.o
