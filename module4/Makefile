CXX := g++

#######################################################
# Problem size definition
NFREQ ?= 16384

#######################################################
# Set NUMA_ALLOCATOR = Y to use NUMA aware memory allocator 
NUMA_ALLOCATOR = N

# Convert NUMA_ALLOCATOR=Y to a preprocessor define
ifeq ($(NUMA_ALLOCATOR),Y)
NUMA_FLAG = -DNUMA_ALLOCATOR
else
NUMA_FLAG =
endif

#######################################################
# Optimization flags are chosen as the last definition.
# Comment out using "#" at the begining of the line or rearrange according to your needs.
#
# Fastest executable (-ffast-math removes checking for NaNs and other things)
OPT=-O3 -ffast-math

# Add thread checking to code
#OPT=-O1 -fsanitize=thread,undefined,bounds -fPIE -pie

# Add memory / address error checking to code
#OPT=-O1 -fsanitize=leak,address,undefined,bounds -fPIE -pie

# Faster compilation time
#OPT=-O1

CXXFLAGS := $(OPT) -Wall -march=native -g -std=c++23 -DNFREQ=$(NFREQ) $(NUMA_FLAG)

default: seq mp1 mp2 mp3 mp4 mp5

seq: seismogram_seq.cpp
	$(CXX) seismogram_seq.cpp $(CXXFLAGS) -o seq

# mp: seismogram_omp1.cpp
# 	$(CXX) seismogram_omp1.cpp $(CXXFLAGS) -fopenmp -o mp

mp1: seismogram_omp1.cpp
	$(CXX) seismogram_omp1.cpp $(CXXFLAGS) -fopenmp -o mp1

mp2: seismogram_omp2.cpp
	$(CXX) seismogram_omp2.cpp $(CXXFLAGS) -fopenmp -o mp2

mp3: seismogram_omp3.cpp
	$(CXX) seismogram_omp3.cpp $(CXXFLAGS) -fopenmp -o mp3

mp4: seismogram_omp4.cpp
	$(CXX) seismogram_omp4.cpp $(CXXFLAGS) -fopenmp -o mp4

mp5: seismogram_omp5.cpp
	$(CXX) seismogram_omp5.cpp $(CXXFLAGS) -fopenmp -o mp5

# mpone: seismogram_omp1.cpp 
#  	$(CXX) seismogram_omp1.cpp $(CXXFLAGS) -fopenmp -o mpone

clean:
	rm -fr seq mp1 mp2 mp3 mp4 mp5
