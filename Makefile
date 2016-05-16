#CXX = g++
CXX = icpc
#CXX = mpiFCCpx
#CXX = FCCpx

ifeq ($(CXX),icpc)
MER_FLAGS   = -DDSFMT_MEXP=1279 -DHAVE_SSE2 
#VEC_REPORT  = -vec_report3
endif
ifeq ($(CXX),g++)
MER_FLAGS   = -DDSFMT_MEXP=1279 -msse2  -DHAVE_SSE2 
endif
ifeq ($(CXX),FCCpx)
MER_FLAGS   = -DDSFMT_MEXP=1279
endif
ifeq ($(CXX),mpiFCCpx)
MER_FLAGS   = -DDSFMT_MEXP=1279
endif

EIGEN_FLAGS = -DEIGEN_NO_DEBUG

ifeq ($(CXX),icpc)
OPENMP = -openmp
endif
ifeq ($(CXX),mpiFCCpx)
OPENMP = -Kopenmp
endif
ifeq ($(CXX),FCCpx)
OPENMP = -Kopenmp
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp
endif

DEBUG = -O0 -Wall -Wnon-virtual-dtor -Woverloaded-virtual -g -DDEBUG
PROFILE = -pg -O2

ifeq ($(CXX),icpc)
RELEASE = -O3 -xHOST -ipo -no-prec-div
endif
ifeq ($(CXX),mpiFCCpx)
RELEASE = -Kfast -Nsrc
endif
ifeq ($(CXX),FCCpx)
RELEASE = -Kfast -Nsrc
endif
ifeq ($(CXX),g++)
RELEASE = -O3 
endif

STDCPP11 = -std=c++11

CFLAGS = $(DEBUG)
#CFLAGS = $(RELEASE)
#CFLAGS = $(PROFILE)
CFLAGS += $(MER_FLAGS)
CFLAGS += $(OPENMP)
CFLAGS += $(EIGEN_FLAGS)
CFLAGS += $(STDCPP11)

#LIBRARY = -lstdc++ #This line should be added using intel compiler option "-fast"
#LIBRARY += -lmpi   #This line should be added using SYSTEM B MPI application

OBJECTS = main.o bucket_sorter.o dSFMT.o dpdsystem.o force_calculator.o initializer.o observer.o time_evolver.o chem_manager.o

#OBJECTS = main.o bucket_sorter.o dpdsystem.o force_calculator.o initializer.o observer.o time_evolver.o chem_manager.o

TARGET = shald_dpd.out
all:$(TARGET)

.SUFFIXES:
.SUFFIXES: .cpp .o
.cpp.o:
	$(CXX) $(CFLAGS) -c $< $(LIBRARY) -o $@
.SUFFIXES: .c .o
.c.o:
	$(CXX) $(CFLAGS) -c $< $(LIBRARY) -o $@

$(TARGET): $(OBJECTS)
	$(CXX) $(CFLAGS) $(VEC_REPORT) $(OBJECTS) $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS) $(TARGET) core.* *~
