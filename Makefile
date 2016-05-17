CXX = g++
# CXX = icpc

EIGEN_FLAGS = -DEIGEN_NO_DEBUG

ifeq ($(CXX),icpc)
OPENMP = -openmp
endif
ifeq ($(CXX),g++)
OPENMP = -fopenmp
endif

DEBUG = -O0 -g -DDEBUG
WARNINGS = -Wall -Wextra -Wnon-virtual-dtor -Woverloaded-virtual 

ifeq ($(CXX),icpc)
RELEASE = -O3 -xHOST -ipo -no-prec-div
endif
ifeq ($(CXX),g++)
RELEASE = -O3 -ffast-math -funroll-loops
endif

STDCPP11 = -std=c++11

# CXXFLAGS = $(DEBUG)
CXXFLAGS = $(RELEASE)

CXXFLAGS += $(WARNINGS)
CXXFLAGS += $(OPENMP)
CXXFLAGS += $(EIGEN_FLAGS)
CXXFLAGS += $(STDCPP11)

OBJECTS = main.o bucket_sorter.o dpdsystem.o force_calculator.o observer.o time_evolver.o chem_manager.o
OBJECTS_DPD = $(addprefix ./src/,$(OBJECTS))

TARGET = shald_dpd.out config_maker.out

all : $(TARGET)

.SUFFIXES :
.SUFFIXES : .cpp .o
.cpp.o:
	$(CXX) $(CXXFLAGS) -c $< $(LIBRARY) -o $@

shald_dpd.out : $(OBJECTS_DPD)
	$(CXX) $(CXXFLAGS) $(OBJECTS_DPD) $(LIBRARY) -o $@

config_maker.out : ./src/config_maker.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS_DPD) $(TARGET) core.* *~ ./src/*~
