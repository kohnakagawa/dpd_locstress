CXX = icpc

EIGEN_FLAGS = -DEIGEN_NO_DEBUG

OPENMP = -openmp
DEBUG = -O0 -g -DDEBUG
WARNINGS = -Wall -Wextra # -Wnon-virtual-dtor -Woverloaded-virtual 
RELEASE = -O3 -ipo -no-prec-div
STDCPP11 = -std=c++11

# CXXFLAGS = $(DEBUG)
CXXFLAGS = $(RELEASE)

CXXFLAGS += $(WARNINGS)
CXXFLAGS += $(OPENMP)
CXXFLAGS += $(EIGEN_FLAGS)
CXXFLAGS += $(STDCPP11)

OBJECTS = main bucket_sorter dpdsystem force_calculator observer time_evolver chem_manager
OBJECTS_DPD = $(addprefix ./src/,$(OBJECTS))
OBJECTS_DPD_SSE4 = $(addsuffix .sse4_o,$(OBJECTS_DPD))
OBJECTS_DPD_AVX  = $(addsuffix .avx_o,$(OBJECTS_DPD))
OBJECTS_DPD_AVX2 = $(addsuffix .avx2_o,$(OBJECTS_DPD))

LIBRARY = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread # for lapack

TARGET = shald_dpd_sse4.out shald_dpd_avx.out shald_dpd_avx2.out config_maker.out

all : $(TARGET)

.SUFFIXES :
.SUFFIXES : .cpp .sse4_o
.cpp.sse4_o:
	$(CXX) $(CXXFLAGS) -xSSE4.2 -c $< $(LIBRARY) -o $@

.SUFFIXES : .cpp .avx_o
.cpp.avx_o:
	$(CXX) $(CXXFLAGS) -xAVX -c $< $(LIBRARY) -o $@

.SUFFIXES : .cpp .avx2_o
.cpp.avx2_o:
	$(CXX) $(CXXFLAGS) -xCORE-AVX2 -c $< $(LIBRARY) -o $@

shald_dpd_sse4.out : $(OBJECTS_DPD_SSE4)
	$(CXX) $(CXXFLAGS) -xSSE4.2 $(OBJECTS_DPD_SSE4) $(LIBRARY) -o $@

shald_dpd_avx.out : $(OBJECTS_DPD_AVX)
	$(CXX) $(CXXFLAGS) -xAVX $(OBJECTS_DPD_AVX) $(LIBRARY) -o $@

shald_dpd_avx2.out : $(OBJECTS_DPD_AVX2)
	$(CXX) $(CXXFLAGS) -xCORE-AVX2 $(OBJECTS_DPD_AVX2) $(LIBRARY) -o $@

config_maker.out : ./src/config_maker.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

clean:
	rm -f $(OBJECTS_DPD_SSE4) $(OBJECTS_DPD_AVX) $(OBJECTS_DPD_AVX2) $(TARGET) core.* *~ ./src/*~
