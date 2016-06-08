CXX = g++

DEBUG = -O0 -g -DDEBUG
WARNINGS = -Wall -Wextra -Wnon-virtual-dtor -Woverloaded-virtual 

RELEASE = -O3

STDCPP11 = -std=c++11

CXXFLAGS = $(DEBUG)
# CXXFLAGS = $(RELEASE)

CXXFLAGS += $(WARNINGS)
CXXFLAGS += $(STDCPP11)

LIBRARY = -llapacke -llapack

TARGET = do_cosine do_general do_dihedral rp_check tri_force

all : $(TARGET)

do_cosine : cosine_bend.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

do_general : general_bend.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

do_dihedral : dihedral.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

rp_check : check_rp.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

tri_force : tri_force_pov.cpp
	$(CXX) $(CXXFLAGS) $< $(LIBRARY) -o $@

clean:
	rm -f $(TARGET) core.* *~
