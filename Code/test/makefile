CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags
LDFLAGS := -L /home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/lib -L ~/lib # Library paths
INCLUDES := -I /home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/include -I ~/include # Include paths
LIBS := -ldivsufsort -lsdsl # Libraries

TARGETS := chk # Change this to your desired executable name 
SOURCES := chk.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

chk: chk.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		rm -f $(TARGETS)