CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags
LDFLAGS := -L ../libdivsufsort/build/lib -L ~/lib # Library paths
INCLUDES := -I ../libdivsufsort/build/include -I ~/include # Include paths
LIBS := -ldivsufsort -lsdsl # Libraries

TARGETS := summary # Change this to your desired executable name 
SOURCES := summary.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

summary: summary.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		rm -f $(TARGETS)