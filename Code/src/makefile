CXX := g++ # Compiler
CXXFLAGS := -std=c++17 -O3 # Compiler flags
LDFLAGS := -L /home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/lib -L ~/lib # Library paths
INCLUDES := -I /home/daanish/projects/DiploidGenomeAssembly/libdivsufsort/build/include -I ~/include # Include paths
LIBS := -ldivsufsort -lsdsl # Libraries

# TARGETS := intradouble intratriple interdouble # Change this to your desired executable name 
# SOURCES := intradouble_c.cpp intratriple.cpp interdouble.cpp # Assumes your source files have a .cpp extension

TARGETS := intratriple # Change this to your desired executable name 
# SOURCES := intratriple.cpp # Assumes your source files have a .cpp extension

all: $(TARGETS)

# intradouble: intradouble_c.cpp
# 		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

intratriple: intratriple.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)
		
# interdouble: interdouble.cpp
# 		$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		rm -f $(TARGETS)