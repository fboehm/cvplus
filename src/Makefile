
OUTPUTD = cvplus

# Put C++ complier 
#CXX = g++
CXX = clang++-9
CPP_FILES = *.cpp
HPP_FILES = *.hpp 



# Set complier flags 
CXXFLAG = -std=c++11 -Wall -ggdb -O0 -fopenmp -larmadillo -lblas
all: $(OUTPUTD)
$(OUTPUTD): $(CPP_FILES) $(HPP_FILES)
	$(CXX) $(CPP_FILES) $(CXXFLAG) -v -o $(OUTPUTD) 
