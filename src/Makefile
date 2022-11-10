
OUTPUT1 = cvplus1
OUTPUT2 = cvplus2

# Put C++ compiler 
#CXX = g++
CXX = clang++-9
CPP_FILES1 = main.cpp dtpr.cpp helpers.cpp parse_args.cpp read_inputs.cpp standardize.cpp
HPP_FILES = *.hpp 
CPP_FILES2 = main2.cpp dtpr.cpp helpers.cpp parse_args.cpp read_inputs.cpp


# Set compiler flags 
CXXFLAG = -std=c++11 -Wall -ggdb -O0 -fopenmp -larmadillo -lblas
all: $(OUTPUT1) $(OUTPUT2)

$(OUTPUT1): $(CPP_FILES1) $(HPP_FILES)
	$(CXX) $(CPP_FILES1) $(CXXFLAG) -v -o $(OUTPUT1) 

$(OUTPUT2): $(CPP_FILES2) $(HPP_FILES)
	$(CXX) $(CPP_FILES2) $(CXXFLAG) -v -o $(OUTPUT2) 

