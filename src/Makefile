
OUTPUTD = cvplus

# Set the path of library
ARMALIB = /net/mulan/home/yasheng/Cpp/arma/lib

# Put C++ complier 
#CXX = g++
CXX = clang++-9
CPP_FILES = *.cpp
HPP_FILES = *.hpp 



# Set complier flags 
CXXFLAG = -std=c++11 -Wall
all: $(OUTPUTD)
$(OUTPUTD): $(CPP_FILES) $(HPP_FILES)
	$(CXX) $(CPP_FILES) $(CXXFLAG) -L $(ARMALIB) -v -o $(OUTPUTD) 
