CXX = g++
CXXFLAGS = -std=c++11 -O3  -Wall -mavx
SRCS = main.cpp\
	

OBJS = $(SRCS:.cpp=.o)

nbody : $(OBJS)
	$(CXX) -std=c++11 -mavx -o SpaceHub $(OBJS)

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) 
