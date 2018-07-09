CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -mavx -DKAHAN_SUMMATION
SRCS = main.cpp\
	

OBJS = $(SRCS:.cpp=.o)

nbody : $(OBJS)
	$(CXX) -std=c++11 -o SpaceX $(OBJS)

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) 
