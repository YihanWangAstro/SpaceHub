CXX = g++
CXXFLAGS = -std=c++11 -O3 -mtune=native -Wall -mavx2 -lprofiler -ltcmalloc
SRCS = main.cpp\
	

OBJS = $(SRCS:.cpp=.o)

nbody : $(OBJS)
	$(CXX) -std=c++11 -o SpaceX $(OBJS)

$(OBJS) : $(INCL)
.PHONY : clean	
clean:
	-rm $(OBJS) 
