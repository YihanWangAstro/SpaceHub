CXX=c++
CXXFLAGS=-Wall -std=c++17 -O3 -pthread

all: $(patsubst %.cpp, %.out, $(wildcard *.cpp))

%.out: %.cpp Makefile
	$(CXX) $(CXXFLAGS) $< -o $(@:.out=)

clean: $(patsubst %.cpp, %.clean, $(wildcard *.cpp))

%.clean:
	rm -f $(@:.clean=)

