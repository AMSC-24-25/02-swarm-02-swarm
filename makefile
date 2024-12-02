CXX=g++
INCLUDE = -I ./include

SRCS  = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, build/objects/%.o, $(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(INCLUDE) -o build/executors/$@ $^

build/objects%.o: src/%.cpp
	$(CXX) -c $(INCLUDE) -o $@ $^

clean:
	rm -f *.o main