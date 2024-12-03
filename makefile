CXX=g++
INCLUDE=-I./include
CXXFLAGS=-Wall -Wextra -Wpedantic -Wshadow -Werror
DEBUGFLAGS=-O0 -g
OPTFLAGS=-O3 -DNDEBUG -march=native -mtune=native

SRCS = $(wildcard src/*.cpp)
OBJS = $(patsubst src/%.cpp, build/%.o, $(SRCS))

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o build/$@ $^

build/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $(INCLUDE) -o $@ $^

format:
	clang-format --style=file -i src/*.cpp include/*.hpp

clean:
	rm -rf build
