CXX=g++
INCLUDE=-I./include
CXXFLAGS=-Wall -Wextra -Wpedantic -Werror
DEBUGFLAGS=-O0 -g -fopenmp
OPTFLAGS=-O3 -DNDEBUG -march=native -mtune=native -fopenmp

SRCS = $(wildcard src/*.cpp)
DEBUG_OBJS = $(patsubst src/%.cpp, build/debug-%.o, $(SRCS))
RELEASE_OBJS = $(patsubst src/%.cpp, build/release-%.o, $(SRCS))

all: debug

debug: $(DEBUG_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) -o build/main $^

release: $(RELEASE_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) -o build/main $^

build/debug-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c $(INCLUDE) -o $@ $^

build/release-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $(INCLUDE) -o $@ $^

format:
	clang-format --style=file -i src/*.cpp include/*.hpp

clean:
	rm -rf build/*
