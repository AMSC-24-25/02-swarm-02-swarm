CXX=g++
INCLUDE=-I./include
CXXFLAGS=-Wall -Wextra -Wpedantic -Werror -fsanitize=address
DEBUGFLAGS=-O0 -g
OPTFLAGS=-O3 -DNDEBUG -march=native -mtune=native

SRCS = $(wildcard src/*.cpp)
DEBUG_OBJS = $(patsubst src/%.cpp, build/debug-%.o, $(SRCS))
RELEASE_OBJS = $(patsubst src/%.cpp, build/release-%.o, $(SRCS))

all: debug

debug: $(DEBUG_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/main.cpp -o build/main $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/test.cpp -o build/test $^

release: $(RELEASE_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/main.cpp -o build/main $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/test.cpp -o build/main $^

build/debug-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c $(INCLUDE) -o $@ $^

build/release-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $(INCLUDE) -o $@ $^

format:
	clang-format --style=file -i src/*.cpp include/*.hpp


clean:
	rm -rf build/*
