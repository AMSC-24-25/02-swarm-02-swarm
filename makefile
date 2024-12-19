CXX=g++
INCLUDE=-I./include
CXXFLAGS=-Wall -Wextra -Wpedantic -Werror
DEBUGFLAGS=-O0 -g -fopenmp
OPTFLAGS=-O3 -DNDEBUG -march=native -mtune=native -fopenmp
BENCHFLAGS = -lbenchmark -lpthread

SRCS = $(wildcard src/*.cpp)
DEBUG_OBJS = $(patsubst src/%.cpp, build/debug-%.o, $(SRCS))
RELEASE_OBJS = $(patsubst src/%.cpp, build/release-%.o, $(SRCS))
BENCH_OBJS = $(patsubst src/%.cpp, build/debug-%.o, $(SRCS))

bench         = $(wildcard benchmark/*.cpp)
bench_targets = $(patsubst benchmark/%.cpp,build/benchmark/%,$(bench))

all: debug 

debug: $(DEBUG_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/main.cpp -o build/main $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/test.cpp -o build/test $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/TestRosenbrock.cpp -o build/TestRosenbrock $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(DEBUGFLAGS) exec/TestRastrigin.cpp -o build/TestRastrigin $^  # Aggiunta per TestRastrigin

release: $(RELEASE_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/main.cpp -o build/main $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/test.cpp -o build/test $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/TestRosenbrock.cpp -o build/TestRosenbrock $^
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(OPTFLAGS) exec/TestRastrigin.cpp -o build/TestRastrigin $^  # Aggiunta per TestRastrigin


benchmark: $(RELEASE_OBJS)
	$(CXX)  bench/benchmark_test.cpp $(INCLUDE) $(OPTFLAGS) -isystem benchmark/include   -Lbenchmark/build/src -lbenchmark -lpthread -o build/mm $^


build/debug-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) -c $(INCLUDE) -o $@ $^

build/release-%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -c $(INCLUDE) -o $@ $^

format:
	clang-format --style=file -i src/*.cpp include/*.hpp exec/*.cpp bench/*.cpp

clean:
	rm -rf build/*
