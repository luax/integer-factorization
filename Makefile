FLAGS = -std=c++0x -O3 -Wall 
GMP   = -lgmpxx -lgmp

all: bin/Main.out bin/Test.out

bin/Main.out:  src/Main.cpp src/Functions.hpp
	g++ $(FLAGS) src/Main.cpp src/Functions.hpp -o bin/Main.out $(GMP)

bin/Test.out: src/Tests.cpp src/Functions.hpp
	cxxtestgen --error-printer -o  .Tests.cxxtestrunner.cpp src/Tests.cpp
	g++ $(FLAGS) src/Functions.hpp .Tests.cxxtestrunner.cpp -o bin/Test.out $(GMP)
