FLAGS = -std=c++11 -O3 -Wall -m64 -pthread -g
GMP   = -lgmpxx -lgmp 
FILES = src/Main.cpp src/Functions.hpp src/QuadraticSieve.hpp

all: bin/Main.out bin/Test.out

bin/Main.out: $(FILES)
	g++ $(FILES) -o bin/Main.out $(FLAGS) $(GMP)

bin/Main.exe: $(FILES)
	x86_64-w64-mingw32-g++ $(FILES) -o bin/Main.exe $(FLAGS) -static  $(GMP) 

bin/Test.out: $(FILES) src/Tests.cpp
	cxxtestgen --error-printer -o  .Tests.cxxtestrunner.cpp src/Tests.cpp
	g++ $(FLAGS) src/Functions.hpp .Tests.cxxtestrunner.cpp -o bin/Test.out $(GMP)