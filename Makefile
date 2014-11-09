FLAGS = -std=c++0x -O3 -Wall -m64 -pthread
GMP   = -lgmpxx -lgmp 
FILES = src/Main.cpp src/Functions.hpp

all: bin/Main.out bin/Test.out

bin/Main.out: $(FILES)
	g++ $(FLAGS) $(FILES) -o bin/Main.out $(GMP)

bin/Main.exe: $(FILES)
	x86_64-w64-mingw32-g++ $(FLAGS) $(FILES) -o bin/Main.exe $(GMP) -static

bin/Test.out: $(FILES) src/Tests.cpp
	cxxtestgen --error-printer -o  .Tests.cxxtestrunner.cpp src/Tests.cpp
	g++ $(FLAGS) src/Functions.hpp .Tests.cxxtestrunner.cpp -o bin/Test.out $(GMP)