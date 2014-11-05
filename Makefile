all: bin/Main.out

bin/Main.out: src/Main.cpp
	g++ -O3 -std=c++0x -g -Wall src/Main.cpp -o bin/Main.out -lgmp