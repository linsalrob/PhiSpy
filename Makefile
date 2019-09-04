CC = g++
FLAGS = -O3
all:
	$(CC) $(FLAGS) src/repeatFinder.cpp -o bin/repeatFinder
