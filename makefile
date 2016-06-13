CC = g++
FLAGS = -O3
all:
	$(CC) $(FLAGS) source/repeatFinder.cpp -o source/repeatFinder
