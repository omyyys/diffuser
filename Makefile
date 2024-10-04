CC=g++
CPPFLAGS=-ggdb -Wall -Wextra -Wno-missing-field-initializers -std=c++20
OPT=-Ofast -funroll-loops -fno-signed-zeros -fno-trapping-math -mtune=native -march=native
CPPLIBS=-L./raylib/lib -Wl,-rpath ./raylib/lib/ -lraylib -lm

all: diffuser

diffuser: src/main.cpp solver/solver.cpp solver/solver.hpp
	$(CC) $(CPPFLAGS) $(OPT) -o diffuser src/main.cpp solver/solver.cpp $(CPPLIBS) 
