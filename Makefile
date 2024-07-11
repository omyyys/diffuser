CC=g++
CPPFLAGS=-ggdb -fno-stack-protector -Wall -Wextra -Wno-missing-field-initializers -Ofast -funroll-loops -fno-signed-zeros -fno-trapping-math -mtune=native -march=native  -std=c++20 
CPPLIBS=-L./raylib/lib -Wl,-rpath ./raylib/lib/ -lraylib -lm

all: diffuser

diffuser: src/main.cpp solver/solver.cpp solver/solver.hpp
	$(CC) $(CPPFLAGS) -o diffuser src/main.cpp solver/solver.cpp $(CPPLIBS) 
