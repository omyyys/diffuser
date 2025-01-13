CC=g++
CPPFLAGS=-ggdb -Wall -Wextra -Wno-missing-field-initializers -std=c++20
OPT=-Ofast -funroll-loops -fno-signed-zeros -fno-trapping-math -mtune=native -march=native
CPPLIBS=-L./raylib-macos/lib -Wl,-rpath ./raylib-macos/lib/ -lraylib -lm

all: diffuser

diffuser: src/main.cpp src/solver.cpp include/solver.hpp
	$(CC) $(CPPFLAGS) $(OPT) -o diffuser src/main.cpp src/solver.cpp $(CPPLIBS) 
