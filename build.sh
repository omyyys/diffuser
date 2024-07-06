set -xe 
g++ -ggdb -fno-stack-protector -Wall -Wextra -Wno-missing-field-initializers -std=c++20 -O3 -o diffuser src/main.cpp solver/solver.cpp -L./raylib/lib -Wl,-rpath ./raylib/lib/ -lraylib -lm
