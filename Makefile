CXX := g++
CXXFLAGS := -ggdb -Wall -Wextra -Wno-missing-field-initializers -std=c++20 -DNDEBUG
OPTFLAGS := -Ofast -flto -funroll-loops -fno-signed-zeros -fno-trapping-math \
             -fassociative-math -freciprocal-math -fno-math-errno \
             -mtune=native -march=native -fno-math-errno

SRC_DIR := src
INC_DIR := include
BIN := build/diffuser

RAYLIB_DIR := ./raylib-macos/lib
LDFLAGS := -L$(RAYLIB_DIR) -Wl,-rpath,$(RAYLIB_DIR)
LDLIBS := -lraylib -lm

SRCS := $(wildcard $(SRC_DIR)/*.cpp)

all: $(BIN)

$(BIN): $(SRCS)
	@mkdir -p $(dir $(BIN))
	$(CXX) $(CXXFLAGS) $(OPTFLAGS) -I$(INC_DIR) $^ -o $@ $(LDFLAGS) $(LDLIBS)

clean:
	rm -rf $(BIN)

.PHONY: all clean