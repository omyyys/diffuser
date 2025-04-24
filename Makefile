CXX := clang++
CXXFLAGS := -Wall -Wextra -Wno-missing-field-initializers -std=c++20
OPTFLAGS := -O3 -flto=thin -ffast-math -fstrict-aliasing -fomit-frame-pointer \
             -funroll-loops -fvectorize -fslp-vectorize \
             -march=native -mtune=native \
             -fstrict-enums -fstrict-vtable-pointers \
             -fno-exceptions -fno-rtti \
             -ffp-contract=fast -fno-rounding-math \
             -mllvm -inline-threshold=1000 \
             -mllvm -enable-unroll-and-jam=true

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