# Project-specific paths
LIB_DIR = ./lib
INCLUDE_DIR = ./include
SRC_DIR = ./src
BUILD_DIR = ./build

# Compiler and Flags
CC = gcc
CFLAGS = -I$(INCLUDE_DIR) -Wall -O3 -fPIC
LDFLAGS = -L$(LIB_DIR) -Wl,-rpath,$(LIB_DIR) -lrebound

# OS-specific settings
LIBS = -lm
ifeq ($(OS), Linux)
    LIBS += -lrt
endif
ifeq ($(OS), Windows_NT)
    CFLAGS += /D_GNU_SOURCE /Ox /fp:precise
endif

# Targets
EXE = $(BUILD_DIR)/systemRebound
SRC = $(SRC_DIR)/systemRebound.c

.PHONY: all clean

all: $(EXE)

$(EXE): $(SRC)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) $(LIBS)

clean:
	@rm -rf $(BUILD_DIR)
