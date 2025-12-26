CC := cc
CFLAGS := -O3 -march=native -std=c11 -Wall -Wextra -Iinclude
LDFLAGS := -lm

SRC := src/main.c src/instance.c
OBJ := $(SRC:.c=.o)

all: precompute

precompute: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LDFLAGS)

debug: CFLAGS := -O0 -g -std=c11 -Wall -Wextra -Iinclude -fsanitize=address,undefined
debug: LDFLAGS := -lm -fsanitize=address,undefined
debug: clean precompute

clean:
	rm -f $(OBJ) precompute
