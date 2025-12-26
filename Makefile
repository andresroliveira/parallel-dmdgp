CC := cc
CFLAGS := -O3 -march=native -std=c11 -Wall -Wextra -Iinclude -fopenmp
LDFLAGS := -lm

BUILD := build

COMMON_SRC := src/instance.c src/mat4.c src/geom_mat4.c src/score.c src/search_omp.c
COMMON_OBJ := $(COMMON_SRC:src/%.c=$(BUILD)/%.o)

all: $(BUILD)/precompute $(BUILD)/points $(BUILD)/search

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/%.o: src/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c -o $@ $<

$(BUILD)/precompute: $(BUILD)/precompute_main.o $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD)/points: $(BUILD)/points_main.o $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILD)/search: $(BUILD)/search_main.o $(COMMON_OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

debug: CFLAGS := -O0 -g -std=c11 -Wall -Wextra -Iinclude -fopenmp -fsanitize=address,undefined
debug: LDFLAGS := -lm -fsanitize=address,undefined
debug: clean all

clean:
	rm -rf $(BUILD)

.PHONY: all debug clean
