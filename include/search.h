#ifndef SEARCH_H
#define SEARCH_H

#include <stdint.h>
#include "instance.h"

typedef struct {
    int found;          // 1 if found
    uint64_t k;         // smallest valid k
    double g;           // g(h(k))
} SearchResult;

// Find the smallest k in [0, 2^(n-3)) such that score <= delta.
// Returns found=1 if exists, else found=0.
SearchResult search_first_k_omp(const Instance *I, double delta);

#endif // SEARCH_H

