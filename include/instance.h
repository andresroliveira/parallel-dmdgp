#ifndef INSTANCE_H
#define INSTANCE_H

#include <stddef.h>

typedef struct {
    int u, v;        // 1-based
    double d;        // distance
    double d2;       // distance^2 (precompute)
} Edge;

typedef struct {
    int n;           // number of vertices
    int m;           // number of edges read
    Edge *E;         // array length m

    // Dense distance matrix + presence flags (1..n)
    double *dist;    // size (n+1)*(n+1)
    unsigned char *has; // same size; 1 if dist present

    // Precomputed arrays (1..n)
    double *theta;    // theta[k] for k>=3
    double *ctheta;   // cos(theta[k])
    double *stheta;   // sin(theta[k])
    double *cw;       // cw[k] = cos(omega_k) for k>=4
    double *abs_sw;   // abs_sw[k] = |sin(omega_k)| for k>=4
    double *bond;     // Bond lengths (2..n): bond[k] = d[k-1][k]
} Instance;

// Load file, allocate matrices, store edges/distances.
// Returns 1 on success, 0 on failure (prints error to stderr).
int instance_load(const char *path, Instance *I);

// Validate DMDGP-required distances for the vertex order 1..n.
// Returns 1 if valid, 0 if invalid (prints the missing requirements).
int instance_validate_dmdgp(const Instance *I);

// Precompute theta, cw, abs_sw.
// Requires instance_validate_dmdgp() to be true.
// Returns 1 on success, 0 on numerical failure (prints details).
int instance_precompute(Instance *I);

// Free memory.
void instance_free(Instance *I);

#endif // INSTANCE_H
