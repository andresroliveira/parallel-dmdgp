#include "instance.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static inline size_t IDX(int n, int i, int j) {
    return (size_t)i * (size_t)(n + 1) + (size_t)j;
}

static inline int has_d(const Instance *I, int a, int b) {
    return I->has[IDX(I->n, a, b)] != 0;
}

static inline double get_d(const Instance *I, int a, int b) {
    return I->dist[IDX(I->n, a, b)];
}

static void set_d(Instance *I, int a, int b, double d) {
    I->dist[IDX(I->n, a, b)] = d;
    I->has[IDX(I->n, a, b)] = 1;
}

static int alloc_mats(Instance *I) {
    int n = I->n;
    size_t sz = (size_t)(n + 1) * (size_t)(n + 1);

    I->dist = (double *)calloc(sz, sizeof(double));
    I->has = (unsigned char *)calloc(sz, sizeof(unsigned char));
    I->theta = (double *)calloc((size_t)(n + 1), sizeof(double));
    I->cw = (double *)calloc((size_t)(n + 1), sizeof(double));
    I->abs_sw = (double *)calloc((size_t)(n + 1), sizeof(double));

    if (!I->dist || !I->has || !I->theta || !I->cw || !I->abs_sw)
        return 0;
    return 1;
}

int instance_load(const char *path, Instance *I) {
    memset(I, 0, sizeof(*I));

    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "ERROR: cannot open file: %s\n", path);
        return 0;
    }

    if (fscanf(f, "%d %d", &I->n, &I->m) != 2) {
        fprintf(stderr, "ERROR: failed to read 'n m' header\n");
        fclose(f);
        return 0;
    }
    if (I->n < 4) {
        fprintf(stderr, "ERROR: n must be >= 4 (got %d)\n", I->n);
        fclose(f);
        return 0;
    }
    if (I->m <= 0) {
        fprintf(stderr, "ERROR: m must be > 0 (got %d)\n", I->m);
        fclose(f);
        return 0;
    }

    if (!alloc_mats(I)) {
        fprintf(stderr, "ERROR: out of memory allocating matrices\n");
        fclose(f);
        return 0;
    }

    I->E = (Edge *)malloc((size_t)I->m * sizeof(Edge));
    if (!I->E) {
        fprintf(stderr, "ERROR: out of memory allocating edges\n");
        fclose(f);
        return 0;
    }

    for (int e = 0; e < I->m; e++) {
        int a, b;
        double d;
        if (fscanf(f, "%d %d %lf", &a, &b, &d) != 3) {
            fprintf(stderr, "ERROR: failed to read edge line %d\n", e + 1);
            fclose(f);
            return 0;
        }
        if (a < 1 || a > I->n || b < 1 || b > I->n || a == b) {
            fprintf(stderr, "ERROR: invalid edge (%d,%d) for n=%d\n", a, b,
                    I->n);
            fclose(f);
            return 0;
        }
        if (!(d > 0.0)) {
            fprintf(stderr,
                    "ERROR: non-positive distance for edge (%d,%d): %g\n", a, b,
                    d);
            fclose(f);
            return 0;
        }

        I->E[e].u = a;
        I->E[e].v = b;
        I->E[e].d = d;
        I->E[e].d2 = d * d;

        set_d(I, a, b, d);
        set_d(I, b, a, d);
    }

    fclose(f);
    return 1;
}

int instance_validate_dmdgp(const Instance *I) {
    int ok = 1;
    int n = I->n;

    // Needed for theta[k], k=3..n:
    for (int k = 3; k <= n; k++) {
        int a = k - 2, b = k - 1, c = k;
        if (!has_d(I, b, c) || !has_d(I, a, b) || !has_d(I, a, c)) {
            fprintf(stderr,
                    "DMDGP validation failed for theta[%d]: need d(%d,%d), "
                    "d(%d,%d), d(%d,%d)\n",
                    k, b, c, a, b, a, c);
            ok = 0;
        }
    }

    // Needed for cw[k], k=4..n:
    for (int k = 4; k <= n; k++) {
        int i = k - 3, j = k - 2, l = k; // (i,j,k-1,l) pattern in your code
        int k1 = k - 1;                  // k-1

        // rij=d(i,j), rjk=d(j,k1), rkl=d(k1,l), rjl=d(j,l), ril=d(i,l)
        if (!has_d(I, i, j) || !has_d(I, j, k1) || !has_d(I, k1, l) ||
            !has_d(I, j, l) || !has_d(I, i, l)) {
            fprintf(stderr,
                    "DMDGP validation failed for cw[%d]: need d(%d,%d), "
                    "d(%d,%d), d(%d,%d), d(%d,%d), d(%d,%d)\n",
                    k, i, j, j, k1, k1, l, j, l, i, l);
            ok = 0;
        }
    }

    return ok;
}

int instance_precompute(Instance *I) {
    int n = I->n;

    // theta[k] for k=3..n
    for (int k = 3; k <= n; k++) {
        double d_k1_k = get_d(I, k - 1, k);
        double d_k2_k1 = get_d(I, k - 2, k - 1);
        double d_k2_k = get_d(I, k - 2, k);

        double num = d_k1_k * d_k1_k + d_k2_k1 * d_k2_k1 - d_k2_k * d_k2_k;
        double den = 2.0 * d_k1_k * d_k2_k1;
        double c = num / den;

        // clamp for numerical safety
        if (c > 1.0)
            c = 1.0;
        if (c < -1.0)
            c = -1.0;

        I->theta[k] = acos(c);
    }

    // cw[k] for k=4..n (matches your Python)
    for (int k = 4; k <= n; k++) {
        double rij = get_d(I, k - 3, k - 2);
        double rjk = get_d(I, k - 2, k - 1);
        double rkl = get_d(I, k - 1, k);
        double rjl = get_d(I, k - 2, k);
        double ril = get_d(I, k - 3, k);

        double st = sin(I->theta[k - 1]);
        double ct = cos(I->theta[k - 1]);

        double inner = 4.0 * rjl * rjl * rjk * rjk -
                       (rjl * rjl + rjk * rjk - rkl * rkl) *
                           (rjl * rjl + rjk * rjk - rkl * rkl);
        if (inner <= 0.0 || fabs(st) < 1e-15 || fabs(rjk) < 1e-15) {
            fprintf(stderr,
                    "ERROR: numerical invalid denominator components at k=%d "
                    "(inner=%g, st=%g, rjk=%g)\n",
                    k, inner, st, rjk);
            return 0;
        }

        double num = rij * rij + rjl * rjl - ril * ril -
                     rij * (rjl * rjl + rjk * rjk - rkl * rkl) * ct / rjk;

        double den = rij * sqrt(inner) * st / rjk;

        double cw = num / den;

        // clamp to [-1,1] for safety (should be close)
        if (cw > 1.0)
            cw = 1.0;
        if (cw < -1.0)
            cw = -1.0;

        I->cw[k] = cw;

        double sw2 = 1.0 - cw * cw;
        if (sw2 < 0.0)
            sw2 = 0.0;
        I->abs_sw[k] = sqrt(sw2);
    }

    return 1;
}

void instance_free(Instance *I) {
    if (!I)
        return;
    free(I->E);
    free(I->dist);
    free(I->has);
    free(I->theta);
    free(I->cw);
    free(I->abs_sw);
    memset(I, 0, sizeof(*I));
}
