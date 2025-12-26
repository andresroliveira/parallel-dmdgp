// #include "geom.h"
#include "instance.h"
#include "search.h"

#include <stdio.h>
#include <stdlib.h>

static void usage(const char *prog) {
    fprintf(stderr, "usage: %s <instance_file> <delta>\n", prog);
    fprintf(stderr, "example: %s data/7_18.in 1e-4\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char *path = argv[1];
    double delta = strtod(argv[2], NULL);

    Instance I;
    if (!instance_load(path, &I)) {
        instance_free(&I);
        return 1;
    }
    if (!instance_validate_dmdgp(&I)) {
        fprintf(stderr, "ERROR: instance is not a DMDGP for vertex order 1..n. "
                        "Aborting.\n");
        instance_free(&I);
        return 1;
    }
    if (!instance_precompute(&I)) {
        fprintf(stderr, "ERROR: precompute failed.\n");
        instance_free(&I);
        return 1;
    }

    SearchResult R = search_first_k_omp(&I, delta);

    if (!R.found) {
        printf("NO SOLUTION: no k with g <= %.12g\n", delta);
        instance_free(&I);
        return 0;
    }

    printf("FOUND: k=%llu  g=%.12g  (delta=%.12g)\n", (unsigned long long)R.k,
           R.g, delta);

    // print points for the found k (helps debugging)
    // Vec3 *x = calloc((size_t)I.n + 1, sizeof(Vec3));
    // geom_build_points_mat4(&I, R.k, x);
    // for (int i = 1; i <= I.n; i++)
    //     printf("%d %.12f %.12f %.12f\n", i, x[i].x, x[i].y, x[i].z);
    // free(x);

    instance_free(&I);
    return 0;
}
