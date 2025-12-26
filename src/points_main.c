#include "geom.h"
#include "instance.h"
#include <stdio.h>
#include <stdlib.h>

static void print_points(const Instance *I, const Vec3 *x) {
    for (int i = 1; i <= I->n; i++) {
        printf("%d %.12f %.12f %.12f\n", i, x[i].x, x[i].y, x[i].z);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <instance_file> <k_decimal>\n", argv[0]);
        return 1;
    }

    const char *path = argv[1];
    uint64_t k = (uint64_t)strtoull(argv[2], NULL, 10);

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

    Vec3 *x = (Vec3 *)calloc((size_t)I.n + 1, sizeof(Vec3));
    if (!x) {
        fprintf(stderr, "ERROR: out of memory\n");
        instance_free(&I);
        return 1;
    }

    geom_build_points_mat4(&I, k, x);
    print_points(&I, x);

    free(x);
    instance_free(&I);
    return 0;
}
