#include "instance.h"
#include <stdio.h>

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <instance_file>\n", argv[0]);
        return 1;
    }

    Instance I;
    if (!instance_load(argv[1], &I)) {
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

    printf("search: TODO (next step: implement g() + first-k OpenMP)\n");

    instance_free(&I);
    return 0;
}
