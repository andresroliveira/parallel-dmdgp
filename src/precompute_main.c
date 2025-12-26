#include "instance.h"
#include <stdio.h>

static void dump_precompute(const Instance *I) {
    int n = I->n;

    printf("n=%d m=%d\n", I->n, I->m);

    printf("\n Bond lengths (k=2..n):\n");
    for (int k = 2; k <= n; k++) {
        printf("bond[%d]=%.12f\n", k, I->bond[k]);
    }

    printf("\nTheta (k=3..n):\n");
    for (int k = 3; k <= n; k++) {
        printf("theta[%d]=%.12f\n", k, I->theta[k]);
    }

    printf("\nCos(omega) and |sin(omega)| (k=4..n):\n");
    for (int k = 4; k <= n; k++) {
        printf("cw[%d]=%.12f  abs_sw[%d]=%.12f\n", k, I->cw[k], k,
               I->abs_sw[k]);
    }
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <instance_file>\n", argv[0]);
        return 1;
    }

    const char *path = argv[1];

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
        fprintf(stderr, "ERROR: precompute failed (numerical issue).\n");
        instance_free(&I);
        return 1;
    }

    dump_precompute(&I);
    instance_free(&I);
    return 0;
}
