#include "geom.h"
#include <math.h>
#include <string.h>

static inline double d_ij(const Instance *I, int i, int j) {
    // dist is stored in I->dist (dense), 1-based indexing
    size_t n = (size_t)I->n;
    return I->dist[(size_t)i * (n + 1) + (size_t)j];
}

void geom_build_points_mat4(const Instance *I, uint64_t k, Vec3 *x_out) {
    const int n = I->n;

    // zero everything (optional, but keeps deterministic
    for (int i = 0; i <= n; i++)
        x_out[i] = (Vec3){0.0, 0.0, 0.0};

    // Base points (match Python)
    // x[1] = (0,0,0)
    double d12 = d_ij(I, 1, 2);
    double d23 = d_ij(I, 2, 3);

    x_out[2] = (Vec3){-d12, 0.0, 0.0};

    double th3 = I->theta[3];
    x_out[3] = (Vec3){-d12 + d23 * cos(th3), d23 * sin(th3), 0.0};

    // Build B2
    Mat4 B2, B3, B;
    mat4_identity(&B2);
    mat4_identity(&B3);
    mat4_identity(&B);

    // B2 = [[-1,0,0,-d12],[0,1,0,0],[0,0,-1,0],[0,0,0,1]]
    B2.a[0][0] = -1.0;
    B2.a[0][3] = -d12;
    B2.a[1][1] = 1.0;
    B2.a[2][2] = -1.0;
    B2.a[3][3] = 1.0;

    // B3 depends on theta[3] and d23
    double ct = cos(th3);
    double st = sin(th3);
    double di = d23;

    // B3 = [[-ct, -st, 0, -di*ct],
    //       [ st, -ct, 0,  di*st],
    //       [  0,   0, 1,      0],
    //       [  0,   0, 0,      1]]
    B3.a[0][0] = -ct;
    B3.a[0][1] = -st;
    B3.a[0][2] = 0.0;
    B3.a[0][3] = -di * ct;
    B3.a[1][0] = st;
    B3.a[1][1] = -ct;
    B3.a[1][2] = 0.0;
    B3.a[1][3] = di * st;
    B3.a[2][0] = 0.0;
    B3.a[2][1] = 0.0;
    B3.a[2][2] = 1.0;
    B3.a[2][3] = 0.0;
    B3.a[3][0] = 0.0;
    B3.a[3][1] = 0.0;
    B3.a[3][2] = 0.0;
    B3.a[3][3] = 1.0;

    // B = B2 * B3
    mat4_mul(&B2, &B3, &B);

    // Iterate atoms t=4..n
    for (int t = 4; t <= n; t++) {
        // Match Python bit convention:
        // if i & (1 << (n-t)) == 0 => sgn=0 => sw = (+)abs_sw
        // else => sw = (-)abs_sw
        int bit_index = n - t;
        uint64_t bit = (k >> bit_index) & 1ULL;

        double th = I->theta[t];
        double ct_t = cos(th);
        double st_t = sin(th);
        double di_t = d_ij(I, t - 1, t);

        double cw = I->cw[t];
        double abs_sw = I->abs_sw[t];
        double sw = bit ? -abs_sw : +abs_sw;

        Mat4 A, C;
        mat4_identity(&A);
        mat4_identity(&C);

        // A = [[-ct, -st, 0, -di*ct],
        //      [ st*cw, -ct*cw, -sw, di*st*cw],
        //      [ st*sw, -ct*sw,  cw, di*st*sw],
        //      [ 0, 0, 0, 1 ]]
        A.a[0][0] = -ct_t;
        A.a[0][1] = -st_t;
        A.a[0][2] = 0.0;
        A.a[0][3] = -di_t * ct_t;
        A.a[1][0] = st_t * cw;
        A.a[1][1] = -ct_t * cw;
        A.a[1][2] = -sw;
        A.a[1][3] = di_t * st_t * cw;
        A.a[2][0] = st_t * sw;
        A.a[2][1] = -ct_t * sw;
        A.a[2][2] = cw;
        A.a[2][3] = di_t * st_t * sw;
        A.a[3][0] = 0.0;
        A.a[3][1] = 0.0;
        A.a[3][2] = 0.0;
        A.a[3][3] = 1.0;

        // C = B * A
        mat4_mul(&B, &A, &C);

        // position is translation part of C
        x_out[t] = mat4_position(&C);

        // B = C
        mat4_copy(&B, &C);
    }
}
