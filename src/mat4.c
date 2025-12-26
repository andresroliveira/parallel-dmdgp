#include "mat4.h"

void mat4_identity(Mat4 *M) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
            M->a[i][j] = 0.0;
        M->a[i][i] = 1.0;
    }
}

void mat4_copy(Mat4 *dst, const Mat4 *src) {
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            dst->a[i][j] = src->a[i][j];
}

void mat4_mul(const Mat4 *A, const Mat4 *B, Mat4 *C) {
    // C = A*B
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            double s = 0.0;
            for (int k = 0; k < 4; k++)
                s += A->a[i][k] * B->a[k][j];
            C->a[i][j] = s;
        }
    }
}

Vec3 mat4_position(const Mat4 *T) {
    Vec3 p;
    p.x = T->a[0][3];
    p.y = T->a[1][3];
    p.z = T->a[2][3];
    return p;
}
