#ifndef MAT4_H
#define MAT4_H

typedef struct { double x, y, z; } Vec3;

typedef struct {
    double a[4][4];
} Mat4;

void mat4_identity(Mat4 *M);
void mat4_copy(Mat4 *dst, const Mat4 *src);
void mat4_mul(const Mat4 *A, const Mat4 *B, Mat4 *C);   // C = A*B
Vec3 mat4_position(const Mat4 *T);                      // (T[0][3],T[1][3],T[2][3])

#endif // MAT4_H
