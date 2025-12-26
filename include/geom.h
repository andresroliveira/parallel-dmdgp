#ifndef GEOM_H
#define GEOM_H

#include <stdint.h>
#include "instance.h"
#include "mat4.h"

// Build points x[1..n] for a given k, using the 4x4 matrix method 
// Convention: bit index for atom t is (n - t).
void geom_build_points_mat4(const Instance *I, uint64_t k, Vec3 *x_out);

#endif // GEOM_H
