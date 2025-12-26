#include "geom.h"
#include "score.h"
#include "search.h"

#include <omp.h>
#include <stdlib.h>

SearchResult search_first_k_omp(const Instance *I, double delta) {
    SearchResult R;
    R.found = 0;
    R.k = 0;
    R.g = 0.0;

    const int n = I->n;
    const int m_bits = n - 3;
    if (m_bits <= 0 || m_bits >= 63)
        return R;

    const uint64_t total = 1ULL << m_bits;

    // Shared flag and result
    int found = 0;
    uint64_t found_k = 0;
    double found_g = 0.0;

#pragma omp parallel
    {
        Vec3 *x = (Vec3 *)malloc(((size_t)n + 1) * sizeof(Vec3));
        if (!x) {
            // if allocation fails, do nothing
        } else {
#pragma omp for schedule(dynamic, 1024)
            for (uint64_t k = 0; k < total; k++) {

                // Early exit check (benign race; OK)
                if (found)
                    continue;

                geom_build_points_mat4(I, k, x);
                double g = score_g_no_sqrt(I, x);

                if (g <= delta) {
// Try to publish result once
#pragma omp critical
                    {
                        if (!found) {
                            found = 1;
                            found_k = k;
                            found_g = g;
                        }
                    }
                }
            }
            free(x);
        }
    }

    if (found) {
        R.found = 1;
        R.k = found_k;
        R.g = found_g;
    }
    return R;
}
