#include "geom.h"
#include "score.h"
#include "search.h"

#include <limits.h>
#include <omp.h>
#include <stdlib.h>

SearchResult search_first_k_omp(const Instance *I, double delta) {
    SearchResult R;
    R.found = 0;
    R.k = ULLONG_MAX;
    R.g = 0.0;

    const int n = I->n;
    const int m_bits = n - 3;
    if (m_bits <= 0 || m_bits >= 63) {
        // for now we assume small n; extend later if needed
        return R;
    }

    const uint64_t total = 1ULL << m_bits;

    // shared best_k
    uint64_t best_k = ULLONG_MAX;
    double best_g = 0.0;

#pragma omp parallel
    {
        Vec3 *x = (Vec3 *)calloc((size_t)n + 1, sizeof(Vec3));
        if (!x) {
            // If allocation fails, just skip this thread.
            // (You could also abort, but this is fine for now.)
        } else {
#pragma omp for schedule(static)
            for (uint64_t k = 0; k < total; k++) {

                // Early skip: if we already have a best_k smaller than current
                // k, no need to continue. Note: benign race (may do extra
                // work); correctness preserved.
                if (best_k != ULLONG_MAX && k >= best_k)
                    continue;

                geom_build_points_mat4(I, k, x);
                double g = score_g_no_sqrt(I, x);

                if (g <= delta) {
#pragma omp critical
                    {
                        if (k < best_k) {
                            best_k = k;
                            best_g = g;
                        }
                    }
                }
            }
            free(x);
        }
    }

    if (best_k != ULLONG_MAX) {
        R.found = 1;
        R.k = best_k;
        R.g = best_g;
    }
    return R;
}
