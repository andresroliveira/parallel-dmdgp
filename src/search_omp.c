#include "geom.h"
#include "score.h"
#include "search.h"

#include <omp.h>
#include <stdatomic.h>
#include <stdlib.h>

SearchResult search_first_k_omp(const Instance *I, double delta) {
    SearchResult R = {0, 0, 0.0};

    const int n = I->n;
    const int m_bits = n - 3;
    if (m_bits <= 0 || m_bits >= 63)
        return R;

    const uint64_t total = 1ULL << m_bits;

    // Shared "found" flag with atomic visibility
    atomic_int found;
    atomic_init(&found, 0);

    // Shared result (written once under critical)
    uint64_t found_k = 0;
    double found_g = 0.0;

#pragma omp parallel
    {
        Vec3 *x = (Vec3 *)malloc(((size_t)n + 1) * sizeof(Vec3));
        if (!x) {
            // skip thread if allocation fails
        } else {
#pragma omp for schedule(dynamic, 1024)
            for (uint64_t k = 0; k < total; k++) {

                // If someone already found, stop doing work
                if (atomic_load_explicit(&found, memory_order_relaxed)) {
                    continue;
                }

                geom_build_points_mat4(I, k, x);
                double g = score_g_no_sqrt(I, x);

                if (g <= delta) {
// Publish exactly once
#pragma omp critical
                    {
                        if (!atomic_load_explicit(&found,
                                                  memory_order_relaxed)) {
                            found_k = k;
                            found_g = g;
                            atomic_store_explicit(&found, 1,
                                                  memory_order_relaxed);
                        }
                    }
                }
            }
            free(x);
        }
    }

    if (atomic_load_explicit(&found, memory_order_relaxed)) {
        R.found = 1;
        R.k = found_k;
        R.g = found_g;
    }

    return R;
}
