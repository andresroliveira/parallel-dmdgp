#include "score.h"

double score_g_no_sqrt(const Instance *I, const Vec3 *x) {
    double s = 0.0;

    for (int e = 0; e < I->m; e++) {
        int a = I->E[e].u;
        int b = I->E[e].v;

        double dx = x[a].x - x[b].x;
        double dy = x[a].y - x[b].y;
        double dz = x[a].z - x[b].z;

        double dist2 = dx * dx + dy * dy + dz * dz;
        double diff = dist2 - I->E[e].d2;

        s += diff * diff;
    }

    return s;
}
