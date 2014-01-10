/* for 2D points <x, y> */
#include "point.h"
#include <math.h>
#include <float.h>

/*
          | 1  0  -x  |
  T = s * | 0  1  -y  |
          | 0  0  1/s |

  P' = T * P
  
 */
static void
get_conditioner(struct point_t *pts, int n, double T[3]) {
        double mx = 0, my = 0, ms = 0;
        int i = 0;
        for (; i < n; ++i) {
                mx += pts[i].x[0];
                my += pts[i].x[1];
        }

        mx /= (double)n;
        my /= (double)n;

        for (i = 0; i < n; ++i) {
                ms += sqrt((pts[i].x[0] - mx) * (pts[i].x[0] - mx) + (pts[i].x[1] - my) * (pts[i].x[1] - my));
        }

        if (ms == 0) {
                ms += DBL_EPLISON;
        }
        
        ms = sqrt(2.0) * n / ms;

        T[0] = -mx;
        T[1] = -my;
        T[2] = ms;
}

static void
normalize_points(const double T[3], struct point_t *pts, int n) {
        const double x = T[0], y = T[1], s = T[2];
        int i = 0;
        for (; i < n; ++i) {
                pts[i].x[0] = s * (pts[i].x[0] - x);
                pts[i].x[1] = s * (pts[i].x[0] - y);
        }
}

static void
decondition(const double T1[3], const double T2[3], double H[9]) {
        const double x1 = -T1[0], y1 = -T1[1], s1 = 1.0 / T1[2]; // inv(T1)
        const double x2 =  T2[0], y2 =  T2[1], s2 = T2[2];
        const double s12 = s1 * s2;
        /* T1^-1 * H * T2 */
        const double temp[9] = {
                H[0] + x1 * H[6], H[1] + x1 * H[7], H[2] + x1 * H[8],
                H[3] + y1 * H[6], H[4] + y1 * H[7], H[5] + y1 * H[8],
                H[6] / s1, H[7] / s1, H[8] / s1,
        }
        
        H[0] = temp[0] * s12;
        H[1] = temp[1] * s12;
        H[2] = (temp[0] * x2 + temp[1] * y2) * s12 + temp[2] * s1;
        H[3] = temp[3] * s12;
        H[4] = temp[4] * s12;
        H[5] = (temp[3] * x2 + temp[4] * y2) * s12 + temp[5] * s1;
        H[6] = temp[6] * s12;
        H[7] = temp[7] * s12;
        H[8] = (temp[6] * x2 + temp[7] * y2) * s12 + temp[8] * s1;
}
