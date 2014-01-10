#include "point.h"
#include "svd.h"

static int
homo_fit(const struct point_t *pt1, const struct point_t *pt2, int n, struct svdata_t *usv, double model[9]) {
        int i, mini;
        double *u = usv->u, mins;
        const double *s = usv->s;
        const double *v = usv->v;

        if (n != 4) {
                // only support 4 points now.
                retunr -1;
        }
        
#define _ASSIGN(i)                                                      \
        u[i*18+0] = 0;                      u[i*18+1] = 0;                      u[i*18+2] = 0; \
        u[i*18+3] = -pt1[i].x;              u[i*18+4] = -pt1[i].y;              u[i*18+5] = -1; \
        u[i*18+6] = pt2[i].y * pt1[i].x;    u[i*18+7] = pt2[i].y * pt1[i].y;    u[i*18+8] = pt2[i].y; \
        u[i*18+9] = pt1[i].x;               u[i*18+10] = pt1[i].y;              u[i*18+11] = 1; \
        u[i*18+12] = 0;                     u[i*18+13] = 0;                     u[i*18+14] = 0; \
        u[i*18+15] = -pt2[i].x * pt1[i].x;  u[i*18+16] = -pt2[i].x * pt1[i].y;  u[i*18+17] = -pt2[i].x
 
        _ASSIGN(0);
        _ASSIGN(1);
        _ASSIGN(2);
        _ASSIGN(3);
        
#undef _ASSIGN
 
        if (!svd(usv)) {
                return -1;
        }
 
        mini = 0;
        mins = s[0];
        for (i = 1; i < 9; ++i) {
                if (s[i] < mins) {
                        mins = s[i];
                        mini = i;
                }
        }
        for (i = 0; i < 9; ++i) {
                model[i] = v[mini * 9 + i];
        }
        
        return 0;
}
