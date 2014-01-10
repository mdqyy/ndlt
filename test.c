/* Normalized DLT for Homography estimation */
#include "point.h"
#include "svd.h"

#include "homo-fit.c"
#include "normalize.c"

static void
test(void) {
        struct point_t pts1[4] = {};
        struct point_t pts1[4] = {};
        struct svdata_t *svdata = svdata_new(8, 9);
        double T1[3], T2[3], H[9];
        
        get_conditioner(pts1, 4, T1);
        get_conditioner(pts2, 4, T2);
        normalize_points(T1, pts1, 4); // in place
        normalize_points(T2, pts2, 4);
        homo_fit(pts1, pts2, 4, svdata, H);
        decondition(T1, T2, H);
        svdata_free(svdata);
}

int
main(int argc, char *argv[]) {
        test();
        return 0;
}
