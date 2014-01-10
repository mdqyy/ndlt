#ifndef SVD_H
#define SVD_H

struct svdata_t {
        int m, n;
        double *u;  /* mxn */
        double *s;  /* 1xn */
        double *v;  /* nxn */
        double *rv; /* 1xn, buffer */       
};

struct svdata_t svdata_new(int m, int n);
void svdata_free(struct svdata_t *svdata);
int svd(struct svdata_t *svdata);

#endif
