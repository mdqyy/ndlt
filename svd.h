#ifndef SVD_H
#define SVD_H

struct svdata_t;
struct svdata_t svdata_new(int m, int n);
void svdata_free(struct svdata_t *svdata);
int svd(struct svdata_t *svdata);

#endif
