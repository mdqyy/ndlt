/**
 * Use SVD to find null space for small matrix.
 *
 * Adopted from numerical recipes.
 * 
 * @blackball
 */

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#define ABS(v) ((v) > 0 ? (v) : -(v))
#define SQR(v) ((v)*(v))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) > (b) ? (b) : (a))
#define SIGN(a, b) ((b) >= 0 ? ((a) >= 0 ? (a) : -(a)) : ((a) >= 0 ? -(a) : (a)))

struct svdata_t {
        int m, n;
        double *u;  /* mxn */
        double *s;  /* 1xn */
        double *v;  /* nxn */
        double *rv; /* 1xn, buffer */       
};

struct svdata_t*
svdata_new(int m, int n) {
        const size_t sz = sizeof(struct svdata_t) + sizeof(double) * (m*n + n + n*n + n);
        struct svdata_t *sd = (struct svdata_t *)malloc(sz);
        sd->m = m;
        sd->n = n;
        sd->u = (double *)(sd + 1);
        sd->s = sd->u + m*n;
        sd->v = sd->s + n;
        sd->rv = sd->v + n*n;
        return sd; 
}

void 
svdata_free(struct svdata_t *sd) {
        free(sd);
}

static double
pythag(const double a, const double b) {
        const double absa = ABS(a), absb = ABS(b);
        return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
                (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

/* Given a MxN matrix A, where m >= n, find the nullspace of it.
 * The smallest eigenvalue -> lambda
 * And the crosspondent eigenvector -> nullv.
 * Return 0 if success, else -1.
 */
int
svd(struct svdata_t *sd) {
        const double eps = DBL_EPSILON;
        int flag,i,its,j,jj,k,l,nm;
        double anorm,c,f,g,h,s,scale,x,y,z;
        const int m = sd->m;
        const int n = sd->n;
        double *u = sd->u;
        double *w = sd->s;
        double *v = sd->v;
        double *rv1 = sd->rv;

        g = scale = anorm = 0.0;

        for (i=0;i<n;i++) {
                l=i+2;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i < m) {
                        for (k=i;k<m;k++) scale += ABS(u[k*n+i]);
                        if (scale != 0.0) {
                                for (k=i;k<m;k++) {
                                        u[k*n+i] /= scale;
                                        s += u[k*n+i]*u[k*n+i];
                                }
                                f=u[i*n+i];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                u[i*n+i]=f-g;
                                for (j=l-1;j<n;j++) {
                                        for (s=0.0,k=i;k<m;k++) s += u[k*n+i]*u[k*n+j];
                                        f=s/h;
                                        for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
                                }
                                for (k=i;k<m;k++) u[k*n+i] *= scale;
                        }
                }
                w[i]=scale *g;
                g=s=scale=0.0;
                if (i+1 <= m && i+1 != n) {
                        for (k=l-1;k<n;k++) scale += ABS(u[i*n+k]);
                        if (scale != 0.0) {
                                for (k=l-1;k<n;k++) {
                                        u[i*n+k] /= scale;
                                        s += u[i*n+k]*u[i*n+k];
                                }
                                f=u[i*n+(l-1)];
                                g = -SIGN(sqrt(s),f);
                                h=f*g-s;
                                u[i*n+(l-1)]=f-g;
                                for (k=l-1;k<n;k++) rv1[k]=u[i*n+k]/h;
                                for (j=l-1;j<m;j++) {
                                        for (s=0.0,k=l-1;k<n;k++) s += u[j*n+k]*u[i*n+k];
                                        for (k=l-1;k<n;k++) u[j*n+k] += s*rv1[k];
                                }
                                for (k=l-1;k<n;k++) u[i*n+k] *= scale;
                        }
                }
                anorm=MAX(anorm,(ABS(w[i])+ABS(rv1[i])));
        }
        for (i=n-1;i>=0;i--) {
                if (i < n-1) {
                        if (g != 0.0) {
                                for (j=l;j<n;j++)
                                        v[j*n+i]=(u[i*n+j]/u[i*n+l])/g;
                                for (j=l;j<n;j++) {
                                        for (s=0.0,k=l;k<n;k++) s += u[i*n+k]*v[k*n+j];
                                        for (k=l;k<n;k++) v[k*n+j] += s*v[k*n+i];
                                }
                        }
                        for (j=l;j<n;j++) v[i*n+j]=v[j*n+i]=0.0;
                }
                v[i*n+i]=1.0;
                g=rv1[i];
                l=i;
        }
        for (i=MIN(m,n)-1;i>=0;i--) {
                l=i+1;
                g=w[i];
                for (j=l;j<n;j++) u[i*n+j]=0.0;
                if (g != 0.0) {
                        g=1.0/g;
                        for (j=l;j<n;j++) {
                                for (s=0.0,k=l;k<m;k++) s += u[k*n+i]*u[k*n+j];
                                f=(s/u[i*n+i])*g;
                                for (k=i;k<m;k++) u[k*n+j] += f*u[k*n+i];
                        }
                        for (j=i;j<m;j++) u[j*n+i] *= g;
                } else for (j=i;j<m;j++) u[j*n+i]=0.0;
                ++u[i*n+i];
        }
        for (k=n-1;k>=0;k--) {
                for (its=0;its<30;its++) {
                        flag=1;
                        for (l=k;l>=0;l--) {
                                nm=l-1;
                                if (l == 0 || ABS(rv1[l]) <= eps*anorm) {
                                        flag=0;
                                        break;
                                }
                                if (ABS(w[nm]) <= eps*anorm) break;
                        }
                        if (flag) {
                                c=0.0;
                                s=1.0;
                                for (i=l;i<k+1;i++) {
                                        f=s*rv1[i];
                                        rv1[i]=c*rv1[i];
                                        if (ABS(f) <= eps*anorm) break;
                                        g=w[i];
                                        h=pythag(f,g);
                                        w[i]=h;
                                        h=1.0/h;
                                        c=g*h;
                                        s = -f*h;
                                        for (j=0;j<m;j++) {
                                                y=u[j*n+nm];
                                                z=u[j*n+i];
                                                u[j*n+nm]=y*c+z*s;
                                                u[j*n+i]=z*c-y*s;
                                        }
                                }
                        }
                        z=w[k];
                        if (l == k) {
                                if (z < 0.0) {
                                        w[k] = -z;
                                        for (j=0;j<n;j++) v[j*n+k] = -v[j*n+k];
                                }
                                break;
                        }
                        if (its == 29) {
                                return -1;
                        }
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=pythag(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                                i=j+1;
                                g=rv1[i];
                                y=w[i];
                                h=s*g;
                                g=c*g;
                                z=pythag(f,h);
                                rv1[j]=z;
                                c=f/z;
                                s=h/z;
                                f=x*c+g*s;
                                g=g*c-x*s;
                                h=y*s;
                                y *= c;
                                for (jj=0;jj<n;jj++) {
                                        x=v[jj*n+j];
                                        z=v[jj*n+i];
                                        v[jj*n+j]=x*c+z*s;
                                        v[jj*n+i]=z*c-x*s;
                                }
                                z=pythag(f,h);
                                w[j]=z;
                                if (z) {
                                        z=1.0/z;
                                        c=f*z;
                                        s=h*z;
                                }
                                f=c*g+s*y;
                                x=c*y-s*g;
                                for (jj=0;jj<m;jj++) {
                                        y=u[jj*n+j];
                                        z=u[jj*n+i];
                                        u[jj*n+j]=y*c+z*s;
                                        u[jj*n+i]=z*c-y*s;
                                }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                }
        }
        
        return 0;
}

#if 0

#define mprint(A, m, n) do{                                             \
                int i, j;                                               \
                for (i = 0; i < m; ++i) {                               \
                        for (j = 0; j < n; ++j) {                       \
                                printf("%lf,", *(A + i * n + j));       \
                        }                                               \
                        printf("\n");                                   \
                }                                                       \
                printf("\n");                                           \
        } while(0)

static void 
svdata_print(const struct svdata_t *sd) {
        mprint(sd->u, sd->m, sd->n);
        mprint(sd->s, 1, sd->n);
        mprint(sd->v, sd->n, sd->n);
}

static void 
test_case() {
        const double A[3*3] = {1,0,0,1,1,1,0,0,1};
        struct svdata_t *sd = svdata_new(3, 3);
        
        memcpy(sd->u, A, sizeof(double) * 9);
        
        if (svd(sd) != 0) {
                return ;
        }

        svdata_print(sd);
        svdata_free(sd);
}

int 
main(int argc, char *argv[]) {
        test_case();
        getchar();
        return 0;
}

#endif
