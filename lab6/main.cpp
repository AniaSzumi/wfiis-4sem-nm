#include <iostream>
#include "/usr/local/include/gsl/gsl_matrix.h"
#include "/usr/local/include/gsl/gsl_vector.h"
#include <cmath>
#define abs(X) ((X)>0? (X):-(X))
#define N 5

double licz_r(gsl_vector* a, gsl_vector* b, int n, double xj);

int main(){
    gsl_vector *a = gsl_vector_calloc(N+1);
    gsl_vector *b = gsl_vector_calloc(N+1);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector_set(a,0,240.0);
    gsl_vector_set(a,1,-196.0);
    gsl_vector_set(a,2,-92.0);
    gsl_vector_set(a,3,33.0);
    gsl_vector_set(a,4,14.0);
    gsl_vector_set(a,5,1.0);

    int n;
    double x0, x1, Rj0, Rj1;
    FILE *fp;
    fp = fopen("out.txt", "w");
    fprintf(fp, "L i           xj           Rj          Rj'\n");
    fprintf(fp, "------------------------------------------\n\n");

    for(int L=1; L<=N; L++){
        n = N-L+1;
        x0 = 0.0;
        for(int i=1; i<30; i++){
            Rj0 = licz_r(a,b,n,x0);
            Rj1 = licz_r(b,c,n-1,x0);
            x1 = x0 - Rj0/Rj1;
            fprintf(fp, "%d %d %12g %12g %12g\n", L, i, x1, Rj0, Rj1);
            if(abs(x1-x0) < 1.0e-7) break;
            x0 = x1;
        }
        for(int i=0; i<n; i++){
            gsl_vector_set(a,i,gsl_vector_get(b,i));
        }
        fprintf(fp, "\n\n");
    }


    fclose(fp);
    gsl_vector_free(a);
    gsl_vector_free(b);
    gsl_vector_free(c);
    return 0;
}


double licz_r(gsl_vector* a, gsl_vector* b, int n, double xj){
    gsl_vector_set(b,n,0.0);
    for(int k=n-1; k>=0; k--){
        gsl_vector_set(b,k,gsl_vector_get(a,k+1)+xj*gsl_vector_get(b,k+1));
    }
    return (gsl_vector_get(a,0) + xj * gsl_vector_get(b,0) );
}
