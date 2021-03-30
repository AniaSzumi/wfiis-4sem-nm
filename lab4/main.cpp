#include <iostream>
#include <gsl/gsl_eigen.h>
#include <cmath>
#define n 200
#define L 10.0
#define deltax static_cast<double>(L/(n+1))

double rho(double x, int alpha){
    return 1.0 + static_cast<double>(4 * alpha) * x * x;
}

double deltak(int i, int j){
    return i == j ? 1.0 : 0.0;
}

double xi(int i){
    return -L/2 + deltax*(i+1);
}
void print_matrix(gsl_matrix *A);

void solve(gsl_matrix *A, gsl_matrix *B, int alpha, gsl_vector *eval, gsl_matrix *evec, gsl_eigen_gensymmv_workspace *w, FILE *f_eval, FILE *f_evec){
    double N = 1;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            gsl_matrix_set(A,i,j,(-deltak(i,j+1) + 2*deltak(i,j) - deltak(i,j-1))/(deltax*deltax));
            gsl_matrix_set(B,i,j,rho(xi(i),alpha)*deltak(i,j)/N);
        }
    }
    gsl_eigen_gensymmv(A, B, eval, evec, w);
    gsl_eigen_gensymmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);
    fprintf(f_eval, "\n%d ", alpha);
    for(int i=0; i<6; i++){
        fprintf(f_eval, "%g ", sqrt(gsl_vector_get(eval,i)));
    }
    if(alpha == 0 || alpha == 100){
        for(int i=0; i<n; i++){
            fprintf(f_evec, "\n%g ", xi(i));
            for(int j=0; j<6; j++){
                fprintf(f_evec, "%g ", gsl_matrix_get(evec,i,j));
            }
        }
        fprintf(f_evec,"\n\n");
    }

}


int main(){
    int alpha = 0;
    gsl_matrix *A = gsl_matrix_calloc(n,n);
    gsl_matrix *B = gsl_matrix_calloc(n,n);
    gsl_vector *eval = gsl_vector_calloc(n);
    gsl_matrix *evec = gsl_matrix_calloc(n,n);
    gsl_eigen_gensymmv_workspace *w = gsl_eigen_gensymmv_alloc(n);

    FILE *f_eval, *f_evec;
    f_eval = fopen("eval.txt", "w");
    f_evec = fopen("evec.txt", "w");
    for(;alpha <= 100; alpha += 2){
        solve(A,B,alpha,eval,evec,w,f_eval,f_evec);
    }
    fclose(f_eval);
    fclose(f_evec);

    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);

    return 0;
}

void print_matrix(gsl_matrix *A){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            std::cout<<" "<<gsl_matrix_get(A,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}