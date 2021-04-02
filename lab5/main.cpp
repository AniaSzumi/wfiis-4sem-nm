#include <iostream>
#include "/usr/local/include/gsl/gsl_matrix.h"
#include "/usr/local/include/gsl/gsl_vector.h"
#include <cmath>
#define abs(X) ((X)>0? (X):-(X))
#define n 7
#define IT_MAX 12

void print_matrix(gsl_matrix *A);
gsl_vector* matrix_vector_multiplication(gsl_matrix *A, gsl_vector *x);
double scalar_multiplication(gsl_vector *a, gsl_vector *b);
void vector_sum(gsl_vector *x, gsl_vector *a, gsl_vector *b, double alpha);
void save_matrix(gsl_matrix *A, FILE *fp);
void matrix_multiplication(gsl_matrix *X, gsl_matrix *A, gsl_matrix *B);
void matrix_transpose(gsl_matrix *X_inv, gsl_matrix *X);

int main(){
    gsl_matrix *A = gsl_matrix_calloc(n,n);
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            gsl_matrix_set(A,i,j,1.0/(sqrt(2+abs(i-j))));
        }
    }
    gsl_matrix *X = gsl_matrix_calloc(n,n);
    gsl_vector *xi = gsl_vector_calloc(n);
    gsl_vector *xi1 = gsl_vector_calloc(n);
    gsl_vector *xj = gsl_vector_calloc(n);
    double lambda = 0;
    double x_norm = 0;
    FILE *fp_lambda;
    fp_lambda = fopen("lambda.txt", "w");
    // fprintf(fp_lambda, "Wartosci wlasne: \n");

    for(int k=0; k<n; k++){
        for(int i=0; i<n; i++){
            gsl_vector_set(xi,i,1);
        }
        for(int i=0; i<IT_MAX; i++){
            fprintf(fp_lambda, "%d %g\n", i, lambda);
            xi1 = matrix_vector_multiplication(A, xi);
            for(int j=0; j<k; j++){
                for(int a=0; a<n; a++){
                    gsl_vector_set(xj,a,gsl_matrix_get(X,a,j));
                }
                vector_sum(xi1,xi1,xj,-scalar_multiplication(xi1,xj));
            }
            lambda = scalar_multiplication(xi1,xi)/scalar_multiplication(xi,xi);
            x_norm = sqrt(scalar_multiplication(xi1,xi1));
            for(int a=0; a<n; a++){
                gsl_vector_set(xi,a,gsl_vector_get(xi1,a)/x_norm);

            }
        }
        for(int a=0; a<n; a++){
            gsl_matrix_set(X,a,k,gsl_vector_get(xi,a));
        }
        fprintf(fp_lambda, "\n\n");

    }
    FILE *fp_vector;
    fp_vector = fopen("vector.txt", "w");

    fprintf(fp_vector, "\nWektory wlasne: \n");
    save_matrix(X, fp_vector);
    gsl_matrix *X_inv = gsl_matrix_calloc(n,n);
    gsl_matrix *D = gsl_matrix_calloc(n,n);
    gsl_matrix *help = gsl_matrix_calloc(n,n);
    matrix_transpose(X_inv,X);
    matrix_multiplication(help,X_inv,A);
    matrix_multiplication(D,help,X);
    fprintf(fp_vector, "\nD: \n");
    save_matrix(D, fp_vector);

    fclose(fp_lambda);
    fclose(fp_vector);
    gsl_matrix_free(A);
    gsl_matrix_free(X);
    gsl_matrix_free(X_inv);
    gsl_matrix_free(D);
    gsl_vector_free(xi);
    gsl_vector_free(xi1);
    gsl_vector_free(xj);
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

gsl_vector* matrix_vector_multiplication(gsl_matrix *A, gsl_vector *x){
    gsl_vector *y = gsl_vector_calloc(n);
    int jmin, jmax;
    for(int i=0; i<n; i++){
        gsl_vector_set(y,i,0);
        for(int j=0; j<n; j++){
            gsl_vector_set(y,i,gsl_vector_get(y,i)+gsl_matrix_get(A,i,j)*gsl_vector_get(x,j));
        }
    }
    return y;
}

double scalar_multiplication(gsl_vector *a, gsl_vector *b){
    double sum = 0.0;
    for(int i=0; i<n; i++){
        sum += gsl_vector_get(a,i) * gsl_vector_get(b,i);
    }
    return sum;
}

void vector_sum(gsl_vector *x, gsl_vector *a, gsl_vector *b, double alpha){
    for(int i=0; i<n; i++){
        gsl_vector_set(x,i,(gsl_vector_get(a,i)+alpha*gsl_vector_get(b,i)));
    }
}

void save_matrix(gsl_matrix *A, FILE *fp){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            fprintf(fp, " %g ", gsl_matrix_get(A,i,j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

void matrix_multiplication(gsl_matrix *X, gsl_matrix *A, gsl_matrix *B){
    for (int i=0; i<n ; i++){
        for (int j=0; j<n; j++){
            double sum = 0.0; 
            for (int k=0; k<n ; k++){
                sum += gsl_matrix_get(A,i,k) * gsl_matrix_get(B,k,j);
            }
            gsl_matrix_set(X,i,j,sum);
        }
    }
}

void matrix_transpose(gsl_matrix *X_inv, gsl_matrix *X){
    for (int i=0; i<n ; i++){
        for (int j=0; j<n; j++){
            gsl_matrix_set(X_inv,i,j,gsl_matrix_get(X,j,i));
        }
    }
}