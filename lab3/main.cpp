#include <iostream>
#include "/usr/local/include/gsl/gsl_matrix.h"
#include "/usr/local/include/gsl/gsl_vector.h"
#include "../nr/nrutil.h"
#include "../nr/nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
#include "../nr/gaussj.c"
#include <cmath>
#include <time.h>

#define N 1000
#define m 5
#define abs(X) ((X)>0? (X):-(X))
#define max(X,Y) ((X)>(Y)? (X):(Y))
#define min(X,Y) ((X)<(Y)? (X):(Y))

gsl_vector* matrix_vector_multiplication(gsl_matrix *A, gsl_vector *x);
double scalar_multiplication(gsl_vector *a, gsl_vector *b);
void vector_sum(gsl_vector *a, gsl_vector *b, gsl_vector *x, double alpha);
void save_matrix(gsl_matrix *A, FILE *fp);
void print_matrix(gsl_matrix *A);
void print_vector(gsl_vector *v);

int main(){
    clock_t start = clock();
    //ALLOCATING MEMORY
    gsl_matrix *A = gsl_matrix_calloc(N,N);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            if(abs(i-j)>m){
                gsl_matrix_set(A,i,j,0);
            }
            else{
                gsl_matrix_set(A,i,j,1.0/(1+abs(i-j)));
            }
        }
    }   

    gsl_vector *b = gsl_vector_calloc(N);
    for(int i=0; i<N; i++){
        gsl_vector_set(b,i,static_cast<double>(i)+1.0);
    }
    gsl_vector *x = gsl_vector_calloc(N);
    for(int i=0; i<N; i++){
        gsl_vector_set(x,i,0);
    }
    gsl_vector *r = gsl_vector_calloc(N);
    gsl_vector *v = gsl_vector_calloc(N);

    vector_sum(r,b,matrix_vector_multiplication(A,x),-1.0);
    vector_sum(v,b,matrix_vector_multiplication(A,x),-1.0);

    double alpha, beta, scal_r = scalar_multiplication(r,r), scal_x = scalar_multiplication(x,x);
    int k=0;
    gsl_vector *Av = gsl_vector_calloc(N);
    FILE *fp;
    fp = fopen("out.txt", "w");
    fprintf(fp, "  k              rr           alpha            beta              xx\n");
    fprintf(fp, "%3d    %12g              -               -     %12g\n", k++, sqrt(scal_r), sqrt(scal_x));
    //ITERATING
    while(sqrt(scal_r) > 1e-6){
        scal_r = scalar_multiplication(r,r);
        scal_x = scalar_multiplication(x,x);
        Av = matrix_vector_multiplication(A,v);
        alpha = scal_r/scalar_multiplication(v,Av);
        vector_sum(x,x,v,alpha);
        vector_sum(r,r,Av,-alpha);
        beta = scalar_multiplication(r,r)/scal_r;
        vector_sum(v,r,v,beta);
        fprintf(fp, "%3d    %12g    %12g    %12g    %12g\n", k, sqrt(scalar_multiplication(r,r)), alpha, beta, sqrt(scalar_multiplication(x,x)));
        k++;
    }
    clock_t end = clock();
    double time = (double)(end - start) / CLOCKS_PER_SEC;
    std::cout << "Metoda sprzezonych gradientow: " << time << std::endl;
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_vector_free(r);
    gsl_vector_free(Av);


    clock_t startj = clock();
    float **Aj, **bj;
    Aj = matrix(1,N,1,N);
    bj = matrix(1,N,1,1);
    for(int i=1; i<=N; i++){
        for(int j=1; j<=N; j++){
            if(abs(i-j)>m){
                Aj[i][j] = 0.0;
            }
            else{
                Aj[i][j] = 1.0/(1+abs(i-j));
            }
        }
    }   

    for(int i=1; i<=N; i++){
        bj[i][1] = static_cast<double>(i)+1.0;
    }

    gaussj(Aj, N, bj, 1);
    clock_t endj = clock();
    double timej = (double)(endj - startj) / CLOCKS_PER_SEC;
    std::cout << "Metoda eliminacji zupelnej: " << timej << std::endl;
}


gsl_vector* matrix_vector_multiplication(gsl_matrix *A, gsl_vector *x){
    gsl_vector *y = gsl_vector_calloc(N);
    int jmin, jmax;
    for(int i=0; i<N; i++){
        jmin = max(0,i-m);
        jmax = min(i+m,N-1);
        gsl_vector_set(y,i,0);
        for(int j=jmin; j<=jmax; j++){
            gsl_vector_set(y,i,gsl_vector_get(y,i)+gsl_matrix_get(A,i,j)*gsl_vector_get(x,j));
        }
    }
    return y;
}

double scalar_multiplication(gsl_vector *a, gsl_vector *b){
    double sum = 0.0;
    for(int i=0; i<N; i++){
        sum += gsl_vector_get(a,i) * gsl_vector_get(b,i);
    }
    return sum;
}

void vector_sum(gsl_vector *x, gsl_vector *a, gsl_vector *b, double alpha){
    for(int i=0; i<N; i++){
        gsl_vector_set(x,i,(gsl_vector_get(a,i)+alpha*gsl_vector_get(b,i)));
    }
}


void save_matrix(gsl_matrix *A, FILE *fp){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fprintf(fp, " %g ", gsl_matrix_get(A,i,j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

void print_matrix(gsl_matrix *A){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            std::cout<<" "<<gsl_matrix_get(A,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_vector(gsl_vector *v){
    std::cout<<"[ ";
    for(int i=0; i<N; i++){
        std::cout<<gsl_vector_get(v,i)<<" ";
    }
    std::cout<<"]"<<std::endl;
}