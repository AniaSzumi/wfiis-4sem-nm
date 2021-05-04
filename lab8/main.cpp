#include <iostream>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <vector>

double f1(double x){
    return 1/(1+x*x);
}

double f2(double x){
    return cos(2*x);
}

void print_matrix(gsl_matrix *A, int n);
void print_vector(gsl_vector *v, int n);

void wyznacz_M(gsl_vector *xw ,gsl_vector *yw , gsl_vector *m, int n , double alfa , double beta);

double wyznacz_Sx(gsl_vector *xw ,gsl_vector *yw , gsl_vector *m, int n, double x);

void solve(int n, FILE *f1p, FILE *f2p);

void poch(FILE *fpp);

int main(){
    std::vector<int> size = {5, 8, 21};

    FILE *f1p, *f2p, *fpp;
    f1p = fopen("f1.txt", "w");
    f2p = fopen("f2.txt", "w");
    fpp = fopen("pochodne.txt", "w");

    for (int n : size){
        solve(n, f1p, f2p);
    }

    poch(fpp);

    fclose(f1p);
    fclose(f2p);
    fclose(fpp);
    return 0;
}

double dfdx(double x){
    double dx = 0.01;
    return (f1(x-dx)-2*f1(x)+f1(x+dx))/(dx*dx);
}

void poch(FILE *fpp){
    int n = 10;
    double x_min = -5.0, x_max = 5.0;
    double dx = (x_max-x_min)/static_cast<double>(n-1);

    gsl_vector *xw = gsl_vector_calloc(n);
    gsl_vector *yw1 = gsl_vector_calloc(n);
    gsl_vector *m1 = gsl_vector_calloc(n);

    for (int i=0; i<n; ++i){
        gsl_vector_set(xw,i,x_min+dx*i);
        gsl_vector_set(yw1,i,f1(gsl_vector_get(xw,i)));
    }
    wyznacz_M(xw, yw1, m1, n, 0.0, 0.0);

    for (int i=0; i<n; i++){
        fprintf(fpp, "%g %g %g\n", gsl_vector_get(xw,i), gsl_vector_get(m1,i), dfdx(gsl_vector_get(xw,i)));
    }

    gsl_vector_free(xw);
    gsl_vector_free(yw1);
    gsl_vector_free(m1);
}

void solve(int n, FILE *f1p, FILE *f2p){
    double x_min = -5.0, x_max = 5.0;
    double dx = (x_max-x_min)/static_cast<double>(n-1);

    gsl_vector *xw = gsl_vector_calloc(n);
    gsl_vector *yw1 = gsl_vector_calloc(n);
    gsl_vector *yw2 = gsl_vector_calloc(n);
    gsl_vector *m1 = gsl_vector_calloc(n);
    gsl_vector *m2 = gsl_vector_calloc(n);

    for (int i=0; i<n; ++i){
        gsl_vector_set(xw,i,x_min+dx*i);
        gsl_vector_set(yw1,i,f1(gsl_vector_get(xw,i)));
        gsl_vector_set(yw2,i,f2(gsl_vector_get(xw,i)));
    }
    // print_vector(xw, n);
    wyznacz_M(xw, yw1, m1, n, 0.0, 0.0);
    wyznacz_M(xw, yw2, m2, n, 0.0, 0.0);
    
    for(double x = x_min; x <= x_max; x+=0.01){
        fprintf(f1p, "%g %g\n", x, wyznacz_Sx(xw, yw1, m1, n, x));
        fprintf(f2p, "%g %g\n", x, wyznacz_Sx(xw, yw2, m2, n, x));
    }

    fprintf(f1p, "\n\n");
    fprintf(f2p, "\n\n");

    gsl_vector_free(xw);
    gsl_vector_free(yw1);
    gsl_vector_free(yw2);
    gsl_vector_free(m1);
    gsl_vector_free(m2);
}


void wyznacz_M(gsl_vector *xw ,gsl_vector *yw , gsl_vector *m, int n , double alfa , double beta){
    gsl_matrix *A = gsl_matrix_calloc(n,n);
    gsl_vector *d = gsl_vector_calloc(n);

    gsl_matrix_set(A,0,0,1);
    gsl_matrix_set(A,n-1,n-1,1);
    gsl_vector_set(d,0,alfa);
    gsl_vector_set(d,n-1,beta);

    double lambda_i, di, hi, hi1 = gsl_vector_get(xw,1) - gsl_vector_get(xw,0);

    for(int i=1; i<n-1; ++i){
        hi = hi1;
        hi1 = gsl_vector_get(xw,i+1) - gsl_vector_get(xw,i);
        di = 6*((gsl_vector_get(yw,i+1)-gsl_vector_get(yw,i))/hi1 - (gsl_vector_get(yw,i)-gsl_vector_get(yw,i-1))/hi)/(hi+hi1);
        gsl_vector_set(d,i,di);
        lambda_i = hi1/(hi+hi1);
        // std::cout<<lambda_i<<std::endl;
        gsl_matrix_set(A,i,i,2);
        gsl_matrix_set(A,i,i-1,1-lambda_i);
        gsl_matrix_set(A,i,i+1,lambda_i);
    }
    // print_matrix(A, n);
    // print_vector(d,n);

    gsl_linalg_HH_svx(A, d);
    gsl_vector_memcpy(m,d);
    // print_vector(m,n);
}

double wyznacz_Sx(gsl_vector *xw ,gsl_vector *yw , gsl_vector *m, int n, double x){ 
    // znajdz pierwszy podprzedział ( i −1): xw[i − 1] <= x <= xw[i] 
    // int i = 1;
    // // std::cout << x << " ";
    // while(gsl_vector_get(xw,i-1) > x){
    //     ++i;
    //     std::cout << i << " ";
    // }
    // i--;
    int i = n-1;
    for(int k=1; k<n-1; ++k){
        if(gsl_vector_get(xw,k-1) <= x && gsl_vector_get(xw,k) >= x){
            i = k;
        }
    }

    double hi = gsl_vector_get(xw,1) - gsl_vector_get(xw,0);
    double Ai = (gsl_vector_get(yw,i) - gsl_vector_get(yw,i-1))/hi - hi*(gsl_vector_get(m,i) - gsl_vector_get(m,i-1))/6;
    double Bi = gsl_vector_get(yw,i-1) - gsl_vector_get(m,i-1)*hi*hi/6;
    double Sx = gsl_vector_get(m,i-1)*pow((gsl_vector_get(xw,i)-x),3)/(6*hi) + 
                gsl_vector_get(m,i)*pow(x-(gsl_vector_get(xw,i-1)),3)/(6*hi) +
                Ai*(x - gsl_vector_get(xw,i-1)) + Bi;
    return Sx;
}

void print_matrix(gsl_matrix *A, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            std::cout<<" "<<gsl_matrix_get(A,i,j)<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

void print_vector(gsl_vector *v, int n){
    for(int i=0; i<n; i++){
        std::cout<<gsl_vector_get(v,i);
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

