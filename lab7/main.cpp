#include <iostream>
// #include "/usr/local/include/gsl/gsl_matrix.h"
// #include "/usr/local/include/gsl/gsl_vector.h"
#include <cmath>
#include <vector>
#define abs(X) ((X)>0? (X):-(X))

double f(double x){
    return exp(-x*x);
}

double Lagrange(std::vector<double>& x, std::vector<double>& y, int n, double xi);

int main(){
    int n=5;
    double x_min = -5.0, x_max = 5.0;
    std::vector<double> xr, xc, yr, yc;
    FILE *fp1, *fp2;
    fp1 = fopen("zad_1.txt", "w");
    fp2 = fopen("zad_2.txt", "w");
    for(int n=5; n<21; n+=5){
        double h = (x_max-x_min)/static_cast<double>(n);
        xr.reserve(n+1);
        xc.reserve(n+1);
        yr.reserve(n+1);
        yc.reserve(n+1);
        for(int m=0; m<=n; m++){
            xr[m] = x_min + m*h;
            yr[m] = f(xr[m]);
            xc[m] = ((x_max-x_min)*cos(M_PI*(2*m+1)/(2*n+2))+(x_max+x_min))/2;
            yc[m] = f(xc[m]);
            std::cout << xc[m] << " " << yc[m] << std::endl;
        }
        xr[n] = x_max;
        for(double xi=x_min; xi<=x_max; xi+=0.01){
            fprintf(fp1, "%g %g\n", xi, Lagrange(xr, yr, n, xi));
            fprintf(fp2, "%g %g\n", xi, Lagrange(xc, yc, n, xi));
        }
        
        std::cout << std::endl;
        fprintf(fp1, "\n\n");
        fprintf(fp2, "\n\n");
        xr.clear();
        xc.clear();
        yr.clear();
        yc.clear();
    }
    fclose(fp1);
    fclose(fp2);
    return 0;
}

double Lagrange(std::vector<double>& x, std::vector<double>& y, int n, double xi){
    double w=0.0;
    for(int j=0; j<=n; j++){
        double p=1.0;
        for(int k=0; k<=n; k++){
            if(k!=j){
                p *= ((xi-x[k])/(x[j]-x[k]));
            }
        }
        w += y[j] * p;
    }
    return w;
}