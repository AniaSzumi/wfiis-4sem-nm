#include <iostream>
#include <gsl/gsl_linalg.h>
#include <cmath>
#include <vector>
#define x0 2
#define sigma 4

double delta(double x, double alfa){
    double U = rand()/(RAND_MAX+1.0);
    return alfa*(U-0.5);
}

double g(double x, double alfa){
    double a0 = -pow(x0,2)/(2*pow(sigma,2)), a1 = x0/pow(sigma,2), a2 = -1/(2*pow(sigma,2));
    return exp(a0+a1*x+a2*x*x)*(1+delta(x,alfa));
}

double newG(double x, gsl_vector *b){
    return exp(gsl_vector_get(b,0)+gsl_vector_get(b,1)*x+gsl_vector_get(b,2)*pow(x,2)+gsl_vector_get(b,3)*pow(x,3));
}

void solve(int N, double alfa, FILE *pkt, FILE *Gf, FILE *ab){
    double x_min = -3*sigma+x0, x_max = 3*sigma+x0;
    double dx = (x_max-x_min)/static_cast<double>(N-1);

    gsl_vector *xw = gsl_vector_calloc(N);
    gsl_vector *gw = gsl_vector_calloc(N);
    gsl_vector *fw = gsl_vector_calloc(N);
    double x = x_min;
    
    for(int i=0; i<N; i++){
        gsl_vector_set(xw,i,x);
        gsl_vector_set(gw,i,g(x,alfa));
        gsl_vector_set(fw,i,log(gsl_vector_get(gw,i)));
        fprintf(pkt, "%g %g\n", gsl_vector_get(xw,i), gsl_vector_get(gw,i));
        x += dx;
    }
    fprintf(pkt, "\n\n");
    int m = 4;
    gsl_vector *r = gsl_vector_calloc(m);
    gsl_matrix *G = gsl_matrix_calloc(m,m);
    double rk, gik;

    for(int k=0; k<m; ++k){
        rk = 0.0;
        for(int j=0; j<N; ++j){
            rk += gsl_vector_get(fw,j)*pow(gsl_vector_get(xw,j),k);
        }
        gsl_vector_set(r,k,rk);
        for(int i=0; i<m; ++i){
            gik = 0.0;
            for(int j=0; j<N; ++j){
                gik += pow(gsl_vector_get(xw,j),i+k);
            }
            gsl_matrix_set(G,i,k,gik);
        }
    }

    gsl_linalg_HH_svx(G, r);
    
    for(double x=x_min; x<=x_max; x+=0.01){
        fprintf(Gf, "%g %g\n", x, newG(x,r));
    }
    fprintf(Gf, "\n\n");

    double a0 = -pow(x0,2)/(2*pow(sigma,2)), a1 = x0/pow(sigma,2), a2 = -1/(2*pow(sigma,2));
    fprintf(ab, "N = %d:\n", N);
    fprintf(ab, "Analityczne    Numeryczne\n");
    fprintf(ab, "%11f    %10f\n",a0,gsl_vector_get(r,0));
    fprintf(ab, "%11f    %10f\n",a1,gsl_vector_get(r,1));
    fprintf(ab, "%11f    %10f\n",a2,gsl_vector_get(r,2));
    fprintf(ab, "               %10f\n",gsl_vector_get(r,3));

    gsl_vector_free(xw);
    gsl_vector_free(gw);
    gsl_vector_free(fw);
    gsl_matrix_free(G);
    gsl_vector_free(r);

}

int main(){
    FILE *pkt, *G, *ab;
    pkt = fopen("pkt.txt", "w");
    G = fopen("G.txt", "w");
    ab = fopen("ab.txt", "w");

    solve(11,0.0,pkt,G,ab);
    solve(11,0.5,pkt,G,ab);
    solve(101,0.5,pkt,G,ab);

    fclose(pkt);
    fclose(G);
    fclose(ab);


    return 0;
}