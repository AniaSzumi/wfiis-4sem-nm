# include <iostream>
// # include "/usr/local/include/gsl/gsl_math.h"
# include "/usr/local/include/gsl/gsl_linalg.h"
# define N 4

void print_matrix(gsl_matrix *A);
void save_matrix(gsl_matrix *A, FILE *fp);
float max_value(gsl_matrix *A);

int main(){
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_permutation *p = gsl_permutation_calloc(N);
    int signum;
    double delta = 2.0;
    std::cout.precision(2);
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            gsl_matrix_set(A, i, j, (1.0/(i+j+delta)));
        }
    }
    print_matrix(A);
    FILE *fp;
    fp = fopen("out.txt", "w");
    if (!fp) return -1;
    fprintf(fp, "Macierz: \n");
    save_matrix(A, fp);

    gsl_linalg_LU_decomp(A, p, &signum);
    fprintf(fp,"Macierz LU: \n");
    save_matrix(A, fp);

    float wA = static_cast<float>(signum);
    fprintf(fp,"Diagonala macierzy: \n");
    for(int i=0; i<N; i++){
        wA *= gsl_matrix_get(A,i,i);
        fprintf(fp," %g \n",gsl_matrix_get(A,i,i));
    }
    fprintf(fp,"Wyznacznik macierzy: %g \n\n", wA);

    gsl_vector *b1 = gsl_vector_calloc(N);
    gsl_vector *b2 = gsl_vector_calloc(N);
    gsl_vector *b3 = gsl_vector_calloc(N);
    gsl_vector *b4 = gsl_vector_calloc(N);

    gsl_vector *x1 = gsl_vector_calloc(N);
    gsl_vector *x2 = gsl_vector_calloc(N);
    gsl_vector *x3 = gsl_vector_calloc(N);
    gsl_vector *x4 = gsl_vector_calloc(N);
    for(int i=0; i<N; i++){
        gsl_vector_set(b1,i,0.0);
        gsl_vector_set(b2,i,0.0);
        gsl_vector_set(b3,i,0.0);
        gsl_vector_set(b4,i,0.0);
    }
    gsl_vector_set(b1, 0, 1.0);
    gsl_vector_set(b2, 1, 1.0);
    gsl_vector_set(b3, 2, 1.0);
    gsl_vector_set(b4, 3, 1.0);

    gsl_linalg_LU_solve(A, p, b1, x1);
    gsl_linalg_LU_solve(A, p, b2, x2);
    gsl_linalg_LU_solve(A, p, b3, x3);
    gsl_linalg_LU_solve(A, p, b4, x4);

    gsl_matrix *A_inv = gsl_matrix_calloc(N,N);
    for(int i=0; i<N; i++){
        gsl_matrix_set(A_inv, i, 0, gsl_vector_get(x1, i));
        gsl_matrix_set(A_inv, i, 1, gsl_vector_get(x2, i));
        gsl_matrix_set(A_inv, i, 2, gsl_vector_get(x3, i));
        gsl_matrix_set(A_inv, i, 3, gsl_vector_get(x4, i));
    }

    fprintf(fp, "Macierz odwrotna: \n");
    save_matrix(A_inv, fp);

    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            gsl_matrix_set(A, i, j, (1.0/(i+j+delta)));
        }
    }

    float sum;
    gsl_matrix *C = gsl_matrix_calloc(N,N);
    for (int i=0; i<N ; i++){
        for (int j=0; j<N; j++){
            sum = 0.0; 
            for (int k=0; k<N ; k++){
                sum += gsl_matrix_get(A,i,k) * gsl_matrix_get(A_inv,k,j);
            }
            gsl_matrix_set(C,i,j,sum);
        }
    }
    print_matrix(C);
    fprintf(fp, "Wynik mnozenia: \n");
    save_matrix(C, fp);

    float cond = 0.0;
    cond = max_value(A) * max_value(A_inv);
	fprintf(fp, "Wspolczynnik uwarunkowania : %g\n", cond);

    fclose(fp);
    gsl_matrix_free(A);
    gsl_matrix_free(A_inv);
    gsl_matrix_free(C);
    gsl_vector_free(b1);
    gsl_vector_free(b2);
    gsl_vector_free(b3);
    gsl_vector_free(b4);
    gsl_vector_free(x1);
    gsl_vector_free(x2);
    gsl_vector_free(x3);
    gsl_vector_free(x4);
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

void save_matrix(gsl_matrix *A, FILE *fp){
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            fprintf(fp, " %g ", gsl_matrix_get(A,i,j));
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

}

float max_value(gsl_matrix *A){
    float max = abs(gsl_matrix_get(A,0,0));
    for (int i = 0 ; i < N; i++){
        for (int j = 0; j < N; j++){
            if (abs(gsl_matrix_get(A,i,j)) > max){
                max = gsl_matrix_get(A,i,j);
            }
        }
    }
    return max;
}