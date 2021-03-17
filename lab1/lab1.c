#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Dyrektywy zakladajace, ze te trzy pliki sa skopiowane do aktualnego katalogu. */
#include "../nr/nrutil.h"
#include "../nr/nrutil.c" // To mozna usunac, jesli plik jest dodany w poleceniu kompilacji.
#include "../nr/gaussj.c" // To tez mozna usunac, jesli plik jest dodany w poleceniu kompilacji.

/* Dyrektywy dla Taurusa (nie wymagaja kopiowania plikow, ale Taurus musi dzialac...) */
// #include "/opt/NR/numerical_recipes.c/nrutil.h"
// #include "/opt/NR/numerical_recipes.c/nrutil.c"
// #include "/opt/NR/numerical_recipes.c/gaussj.c"

#define N 3 // rozmiar macierzy M: NxN

int main(void)
{
	float **M, **b;
	//	Alokacja macierzy
	M = matrix(1, N, 1, N);
	b = matrix(1, N, 1, 1);
	int i=1, j=1;
	float h = 0.1;

	// 	Wypelnienie macierzy M i wektora b
	for (i = 1; i <= N; ++i)
	{
		b[i][1] = 0.0;
		for (j = 1; j <= N; ++j)
			M[i][j] = 0.0;
			if(i == j || i-2 == j){
				M[i][j] = 1.0;
			}
			if(i-1 == j){
				M[i][j] = h*h - 2;
			}
	}

	b[1][1] = 1.0;
	b[2][1] = 0;
	// b[1][1] = -4.0;
	// b[2][1] = 9.0;
	// b[3][1] = -16.0;

	// M[1][1] = 2.0;
	// M[1][2] = 3.0;
	// M[1][3] = 0.0;
	// M[2][1] = -1.0;
	// M[2][2] = 1.0;
	// M[2][3] = 4.0;
	// M[3][1] = 3.0;
	// M[3][2] = 11.0;
	// M[3][3] = 1.0;

	//	Rozwiazanie ukladu rownan Mx=b - wywolanie procedury:
	gaussj(M, N, b, 1);

	//	Wypisanie rozwiazania, ktore procedura gaussj(M, N, b, 1); zapisala w wektorze b.
	for (int i = 1; i <= N; ++i)
		printf("%g\n", b[i][1]);

	//	Zwolnienie pamieci
	free_matrix(M, 1, N, 1, N);
	free_matrix(b, 1, N, 1, 1);

	return 0;
}