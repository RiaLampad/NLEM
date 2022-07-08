#include "local_minimization.h"

double local_minimization::llmult ( double a[], double x[], int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Performs the multiplication B = X[t] * L * L[t] * X
//    where X is a vector and L is a lower triangular Choleski factor.
//    L is input as its traspose L[t].
//
//  Input arguments:
//    A          Upper triangular (Choleski factor L[t]).
//    X          Vector.
//    N          Dimensionality of the problem.
//    LD         Leading (first) dimension of matrix A.
//
//  Output arguments:
//    B          Result of the multiplication (scalar).
//
//  ---------------------------------------------------------------------
{
	int i, j;
	double b, sum;

	b = 0;
	for (i=0; i<n; i++) {
		sum = 0;
		for (j=i; j<n; j++)
			sum += a[idx(i,j,ld)]*x[j];
		b += sum*sum;
	}
	return b;
}
