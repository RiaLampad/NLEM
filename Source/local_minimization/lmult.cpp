#include "local_minimization.h"

void local_minimization::lmult ( double *a, double x[], double b[], int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Performs the multiplication B = L * X.
//    where L is a lower triangular Choleski factor.
//    L is input as its traspose L[t].
//
//  Input arguments:
//    A          Upper triangular (Choleski factor L[t]).
//    X          Vector.
//    N          Dimensionality of the problem.
//    LD         Leading (first) dimension of matrix A.
//
//  Output arguments:
//    B          The resulting left hand vector.
//
//  Notes:
//    X and B normally are two different vectors. If however the calling
//    routine uses the same vector for X and B, the algorithm will work
//    and X will be overwritten.
//
//  ---------------------------------------------------------------------
{
	int i, j;
	double sum;

	for (i=n-1; i>=0; i--) {
		sum = 0;
		for (j=0; j<=i; j++)
			sum += a[idx(j,i,ld)]*x[j];
		b[i] = sum;
	}
}
