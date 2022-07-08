#include "local_minimization.h"

void local_minimization::ltmult ( double *a, double x[], double b[], int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Performs the multiplication B = L[t] * X
//    where L[t] is an upper triangular Cloleski factor.
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

	for (i=0; i<n; i++) {
		sum = 0;
		for (j=i; j<n; j++)
			sum += a[idx(i,j,ld)]*x[j];
		b[i] = sum;
	}
}
