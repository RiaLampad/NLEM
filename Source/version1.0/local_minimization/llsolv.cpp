#include "local_minimization.h"

void local_minimization::llsolv ( double *a, double x[], double b[], int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Solves the set of linear equations A * X = B by direct forward and
//    backward substitution. A is input as its upper triangular Choleski
//    factor L[t].
//
//  Input arguments:
//    A          Upper triangular (Choleski factor L[t]).
//    B          The right hand vector.
//    N          Dimensionality of the problem.
//    LD         Leading (first) dimension of matrix A.
//
//  Output arguments:
//    X          The solution vector.
//
//  Notes:
//    X and B normally are two different vectors. If however the calling
//    routine uses the same vector for X and B, the algorithm will work
//    and B will be overwritten.
//
//  ---------------------------------------------------------------------
{
	int i, j;
	double sum;

	for (i=0; i<n; i++) {
		sum = 0;
		for (j=0; j<=i-1; j++)
			sum += a[idx(j,i,ld)]*x[j];
		x[i] = (b[i]-sum)/a[idx(i,i,ld)];
}
	for (i=n-1; i>=0; i--) {
		sum = 0;
		for (j=i+1; j<n; j++)
			sum += a[idx(i,j,ld)]*x[j];
		x[i] = (x[i]-sum)/a[idx(i,i,ld)];
	}
}
