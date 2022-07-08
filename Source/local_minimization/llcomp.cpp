#include "local_minimization.h"

void local_minimization::llcomp ( double *a, int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Performs the multiplication A = L * L[t]
//    where L is a lower triangular Choleski factor.
//    L is input as the traspose L[t].
//
//  Input arguments:
//    N          Dimensionality of the problem.
//    LD         Leading (first) dimension of matrix A.
//
//  Input / Output arguments:
//    A          On input upper triangular (Choleski factor L[t]).
//               On output the product L * L[t]
//
//  ---------------------------------------------------------------------
{
	int i, j, k;
	double s;

	for (i=n-1; i>=0; i--)
		for (j=n-1; j>=i; j--) {
			s = 0;
			for (k=0; k<=i; k++)
				s += a[idx(k,i,ld)] * a[idx(k,j,ld)];
			a[idx(i,j,ld)] = s;
		}

//  Fill in lower part.
	for (i=1; i<n; i++)
		for (j=0; j<i; j++)
			a[idx(i,j,ld)] = a[idx(j,i,ld)];
}
