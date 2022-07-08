#include "local_minimization.h"

void local_minimization::llupd ( double *a, double v[], double u[], int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Updates the Choleski factors of a matrix with the rank 2 update:
//    L[new] = L + V * U[t]
//    The matrix is input as the upper triangular Choleski factor L[t].
//
//    For more information see:
//
//    Dennis J.E. and Schnabel R.B.
//    Numerical Methods for Unconstrained Optimization and Nonlinear
//    Equations
//    Prentice-Hall, 1983. p. 55-58.
//
//  Input arguments:
//    V          Vector.
//    U          Vector.
//    N          Dimensionality of the problem.
//    LD         Leading (first) dimension of matrix A.
//
//  Input / Output arguments:
//    A          Upper triangular (Choleski factor L[t]).
//               On output contains the updated factor.
//
//  ---------------------------------------------------------------------
{
	int i, j, k;
	double aa, bb;

//  Zero first lower diagonal.
	for (i=1; i<n; i++)
		a[idx(i,i-1,ld)] = 0.0;

//  Search for non-zero U(I).
	for (k=n-1; k>=0; k--)
		if (u[k] != 0.0) goto L20;
	k = 0;

L20:
	for (i=k-1; i>=0; i--) {
		aa = u[i];
		bb = -u[i+1];
		jrot(a,i,aa,bb,n,ld);
		if (aa == 0.0)
			u[i] = abs(bb);
		else
			u[i] = sqrt(aa*aa+bb*bb);
	}

	for (j=0; j<n; j++)
		a[idx(0,j,ld)] += u[0]*v[j];

	for (i=0; i<=k-1; i++) {
		aa = a[idx(i,i,ld)];
		bb = -a[idx(i+1,i,ld)];
		jrot(a,i,aa,bb,n,ld);
	}
}
