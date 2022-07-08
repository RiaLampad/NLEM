#include "local_minimization.h"

void local_minimization::lldec ( double *a, double adiag[], int ld, int n, double oflmax, double *addmax, double comacc )
//  ---------------------------------------------------------------------
//
//  Description:
//    Calculates the perturbed Choleski decomposition of matrix A.
//    It actualy decomposes A+D, where D is a diagonal matrix, with
//    non-negative elements, large enough to allow the decomposition
//    to proceed. For details on the algorithm, see:
//
//      Dennis J.E & Schnabel R.B.
//      Numerical Methods for Unconstrained Optimization and Nonlinear
//      Equations, Prentice-Hall, 1983.
//
//      Gill P.E., Murray W. and Wright M.H.,
//      Practical Optimization, Academic Press, 1989.
//
//  Input arguments:
//    ADIAG      Diagonal elements of A.
//    LD         Leading dimension of matrix A.
//    N          Dimensionality of matrix A.
//    OFLMAX     On output, all elements L(i,j) of the Choleski
//               factor L, are less than, or equal to OFLMAX.
//    COMACC     Machine accuracy.
//
//  Output arguments:
//    ADDMAX     The maximum element, added to the diagonal of A.
//
//  Input / Output arguments:
//    A          On input the lower triangular part of A (excluding the
//               diagonal emements, which are input on ADIAG), contains
//               the matrix to be decomposed.
//               On output, the input part is intact, while the upper
//               triangular part plus the diagonal elements of A, contain
//               the Choleski factor L[t].
//
//  ---------------------------------------------------------------------
{
	int i, j, k;
	double aminl, aminl2, sum, aminjj;

	aminl = sqrt(sqrt(comacc));
	aminl2 = 0;
	if (oflmax == 0) {
		oflmax = 0;
		for (i=0; i<n; i++) 
			oflmax = max(oflmax,abs(adiag[i]));
		aminl2 = sqrt(comacc)*oflmax;
	}
	*addmax = 0;

	for (i=0; i<n; i++)
		for (j=0; j<=i; j++)
			a[idx(j,i,ld)] = 0;

	for (j=0; j<n; j++) {
		sum = 0;
		for (i=0; i<=j-1; i++)
			sum += pow(a[idx(i,j,ld)],2);
		a[idx(j,j,ld)] = adiag[j]-sum;
		aminjj = 0;
		for (i=j+1; i<n; i++) {
			sum = 0;
			for (k=0; k<=j-1; k++)
				sum += a[idx(k,i,ld)]*a[idx(k,j,ld)];
			a[idx(j,i,ld)] = a[idx(i,j,ld)]-sum;
			aminjj = max(aminjj,abs(a[idx(j,i,ld)]));
		}
		aminjj = max(aminjj/oflmax,aminl);
		if (a[idx(j,j,ld)] > aminjj*aminjj)
			a[idx(j,j,ld)] = sqrt(a[idx(j,j,ld)]);
		else {
			if (aminjj < aminl2)
				aminjj = aminl2;
			*addmax = max(*addmax,aminjj*aminjj-a[idx(j,j,ld)]);
			a[idx(j,j,ld)] = aminjj;
		}
		for (i=j+1; i<n; i++)
			a[idx(j,i,ld)] /= a[idx(j,j,ld)];
	}
}
