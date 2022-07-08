#include "local_minimization.h"

void local_minimization::jrot ( double *a, int i, double aa, double bb, int n, int ld )
//  ---------------------------------------------------------------------
//
//  Description:
//    Performs a Jacobi rotation on rows I and I+1 of the NxN matrix A.
//    The parameters of the rotation are:
//      COS(phi) = AA/SQRT(AA**2+BB**2)  SIN(phi) = BB/SQRT(AA**2+BB**2)
//
//    For more information see:
//
//    Dennis J.E. and Schnabel R.B.
//    Numerical Methods for Unconstrained Optimization and Nonlinear
//    Equations
//    Prentice-Hall, 1983. p. 55-58.
//
//  Input arguments:
//    I          Rotation is performed on rows I, I+1.
//    AA         First parameter of the rotation.
//    BB         Second parameter of the rotation.
//    N          Dimension of matrix A.
//    LD         First (leading) dimension of matrix A.
//
//  Input / Output arguments:
//    A          NxN matrix to be rotated.
//
//  ---------------------------------------------------------------------
{
	double one = 1.0;
	int i1, j, k1, k2;
	double c, s, r, y, w;

	if (aa == 0.0) {
		c = 0;
		s = sign(one,bb);
	} else {
		r = sqrt(aa*aa+bb*bb);
		c = aa/r;
		s = bb/r;
	}

	i1 = i+1;
	for (j=i; j<n; j++) {
		k1 = idx(i,j,ld);
		k2 = idx(i1,j,ld);
		y = a[k1];
		w = a[k2];
		a[k1] = c*y-s*w;
		a[k2] = s*y+c*w;
	}
}
