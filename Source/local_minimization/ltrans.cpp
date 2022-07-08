#include "local_minimization.h"

void local_minimization::ltrans ( double *h, int ld, double x[] )
//  ---------------------------------------------------------------------
//
//  Description:
//    This subroutine returns the Choleski factor L[t] (of the Hessian
//    matrix) with respect to a set of transformed variables.
//    This is necessary for minimization methods that handle bounds via
//    transformations.
//    See also subroutines XTOY and YTOX.
//
//  Input arguments:
//    LD         The leading dimension of matrix H.
//    X          The transformed variables.
//
//  Input / Output arguments:
//    H          On input the upper triangular part contains a
//               Choleski factor L[t].
//               On output the upper triangular part contains the
//               Choleski factor L[t] with respect to the transformed
//               variables.
//
//  ---------------------------------------------------------------------
{
	double two = 2;
	int i, mm, ma, nod;
	double con, yy;

	nod = detail.nod;
	mm = 0;
	for (i=0; i<nod; i++) {
		if (detail.ixat[i] != 0) {
			ma = detail.marg[i];
			if (ma == 0)
				con = 1;
			else if (ma == -1)
				con = two*x[i];
			else if (ma == 1)
				con = -two*x[i];
			else {
				yy = x[i];
				con = -two*yy*(detail.xll[i]-detail.xrl[i])*axp(-yy*yy);
			}
			h[idx(mm,mm,ld)] *= con;
			mm ++ ;
		}
	}
}
