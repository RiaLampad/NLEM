#include "local_minimization.h"

void local_minimization::gtrans ( double g[], double x[] )
//  ---------------------------------------------------------------------
//
//  Description:
//    This subroutine returns the gradient vector with respect to
//    a set of transformed variables.
//    This is necessary for minimization methods that handle bounds via
//    transformations.
//    See also subroutines XTOY and YTOX.
//
//  Input arguments:
//    X          The transformed variables.
//
//  Input / Output arguments:
//    G          On input the gradient as calculated by
//               SUBROUTINE DERMER.
//               On output the gradient with respect to the transformed
//               variables.
//
//  ---------------------------------------------------------------------
{
	int i, ma, nod;
	double xi, two = 2.0;

	nod = detail.nod;
	for (i=0; i<nod; i++) {
		if (detail.ixat[i] != 0) {
			ma = detail.marg[i];
			if (ma == -1)
				g[i] = two*x[i]*g[i];
			else if (ma == 1)
				g[i] = -two*x[i]*g[i];
			else if (ma == 2) {
				xi = x[i];
				g[i] = -two*xi*(detail.xrl[i]-detail.xll[i])*axp(-xi*xi)*g[i];
			}
		}
	}
}
