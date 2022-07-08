#include "local_minimization.h"

void local_minimization::ytox ( double x[], int n )
//  ---------------------------------------------------------------------
//
//  Description:
//    Subroutine YTOX cancels the effect of subroutine XTOY, performing
//    the reverse transformations:
//      For a parameter with a left margin (XL) only:
//        XNEW = XL + XOLD**2
//      For a parameter with a right margin (XR) only:
//        XNEW = XR - XOLD**2
//      For a parameter with both margins:
//        XNEW = (XR-XL)*EXP(-X(I)**2) + XL
//    See also the discussion in subroutine XTOY.
//
//  Input arguments:
//    N          Dimensionality of the objective function.
//
//  Input / Output arguments:
//    X          On input an array of transformed parameters. On output
//               an array of normal Merlin parameters.
//
//  Notes:
//    Minimization routines using the XTOY/YTOX transformations should
//    change the contents of array POINT in common block PARAMS, only
//    upon termination.
//
//  ---------------------------------------------------------------------
{
	int i, ma;
	double xl, xr;

	for (i=0; i<n; i++) {
		ma = detail.marg[i];
		xl = detail.xll[i];
		xr = detail.xrl[i];
		if (detail.ixat[i] == 0) {
			ma = 2;
			xl = detail.point[i];
			xr = detail.point[i];
		}
		if (ma == -1)
			x[i] = xl + pow(x[i],2);
		else if (ma == 1)
			x[i] = xr - pow(x[i],2);
		else if (ma == 2)
			x[i] = (xr-xl)*axp(-pow(x[i],2)) + xl;
	}
}
