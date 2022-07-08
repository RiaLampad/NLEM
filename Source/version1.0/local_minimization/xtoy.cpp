#include "local_minimization.h"

void local_minimization::xtoy ( double x[], int n )
//  ---------------------------------------------------------------------
//
//  Description:
//    Subroutine XTOY eliminates margins and fix statuses by transforming
//    the minimization parameters. This is necessary for routines that
//    cannot handle margins and fixed variables directly.
//    The tranformations used are:
//      For a parameter with a left margin (XL) only:
//        XNEW = SQRT(XOLD-XL)
//      For a parameter with a right margin (XR) only:
//        XNEW = SQRT(XR-XOLD)
//      For a parameter with both margins:
//        XNEW = SQRT(LOG((XR-XL)/(XOLD-XL)
//    Subroutine YTOX performs the reverse transformation. For a general
//    discussion on tranformations and some of the pitfalls, see:
//      Sisser F.S., Elimination of bounds in optimization problems by
//      transforming variables.
//      Math. Prog. 20 (1981) 110-121.
//
//  Input arguments:
//    N          Dimensionality of the objective function.
//
//  Input / Output arguments:
//    X          On input an array of Merlin parameters. On output
//               the transformed parameters, as described earlier.
//
//  Notes:
//    Minimization routines using the XTOY/YTOX transformations should
//    change the contents of array POINT in common block PARAMS, only
//    upon termination.
//
//  ---------------------------------------------------------------------
{
	int i, ma;
	double xl, xr, bma, dxl, d;

	for (i=0; i<n; i++) {
		ma = detail.marg[i];
		xl = detail.xll[i];
		xr = detail.xrl[i];
		if (detail.ixat[i] == 0) {
			ma = 2;
			xl = x[i];
			xr = x[i];
		}
		if (ma == -1)
			x[i] = sqrt(x[i]-xl);
		else if (ma == 1)
			x[i] = sqrt(xr-x[i]);
		else if (ma == 2) {
			bma = xr-xl;
			if (bma != 0.0) {
				d = 0.0;
				dxl = x[i]-xl;
				if (dxl == 0.0)
					d = machine_constants.comacc;
				x[i] = sqrt(log(bma/(dxl+d)));
			} else
				x[i] = 1.0;
		}
	}
}
