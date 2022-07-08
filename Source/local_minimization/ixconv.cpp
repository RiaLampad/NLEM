#include "local_minimization.h"

int local_minimization::ixconv ( double point[], double oldpoi[], int nod, double xtol )
//  ---------------------------------------------------------------------
//
//  Description:
//    Checks whether the scaled distance between two successive iterates
//    of a minimization routine are considered close enough (less than
//    XTOL). For details see:
//      Dennis J.E. & Schnabel R.B., Numerical methods for unconstrained
//      optimization and nonlinear equations.
//      Prentice-Hall, 1983. pp. 159-161, 278-280.
//
//  Input arguments:
//    POINT      The current POINT.
//    OLDPOI     The point at the previous iteration of the calling
//               routine.
//    NOD        Dimensionality of the objective function.
//    XTOL       X-tolerance.
//
//  Returned value:
//    IXCONV = 0 -> At least one minimization parameter has a scaled
//                  distance > XTOL, from the previous iterate.
//    IXCONV = 1 -> Scaled distance less than XTOL for all minimization
//                  parameters.
//
//  ---------------------------------------------------------------------
{
	int i;
	double xmax, xc;
	double one = 1.0;

	xmax = 0;
	for (i=0; i<nod; i++) {
		xc = abs((point[i]-oldpoi[i])/max(abs(oldpoi[i]),one));
		xmax = max(xmax,xc);
	}
	return xmax <= xtol;
}
