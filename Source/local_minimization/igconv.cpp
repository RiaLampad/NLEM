#include "local_minimization.h"

int local_minimization::igconv ( double grad[], double point[], int nod, double val, double gtol )
//  ---------------------------------------------------------------------
//
//  Description:
//    Checks whether any component of the (scaled) gradient vector is
//    larger than GTOL. For details see:
//      Dennis J.E. & Schnabel R.B., Numerical methods for unconstrained
//      optimization and nonlinear equations.
//      Prentice-Hall, 1983. pp. 159-161, 278-280.
//
//  Input arguments:
//    GRAD       The current gradient vector.
//    POINT      The current point.
//    NOD        Dimensionality of the objective function.
//    VAL        Current value of the objective function.
//    GTOL       Tolerance for the scaled gradient.
//
//  Returned value:
//    IGCONV = 0 -> At least one scaled gradient component is larger
//                  than GTOL.
//    IGCONV = 1 -> No scaled gradient component is larger than GTOL.
//
//  ---------------------------------------------------------------------
{
	int i;
	double gmax, gc;
	double one = 1.0;

	gmax = 0;
	for (i=0; i<nod; i++) {
		gc = abs(max(abs(point[i]),one)*grad[i]/max(abs(val),one));
		gmax = max(gmax,gc);
	}
	return gmax <= gtol;
}
