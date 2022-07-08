#include "local_minimization.h"

int local_minimization::ifconv ( double val, double oldval, double eps )
//  ---------------------------------------------------------------------
//
//  Description:
//    Checks whether the relative drop of the function value, in
//    two successive iterations of a minimization method is less than
//    EPS.
//
//  Input arguments:
//    VAL        Current value of the objective function.
//    OLDVAL     Value of the objective function at the previous
//               iteration.
//    EPS        Minimum relative drop.
//
//  Returned value:
//    IFCONV = 0 -> Relative drop greater than EPS.
//    IFCONV = 1 -> Relative drop less than or equal to EPS.
//
//  ---------------------------------------------------------------------
{
	double one = 1.0;
	double t1, t2;

	t1 = abs(val-oldval);
	t2 = eps*max(max(abs(val),abs(oldval)),one);
	return t1 <= t2;
}
