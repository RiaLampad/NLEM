#include "local_minimization.h"

double local_minimization::axp ( double z )
//  ---------------------------------------------------------------------
//
//  Description:
//    Calculates EXP(Z) avoiding underflows. If EXP(Z) would cause an
//    underflow, AXP is set to 0.
//
//  Input arguments:
//    Z          Argument for EXP(Z)
//
//  ---------------------------------------------------------------------
{
	if (z <= machine_constants.dwalog)
		return 0.0;
	else
		return exp(z);
}
