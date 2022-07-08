#include "local_minimization.h"

double local_minimization::fcnt ( double x[], int n )
//  ---------------------------------------------------------------------
//
//  Description:
//    This is the objective function that subroutine TRUSMV minimizes.
//
//  Input arguments:
//    X          The minimization variables (transformed).
//    N          Number of variables.
//
//  ---------------------------------------------------------------------
{
      double xtmp[MXV];

//  Back-transform the variables.
	vassign(xtmp,x,n);
	ytox(xtmp,n);
	return acsq(xtmp,n);
}
