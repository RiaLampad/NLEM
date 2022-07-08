#include "local_minimization.h"

double local_minimization::acsq ( double x[], int n )
//  ---------------------------------------------------------------------
//
//  Description:
//    All Merlin minimization methods that operate on general functions
//    (not on sum of squares) call this routine. It contains a counter
//    and updates some common block parameters. Its value is set equal
//    to the user supplied function FUNMIN or subroutine SUBSUB,
//    according to the functional form (GENERAL or SOS).
//
//  Input arguments:
//    X          The minimization variables.
//    N          Dimensionality of the objective function.
//
//  ---------------------------------------------------------------------
{
//  Increment counters.
	itnocs++;
	return (*detail.funmin)(x,n);
}
