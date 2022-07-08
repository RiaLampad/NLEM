#include "local_minimization.h"

double local_minimization::acsq1 ( double xpoint[], int ii, double x )

//  ---------------------------------------------------------------------
//
//  Description:
//    This is the function that the numerical differentiation routines
//    use. It simply sets the II-th component of XPOINT equal to X, and
//    calls the objective function ACSQ.
//
//  Input arguments:
//    XPOINT     The point at which differentiation occurs.
//    II         Minimization variable whose value must be changed.
//    X          New value for variable II.
//
//  ---------------------------------------------------------------------
{
	int nod;
	double xsave, v;

	nod = detail.nod;
	xsave = xpoint[ii];
	xpoint[ii] = x;
	v = acsq(xpoint,nod);
	xpoint[ii] = xsave;
	return v;
}
