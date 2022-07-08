#include "local_minimization.h"

double local_minimization::vvmult ( double x[], double y[], int n )
//  ---------------------------------------------------------------------
//
//  Description:
//    Computes the vector inner product X * Y.
//
//  Input arguments:
//    X          First vecror.
//    Y          Second vector.
//    N          Dimensionality of the vectors.
//
//  Output arguments:
//    XY         The inner product X*Y.
//
//  ---------------------------------------------------------------------
{
	int i;
	double xy;

	xy = 0;
	for (i=0; i<n; i++)
		xy += x[i]*y[i];
	return xy;
}
