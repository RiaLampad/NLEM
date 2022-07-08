#include "local_minimization.h"

void local_minimization::dermer ( double point[], double val, double grad[], int *noc, int *ngr )
//  ---------------------------------------------------------------------
//
//  Description:
//    This is the main Merlin routine that computes approximations
//    to the gradient vector, according to the curent derivative mode.
//    Subroutine DERMER uses an internal cache (GRALOC) to store the
//    computed derivatives. Depending on the setting of IGRADU and the
//    requested operation (MODE), DERMER can return the cached gradient,
//    saving thus unecessary computations.
//
//  Input arguments:
//    POINT      The point where derivatives are requested.
//    VAL        The value of the objective function calculated at POINT.
//    MODE       The requested operation:
//               MODE = 0 -> Normal operation. If the gradient is
//                           available in the cache, it is retured.
//                           Otherwise the derivatives are evaluated
//                           using the current derivative mode, stored
//                           in the gradient cache, and returned to
//                           the calling routine.
//               MODE = 1 -> The gradient vector is calculated. The
//                           gradient cache is not consulted beforehand,
//                           nor is the computed vector stored in the
//                           cache.
//               MODE = 2 -> Returns the cached gradient (does not
//                           evaluate gradient, unless the gradient
//                           cache is empty).
//
//  Output arguments:
//    GRAD       The gradient vector.
//    NOC        Number of function evaluations.
//    NGR        Number of calls to the user supplied subroutine GRANAL.
//
//  ---------------------------------------------------------------------
{
	int n, ii, nocc;
	double feps, con2, der, comacc;
	double one = 1.0, three = 3.0;

	n = detail.nod;

//  Total number of function calls.
	*noc = 0;

//  Total number of GRANAL calls.
	*ngr = 0;

//  FEPS should be initialized to something more appropriate.
//  It is supposed to be the accuracy to which the objective function
//  is calculated.
  	comacc = machine_constants.comacc;
	feps = comacc;
	con2 = pow(feps,(one/three));

//  Loop over all variables.
	for (ii=0; ii<n; ii++)
		if (detail.ixat[ii] == 0)
			grad[ii] = 0;
		else {
			cderi(point,ii,val,con2,&der,&nocc);
			grad[ii] = der;
			*noc = *noc+nocc;
		}
}
