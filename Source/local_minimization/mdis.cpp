#include "local_minimization.h"

void local_minimization::mdis ( double point[], double val, int iter, int noc, int maxnoc, int iprint, int itra )
//  ---------------------------------------------------------------------
//
//  Description:
//    Issues informative messages during a minimization session.
//    Normally, it is called when a lower value is discovered.
//
//  Input arguments:
//    POINT      The current point.
//    VAL        Value of the objective function at POINT.
//    ITER       Number of iterations that have been performed by the
//               calling routine.
//    NOC        Number of function evaluations that have been performed
//               by the calling routine.
//    MAXNOC     The calling routine is allowed MAXNOC calls at most.
//    IPRINT     Printout level:
//               IPRINT = 0 -> Nothing is printed.
//               IPRINT = 1 -> Simple printout: Prints ITER, VAL,
//                             NOC, MAXNOC
//               IPRINT = 2 -> Full printout. Prints in addition the
//                             minimization parameters, contained in
//                             array POINT.
//    ITRA       Important only when IPRINT=2.
//               ITRA = 0  -> The calling method does not use
//                            transformations.
//               ITRA.NE.0 -> The calling uses tranformations. When
//                            IPRINT=2, we should apply the inverse
//                            transformation before printing the
//                            parameters.
//
//  ---------------------------------------------------------------------
{
	int i, nod;
	double xtmp[MXV];

	if ( iprint == 0 ) return;

	cout << scientific << setprecision(14);
	cout << " Iter: " << iter <<
	        "    Lower value: " << val <<
	        "    Calls: " << noc << " of " << maxnoc << endl;

	if ( iprint == 2 ) {
		nod = detail.nod;
		vassign(xtmp,point,nod);
		ytox(xtmp,nod);
		for (i=0; i<nod; i++)
			cout << setw(4) << i << ")  " << setw(21) << xtmp[i] << endl;
	}
	cout << "----------------------------------------------------------------------" << endl;
}
