#include "local_minimization.h"

void local_minimization::trupd ( double step[], int m, double rho, double val, double *h, 
                                 int ld, double *delta, double xtol, int *iret, double *fnew, 
                                 double gs[], double point[], int nod, int ixat[], double xnew[],
                                 int idel )
//  ---------------------------------------------------------------------
//
//  Description:
//    Updates the trust region radius.
//
//  Input arguments:
//    STEP       The step to be taken.
//    M          Number of non-fixed variables.
//    RHO        The rho-line parameter.
//    VAL        Current value of the objective function.
//    H          Choleski factor L[t] of the Hessian matrix.
//    LD         Leading dimension of matrix H.
//    XTOL       X-convergence criterion. Defines the minimum
//               acceptable step length.
//    GS         The gradient vector for the M non-fixed variables.
//    POINT      The current point.
//    NOD        Dimensionality of the objective function.
//    IXAT       Fix statuses.
//    IDEL       IDEL = 0 -> STEP is the full Newton step.
//               IDEL = 1 -> STEP is the double dogleg step.
//
//  Output arguments:
//    IRET       Reason for return:
//               IRET = 0 -> A new, lower function value has been found,
//                           and the trust region radius has been
//                           updated.
//               IRET = 1 -> The step to be taken is too small.
//               IRET = 2 -> Failed to produce a lower function value.
//                           The trust region radius has been updated
//                           however. The calling routine must now
//                           calculate a new STEP and call this routine
//                           once more.
//    FNEW       The new function value. (Not necessarily lower than
//               VAL.)
//    XNEW       The new point.
//
//  Input / Output arguments:
//    DELTA      On input the current trust region radius. On output, the
//               new, updated radius.
//
//  ---------------------------------------------------------------------
{
	int mm, i;
	double stplen, g0, xmax, dd, shs, df, dfpred;
	double one = 1.0;

	stplen = vvmult(step,step,m);
	stplen = sqrt(stplen);
	g0 = vvmult(gs,step,m);

//  XNEW is the new point.
	mm = 0;
	for (i=0; i<nod; i++)
		if (ixat[i] == 0)
			xnew[i] = point[i];
		else {
			xnew[i] = point[i]+step[mm];
			mm ++ ;
		}

// This is normally neccesary only if FNEW is above the rho-line.
// (A backtrack is in need).
	xmax = 0;
	mm = 0;
	for (i=0; i<nod; i++)
		if (ixat[i] != 0) {
			xmax = max(xmax,abs(step[mm])/max(abs(xnew[i]),one));
			mm ++ ;
		}
	if (xmax <= xtol) {
		*fnew = val;
		*iret = 1;
		return;
	}

//  Evaluation of the objective function occurs only in the following
//  statement.
	*fnew = fcnt(xnew,nod);
	df = *fnew-val;
	if (df >= rho*g0) {
//
//  Interpolate for the new DELTA.
		dd = -g0*stplen/(2*(df-g0));
		if (dd < 0.1*(*delta)) {
			*delta = 0.1*(*delta);
		} else if (dd > 0.5*(*delta)) {
			*delta = 0.5*(*delta);
		} else {
			*delta = dd;
		}
		*iret = 2;
	}
	else {

//  Update the trust region radius according to the predicted and
//  actual reductions in the function value.
//  SHS = s[t] * H * s
		shs = llmult(h,step,m,ld);
		dfpred = g0+shs/2;
		if (df > 0.25*dfpred)
			*delta = (*delta)/2;
		else if (df < 0.75*dfpred)
			if (idel != 0)
				*delta = 2*(*delta);
		*iret = 0;
	}
}
