#include "local_minimization.h"

void local_minimization::ttdri ( double *delta, int m, double gs[], double *h, 
                                 int ld, double rho, double *val, double xtol, int *nocs,
                                 double point[], int nod, int ixat[], double step[], 
                                 double sn[], double xnew[], int *iret )
//  ---------------------------------------------------------------------
//
//  Description:
//    Implements the trust region approach. It calculates the Cauchy
//    and Newton points, takes the appropriate step, and updates the
//    trust region radius accordingly.
//    For more information see:
//      Dennis J.E. & Schnabel R.B.
//      Numerical methods for unconstrained optimization and
//      nonlinear equations
//      Prentice Hall, Englewood Cliffs, New Jersey, 1983. pp 139-147.
//
//  Input arguments:
//    M          Number of non-fixed variables. This is used as the
//               dimensionality of the problem.
//    GS         The gradient vector for the M non-fixed variables.
//    H          Contains the Choleski factor L[t] on the upper
//               triangular part.
//    LD         Leading dimension of matrix H.
//    RHO        The rho-line parameter.
//    XTOL       X-convergence criterion.
//    NOD        Dimensionality of the objective function.
//    IXAT       Fix statuses.
//    STEP       Work array of length NOD.
//    SN         Work array of length NOD.
//    XNEW       Work array of length NOD.
//
//  Output arguments:
//    NOCS       Number of function calls.
//    IRET       Reason for return:
//               IRET = 0 -> A new point has been successfully found.
//               IRET = 1 -> The step to be taken is too small. (It will
//                           not result in a new point.)
//
//  Input / Output arguments:
//    DELTA      On input the trust region radius. On output, the new,
//               updated radius.
//    VAL        Value of the objective function.
//    POINT      On input the current point. On output, the new, updated
//               point.
//
//  ---------------------------------------------------------------------
{
	int i, first, idel;
	double snlen, gg, gsn, alc, ghg, sclen, an, fnew;

//  Number of function calls.
	*nocs = 0;

//  Solve H * s = -g for the Newton step SN. The minus sign will be
//  added later.
	llsolv(h,sn,gs,m,ld);

//  SNLEN is the Newton step length.
	snlen = vvmult(sn,sn,m);
	snlen = sqrt(snlen);

      first = 1;
L1000:
	if (snlen <= *delta) {
		for (i=0; i<m; i++)
			step[i] = -sn[i];
		idel = 0;
	}
	else {
		if (first) {
			first = 0;

//  GHG = g[t] * H * g
//  GG = g[t] * g
			ghg = llmult(h,gs,m,ld);
			gg = vvmult(gs,gs,m);

//  Cauchy step: SC = -al * g,  al = (g*g)/(g*H*g)
			alc = gg/ghg;

//  SCLEN is the Cauchy step length.
			sclen = alc*sqrt(gg);
			alc = -alc;
			gsn = vvmult(gs,sn,m);
			an = 0.2 + 0.8*gg*gg/(ghg*gsn);
		}
		dogleg(*delta,gs,step,m,sn,snlen,alc,sclen,an,gsn,gg);
		idel = 1;
	}
	trupd(step,m,rho,*val,h,ld,delta,xtol,iret,&fnew,gs,point,
            nod,ixat,xnew,idel);
	*nocs += 1;
	if (*iret == 2) goto L1000;

	if (*iret == 0 || fnew < *val) {
		for (i=0; i<nod; i++)
			point[i] = xnew[i];
		*val = fnew;
	} 
}
