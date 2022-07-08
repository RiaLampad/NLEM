#include "local_minimization.h"

void local_minimization::trusmv ( double *h, double point[], double grad[], int nod, 
                                  double *val, int nomaxc, double comacc, double gtol, 
                                  int itmax, double xtol, double rho, double ftol, 
                                  int *info, int *iter, int *nocex, int *nogex, int iprint, 
                                  int ixat[], int ihessu, double xs[], double gn[],
                                  double y[], double gs[], double w[], int iradu, double *delta,
                                  int marg[] )
//  ---------------------------------------------------------------------
//
//  Description:
//    Implements the BFGS method, taking the trust region approach,
//    in place of the line search. Uses Choleski factorization
//    of the Hessian. The Choleski factor L[t] is stored in the upper
//    triangular part of matrix H. Fixed variables are removed, thus
//    reducing the dimensionality of the problem. Margins are treated
//    by appropriate transformations.
//    For more information see:
//      Dennis J.E. & Schnabel R.B.
//      Numerical methods for unconstrained optimization and
//      nonlinear equations
//      Prentice Hall, Englewood Cliffs, New Jersey, 1983.
//
//  Input arguments:
//    GRAD       Work array of length NOD.
//    NOD        Dimensionality of the objective function.
//    NOMAXC     Maximum allowed number of calls to the objective
//               function.
//    COMACC     Machine accuracy.
//    GTOL       G-convergence criterion.
//    ITMAX      Maximum allowed number of BFGS iterations.
//               Set ITMAX = -1 for an unlimited number.
//    XTOL       X-convergence criterion.
//    RHO        Parameter that defines the rho-line.
//    FTOL       F-convergence criterion.
//    IPRINT     Printout level:
//               IPRINT = 0 -> Nothing is printed.
//               IPRINT = 1 -> Simple printout: Prints ITER, VAL,
//                             NOC, MAXNOC
//               IPRINT = 2 -> Full printout. Prints in addition the
//                             minimization parameters, contained in
//                             array POINT.
//    IXAT       Fix statuses.
//    IGRADU     IGRADU = 0 -> Ask subroutine DERMER to recalculate
//                             the gradient at the current point. If
//                             a valid cached gradient exists, DERMER
//                             will return it.
//               IGRADU = 1 -> Ask subroutine DERMER to return its
//                             cached gradient. No recalculation
//                             takes place, unless of course the
//                             gradient cache is empty.
//    XS         Work array of length NOD.
//    GN         Work array of length NOD.
//    Y          Work array of length NOD.
//    GS         Work array of length NOD.
//    W          Work array of length NOD.
//    IRADU      IRADU = 0 -> Initialize the trust region radius.
//               IRADU = 1 -> Use the radius from the previous run.
//    MARG       Margin statuses.
//
//  Output arguments:
//    INFO       Reason for termination:
//               INFO = 1 -> Target value has been reached.
//               INFO = 2 -> The gradient criterion is satisfied.
//               INFO = 3 -> All function calls have been exhausted.
//               INFO = 4 -> X-convergence has been achieved.
//               INFO = 5 -> F-convergence has occurred.
//               INFO = 6 -> Division by zero while updating H.
//               INFO = 7 -> Iterations exhausted.
//               INFO = 8 -> All variables are fixed.
//               INFO = 9 -> Further progress not possible.
//    ITER       Number of BFGS iterations that were performed.
//    NOCEX      Number of function evaluations that were performed.
//    NOGEX      Number of gradient evaluations that were performed.
//    DELTA      The trust region radius.
//
//  Input / Output arguments:
//    H          On input the upper triangular part of H may contain an
//               initial approximation to the Choleski L[t] factor of
//               the Hessian matrix. On output contains the upper
//               triangular part contains the Choleski L[t] of the
//               BFGS approximation to the Hessian.
//    POINT      On input, the current point. On output, the new point
//               after the minimization.
//    VAL        Value of the objective function.
//
//  ---------------------------------------------------------------------
{
	int i, j, m, mm, n1, n2, nma, ispd, iret;
	double oldval, yty, dtd, yts, aa, tt;
	double xtmp[MXV];

//  Number of iterations performed so far.
	*iter = 0;

//  Number of function calls so far.
	*nocex = 0;

//  Number of (analytic) gradient calls so far.
	*nogex = 0;

//  Did we reach our target value ?
//      IF (ITCONV(VAL).NE.0) THEN
//         INFO = 1
//         RETURN
//      END IF

//  Count non-fixed variables. M will be used as the dimensionality
//  of the objective function in all subsequent matrix/vector operations.
	m = freeva();
	if (m == 0) {
		*info = 8;
		return;
	}

//  Calculate derivatives at the starting point.
	vassign(xtmp,point,nod);
	ytox(xtmp,nod);
	dermer(xtmp,*val,grad,&n1,&n2);
	gtrans(grad,point);
	*nocex += n1;
	*nogex += n2;

//  Copy derivatives to an array of M elements.
	mm = 0;
	for (i=0; i<nod; i++)
		if (ixat[i] != 0) {
			gs[mm] = grad[i];
			mm ++ ;
		}

//  Check the gradient criterion for the initial point.
//  The numerical differentiation routine DERMER is expected to
//  return GRAD(I) = 0 for all fixed variables. Otherwise the
//  following test will fail, since the whole array GRAD is examined.
	if (igconv(grad,point,nod,*val,gtol) != 0) {
		*info = 2;
		return;
	}

//  Initialize Hessian. (Upper triangular part only.)
	if (ihessu == 0)
		for (j=0; j<m; j++) {
			h[idx(j,j,nod)] = 1;
			for (i=0; i<=j-1; i++)
				h[idx(i,j,nod)] = 0;
		}

//  Check if margins are present.
	nma = 0;
	for (i=0; i<nod; i++)
		if (marg[i] != 0)
			nma++;
	if (nma != 0) {
//  Transform the Hessian to new coordinates.
		ltrans(h,nod,point);
//  Compose the Hessian from its Choleski factors.
		llcomp(h,nod,nod);
//  Decompose the Hessian forcing possitive definiteness.
		hpdd(h,xtmp,nod,nod,comacc,&ispd);
	}

//  Prepare an initial trust region radius.
	if (iradu == 0) {
		*delta = vvmult(gs,gs,nod);
		*delta = sqrt(*delta)/10;
	}
//  ---------------------------------------------------------------------
//  The following loop performs TRUST iterations.
//  ---------------------------------------------------------------------
L1000:
	if (itmax != -1 && (*iter) >= itmax) {
		*info = 7;
		return;
	}
	if (*nocex >= nomaxc) {
		*info = 3;
		return;
	}
	(*iter) ++ ;

//  Save current point and value for the convergence tests.
	oldval = *val;
	vassign(xs,point,nod);

//  Calculate the next point and update the trust region.
	ttdri(delta,m,gs,h,nod,rho,val,xtol,&n1,point,nod,ixat,gn,y,w,&iret);
	*nocex += n1;

//  Any progress ? (IRET = 1, means the step taken was too small.)
	if (iret == 1 && xtol <= comacc) {
		*info = 9;
		return;
	}
	if (*val < oldval)
		mdis(point,*val,*iter,*nocex,nomaxc,iprint,1);

//  Did we reach our target value ?
//	IF (ITCONV(VAL).NE.0) THEN
//		INFO = 1
//		RETURN
//	END IF

//  Check for X-convergence.
	if (ixconv(point,xs,nod,xtol) != 0) {
		*info = 4;
		return;
	}

//  Check for F-convergence.
	if (ifconv(*val,oldval,ftol) != 0) {
		*info = 5;
		return;
	}

//  Calculate derivatives at the new point.
	vassign(xtmp,point,nod);
	ytox(xtmp,nod);
	dermer(xtmp,*val,gn,&n1,&n2);
	gtrans(gn,point);
	*nocex += n1;
	*nogex += n2;

//  Check the gradient criterion.
	if (igconv(gn,point,nod,*val,gtol) != 0) {
		*info = 2;
		return;
	}

//  Calculate the difference in the gradient: Y = g[k+1]-g[k]
//  and in the variables: XS = x[k+1] - x[k]
//  Update derivatives in GRAD.
	mm = 0;
	for (i=0; i<nod; i++) {
		if (ixat[i] != 0) {
			xs[mm] = point[i]-xs[i];
			y[mm] = gn[i]-grad[i];
			gs[mm] = gn[i];
			mm ++ ;
		}
		grad[i] = gn[i];
	}

//  Calculate vector inner products needed to determine whether the
//  Hessian can be updated.
	yty = vvmult(y,y,m);
	dtd = vvmult(xs,xs,m);
	yts = vvmult(y,xs,m);

//  Can we update the Hessian ?
	if (yts > sqrt(comacc*dtd*yty)) {

//  XS = L[t] * s
		ltmult(h,xs,xs,m,nod);

//  TT = s[t] * H * s
		tt = vvmult(xs,xs,m);
		if (tt == 0.0) {
			*info = 6;
			return;
		}

//  GN = H * s
		lmult(h,xs,gn,m,nod);
		aa = sqrt(yts/tt);
		for (i=0; i<m; i++)
			gn[i] = aa*(y[i]-aa*gn[i])/yts;

//  Update Choleski factors.
		llupd(h,gn,xs,m,nod);
	}
	goto L1000;
//  ---------------------------------------------------------------------
//  End of TRUST iterations.
//  ---------------------------------------------------------------------
}
