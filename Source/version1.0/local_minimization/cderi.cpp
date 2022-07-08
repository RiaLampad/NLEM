#include "local_minimization.h"

void local_minimization::cderi ( double point[], int ii, double f0, double stpcon, double *der, int *noc )
//  ---------------------------------------------------------------------
//
//  Description:
//    Subroutine CDERI computes an approximation (DER) to the first
//    partial derivative of the objective function, with respect to
//    variable II. Normally it uses the central difference formula
//    DER = (F(X+H)-F(X-H))/(2*H). If any one of X+H, X-H extend
//    beyond the margins, the 3 point forward or backward difference
//    formulae DER = (4*F(X+H)-3*F(X)-F(X+2*H))/(2*H) and
//    DER = -(4*F(X-H)-3*F(X)-F(X-2*H))/(2*H) are used.
//    If these fail too (X+2*H, X-2*H are outside the margins) DER = 0
//    is returned (the problem is probably very badly scaled).
//
//    For a discussion of numerical differentiation techniques see:
//
//    1) Dennis J.E. and Schnabel R.B.
//       Numerical Methods for Unconstrained Optimization and Nonlinear
//       Equations
//       Prentice-Hall, 1983. pp. 103-105, 322-323.
//    2) Press W.H., Teukolsky S.A., Vetterling W.T., Flannery B.P
//       Numerical recipes, second edition
//       Cambridge University Press, 1992. pp. 180-184.
//
//  Input arguments:
//    POINT      The current point.
//    II         Differentiation is performed with respect to variable
//               II.
//    F0         Current value of the objective function.
//    STPCON     Constant that determines the step size. The cited
//               references discuss the selection of the step size.
//
//  Output arguments:
//    DER        The computed derivative.
//    NOC        Number of calls to the objective function.
//
//  ---------------------------------------------------------------------
{
	int ma;
	double tmp, h;
	double xp, xm, xpp, xmm;
	double fp, fm, fpp, fmm;
	double one = 1.0;

	tmp = point[ii];
	h = stpcon*max(abs(tmp),one);
	ma = detail.marg[ii];
	*noc = 0;

	xp = tmp+h;
	h = xp-tmp;
	xm = tmp-h;

//  Make sure the new "X" lies inside the margins.
	if (ma == -1 || ma == 2)
		if (xm < detail.xll[ii]) {
			xpp = xp+h;
			if (ma == 1 || ma == 2)
				if (xpp > detail.xrl[ii]) {
					*der = 0;
					return;
				}
//  Use the 3 point forward difference formula.
			fpp = acsq1(point,ii,xpp);
			fp = acsq1(point,ii,xp);
			*der = (4*fp-3*f0-fpp)/(2*h);
			*noc = 2;
			return;
		}
	if (ma == 1 || ma == 2)
		if (xp > detail.xrl[ii]) {
			xmm = xm-h;
			if (ma == -1 ||  ma == 2)
				if (xmm < detail.xll[ii]) {
					*der = 0;
					return;
				}
//  Use the 3 point backward difference formula.
			fmm = acsq1(point,ii,xmm);
			fm = acsq1(point,ii,xm);
			*der = -(4*fm-3*f0-fmm)/(2*h);
			*noc = 2;
			return;
		}

	fp = acsq1(point,ii,xp);
	fm = acsq1(point,ii,xm);
	*der = (fp-fm)/(2*h);
	*noc = 2;
}
