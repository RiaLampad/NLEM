#include "local_minimization.h"

void local_minimization::dogleg ( double delta, double gs[], double step[], int m, 
                               double sn[], double snlen, double alc, double sclen, 
                               double an, double gsn, double gg )
//  ---------------------------------------------------------------------
//
// Description:
//    Calculates the double dogleg step.
//
//  Input arguments:
//    DELTA      Number of variables of the objective function.
//    GS         The gradient vector for the M non-fixed variables.
//    M          Number of non-fixed variables.
//    SN         The Newton step.
//    SNLEN      The Newton step length.
//    ALC        The Cauchy step is ALC * GS.
//    SCLEN      The Cauchy step length.
//    AN         Distance along the Newton direction that defines
//               the last part of the double dogleg.
//    GSN        The vector product GS * SN.
//    GG         The vector product GS * GS.
//
//  Output arguments:
//    STEP       The double dogleg step.
//
//  ---------------------------------------------------------------------
{
	int i;
	double usc, uu, c, c1, c2;

	if (an*snlen <= delta) {
//  Take only a partial Newton step (point SN is inside the radius).
		c = -delta/snlen;
		for (i=0; i<m; i++)
			step[i] = c*sn[i];
	}
	else if (sclen >= delta) {
//  Take only a partial steepest descent step (CP is out of the radius).
		c = -delta/sqrt(gg);
		for (i=0; i<m; i++)
			step[i] = c*gs[i];
	}
	else {
//  Take the double dogleg step.
		usc = -alc*(an*gsn+alc*gg);
		uu = 2*an*alc*gsn + pow(an*snlen,2) + alc*alc*gg;
		c = (-usc+sqrt(usc*usc-uu*(sclen*sclen-delta*delta)))/uu;
		c1 = alc*(1.0-c);
		c2 = c*an;
		for (i=0; i<m; i++)
			step[i] = c1*gs[i]-c2*sn[i];
	}
}
