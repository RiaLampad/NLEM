#include "local_minimization.h"

local_minimization::local_minimization ( double (*fptr)(double [], int), int n )
{
	int i;

	if ( n > MXV )
		mabort("Too many variables.","local_minimization");
	detail.funmin = fptr;
	detail.nod = n;

	machine_constants.comacc = DBL_EPSILON;
	machine_constants.bignum = DBL_MAX;
	machine_constants.dwarf = DBL_MIN;
	machine_constants.dwalog = log(machine_constants.dwarf);
	machine_constants.blim = 1.0e20;

	itnocs = 0;

	for (i=0; i<n; i++) {
		detail.ixat[i] = 1;
		detail.marg[i] = 0;
		detail.xll[i] = -machine_constants.bignum;
		detail.xrl[i] = machine_constants.bignum;
		detail.point[i] = 0.0;
	}
}
