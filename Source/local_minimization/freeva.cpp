#include "local_minimization.h"

int local_minimization::freeva ( )
//  ---------------------------------------------------------------------
//
//  Description:
//    Determines the number of free (non-fixed) variables.
//
//  ---------------------------------------------------------------------
{
	int i, nfree;

	nfree = 0;
	for (i=0; i<detail.nod; i++)
		if (detail.ixat[i] != 0)
			nfree++;
	return nfree;
}
