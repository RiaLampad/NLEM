#include "local_minimization.h"

void local_minimization::vassign ( double v1[], double v2[], int n )
//  Performs the vector assignment: v1 = v2
{
	int i;

	for (i=0; i<n; i++)
		v1[i] = v2[i];
}
