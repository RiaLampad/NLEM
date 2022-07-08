#include "local_minimization.h"

double local_minimization::get_minimizer ( int i )
{
	if ( i < 0 || i >= detail.nod )
		mabort("No such variable","get_minimizer");
	return detail.point[i];
}

double local_minimization::get_lower_value ( )
{
	return detail.val;
}

int local_minimization::get_function_counter ( )
{
	return itnocs;
}

void local_minimization::get_gradient ( double grad[] )
{
	int noc, ngr;

	dermer(detail.point,detail.val,grad,&noc,&ngr);
}
