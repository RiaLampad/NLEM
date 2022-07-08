#include "local_minimization.h"

void local_minimization::set_variable ( int i, double x )
{
	if ( i < 0 || i >= detail.nod )
		mabort("No such variable","set_point");
	detail.point[i] = x;
}

void local_minimization::set_lower_bound ( int i, double b )
{
	if ( i < 0 || i >= detail.nod )
		mabort ("No such variable","set_left_bound");
	detail.xll[i] = b;
      if (detail.marg[i] == 0)
		detail.marg[i] = -1;
	else if (detail.marg[i] == 1)
		detail.marg[i] = 2;
}

void local_minimization::set_upper_bound ( int i, double b )
{
	if ( i < 0 || i >= detail.nod )
		mabort ("No such variable","set_right_bound");
	detail.xrl[i] = b;
      if (detail.marg[i] == 0)
		detail.marg[i] = 1;
	else if (detail.marg[i] == -1)
		detail.marg[i] = 2;
}
