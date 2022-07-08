#include "local_minimization.h"

void local_minimization::trtm ( int info, char *mess )
//  ---------------------------------------------------------------------
//
//  Description:
//    Constructs a message, describing the reason for termination of the
//    TRUST method.
//
//  Input arguments:
//    INFO       Reason for termination of the TRUST method, as returned
//               by subroutine TRUSMV.
//
//  Output arguments:
//    MESS       Message, describing the reason for termination.
//
//  ---------------------------------------------------------------------
{
	if (info == 1)
		strcpy(mess,"Target value has been reached");
	else if (info == 2)
		strcpy(mess,"The gradient criterion is satisfied");
	else if (info == 3)
		strcpy(mess,"All function evaluations have been used");
	else if (info == 4)
		strcpy(mess,"The minimization variables have converged");
	else if (info == 5)
		strcpy(mess,"The function value has converged");
	else if (info == 6)
		strcpy(mess,"The Hessian matrix cannot be updated");
	else if (info == 7)
		strcpy(mess,"The specified number of iterations has been reached");
	else if (info == 8)
		strcpy(mess,"All variables are fixed");
	else if (info == 9)
		strcpy(mess,"Further progress is not possible");
	else {
		cout << "trtm: Incorrect \"info\"" << endl;
		exit(1);
	}
}
