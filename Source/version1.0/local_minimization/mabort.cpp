#include "local_minimization.h"

void local_minimization::mabort ( const char msg[], const char sub[] )
{
	cout << endl;
	cout << "**** Routine \"" << sub << "\" aborted." << endl;
	cout << "**** " << msg << endl;
	abort();
}
