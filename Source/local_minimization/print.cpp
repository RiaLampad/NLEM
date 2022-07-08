#include "local_minimization.h"

void local_minimization::print_hessian ( const char msg[], double *h, int n, int ld )
{
	int i, j;

	cout << fixed << setprecision(12);
	for (i=0; i<n; i++) {
		cout << msg << " ";
		for (j=i; j<n; j++)
			 cout << setw(15) << h[idx(i,j,ld)] << " ";
		cout << endl;
	}
}

void local_minimization::print_vector ( const char msg[], double v[], int n)
{
	int i;

	cout << fixed << setprecision(12);
	cout << msg << " ";
	for (i=0; i<n; i++)
		cout << v[i] << " ";
	cout << endl;
}
