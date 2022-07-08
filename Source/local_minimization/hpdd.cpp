#include "local_minimization.h"

void local_minimization::hpdd ( double *a, double adiag[], int ld, int n, double comacc, int *ispd )
//  ---------------------------------------------------------------------
//
//  Description:
//    Subroutine HPDD finds a non-negative constant m, such that
//    A + m*I is safely positive definite, where m=0 if A is already
//    positive definite. It then sets A = A+m*I and calculates the
//    Choleski decomposition of the new A.
//    Actualy only the lower triangular part of A is referenced,
//    while the Cholesky factor L[t] is stored in the upper
//    triangular part of A.
//    For details on the implementation see:
//
//      Dennis J.E & Schnabel R.B.
//      Numerical Methods for Unconstrained Optimization and Nonlinear
//      Equations, Prentice-Hall, 1983.
//
//  Input arguments:
//    LD         Leading dimension of matrix A.
//    N          Dimensionality of matrix A.
//    COMACC     Machine accuracy.
//
//  Output arguments:
//    ADIAG      Diagonal elements of the matrix A+m*I.
//    ISPD       Indicates whether the original matrix was positive
//               definite:
//               ISPD = 0 -> The original matrix A was not positive
//                           definite.
//               ISPD = 1 -> The original matrix was positive definite.
//
//  Input / Output arguments:
//    A          On input an NxN symmetrix matrix. Only the diagonal and
//               lower triangular elements are referenced.
//               On output the diagonal and upper triangular part
//               contains the Choleski factor L[t].
//
//  ---------------------------------------------------------------------
{
	int i, j;
	double dma, dmi, aii, pdm, xmu, offmax, con, ema, emi, soff, comtol;
	double sdd, addmax;
	double zero = 0.0, one = 1.0;

	comtol = sqrt(comacc);
	for (i=0; i<n; i++)
		adiag[i] = a[idx(i,i,ld)];

	dma = adiag[0];
	dmi = adiag[0];
	for (i=1; i<n; i++) {
		aii = adiag[i];
		dma = max(dma,aii);
		dmi = min(dmi,aii);
	}
	pdm = max(zero,dma);
	if (dmi <= comtol*pdm) {
		xmu = 2*(pdm-dmi)*comtol-dmi;
		dma = dma+xmu;
	} else
		xmu = 0.0;
	offmax = 0.0;
	for (i=0; i<n; i++)
		for (j=0; j<=i-1; j++)
			offmax = max(offmax,abs(a[idx(j,i,ld)]));
	con = one+2*comtol;
	if (offmax*con > dma) {
		xmu = xmu-dma+offmax*con;
		dma = dma*con;
	}
	if (dma == 0.0) {
		xmu = 1.0;
		dma = 1.0;
	}
	*ispd = 1;
	if (xmu > 0.0) {
		*ispd = 0;
		for (i=0; i<n; i++)
			adiag[i] = adiag[i]+xmu;
	} 
	offmax = sqrt(max(dma,offmax/n));

	lldec(a,adiag,ld,n,offmax,&addmax,comacc);
	if (addmax > 0) {
		ema = adiag[0];
		emi = adiag[0];
		for (i=0; i<n; i++) {
			soff = 0;
			for (j=0; j<=i-1; j++)
				soff = soff+abs(a[idx(i,j,ld)]);
			for (j=i+1; j<n; j++)
				soff = soff+abs(a[idx(j,i,ld)]);
			ema = max(ema,adiag[i]+soff);
			emi = min(emi,adiag[i]-soff);
		}

		sdd = (ema-emi)*comtol-emi;
		sdd = max(zero,sdd);
		xmu = min(addmax,sdd);
		for (i=0; i<n; i++)
			adiag[i] = adiag[i]+xmu;
		lldec(a,adiag,ld,n,zero,&addmax,comacc);
		*ispd = 0;
	}
}
