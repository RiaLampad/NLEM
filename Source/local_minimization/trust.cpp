#include "local_minimization.h"

void local_minimization::trust ( int noc )
{
	int info, ihessu, iprint, iradu, itmax, nocex, nogex, iter;
	int nod;
	double v[MXV][MXV];
	double xtmp[MXV];
	double xt1[MXV], xt2[MXV], xt3[MXV], xt4[MXV], xt5[MXV], xt6[MXV];
	double xtol, ftol, gtol, delta, rho, comacc;
	char mess[100];

//  Set default values for some panel parameters and process any
//  changes the user enters.
	ftol = 0.0;
	xtol =  machine_constants.comacc;
	gtol = 0.0;
	rho = 1.0e-4;

	itmax = -1;
	ihessu = 0;
	iradu = 0;
	iprint = 0;

	delta = 0.0;
	iradu = 0;
	comacc = machine_constants.comacc;
	nod = detail.nod;
	detail.val = acsq(detail.point,nod);

//  Transform the variables.
	vassign(xtmp,detail.point,nod);
	xtoy(xtmp,nod);

//  Call the minimizer.
	trusmv((double *)v,xtmp,xt1,nod,&detail.val,noc,comacc,gtol,itmax,xtol,
	       rho,ftol,&info,&iter,&nocex,&nogex,iprint,detail.ixat,
	       ihessu,xt2,xt3,xt4,xt5,xt6,iradu,&delta,detail.marg);

//  Inverse transformation.
	ytox(xtmp,nod);
	vassign(detail.point,xtmp,nod);
	detail.val = acsq(detail.point,nod);

//  Construct a message describing the reason for termination.
	trtm(info,mess);

//  Report the reason for termination.
	if ( iprint !=0 )
		cout << "  TRUST: " << info << " " << mess << endl;

}
