//--------------------------------------------------
//functions for available dielectric function models
//--------------------------------------------------
#include <cmath>
#include <iostream>
#include "complex.h"
#include "dielectric_models.h"
using namespace std;

//-----------------------------------------------------
//LORENTZ model functions------------------------------

complex lrntz(double w, double dw, double *f){
  //w: photon energy (eV), f[0-2]: gamma, w0, sigma
  complex u = cmplx(pow(f[1]+dw,2) - pow(w,2), -f[0]*w);
  return inv(u) * pow(f[2],2);
}

complex lorentz(double lamda, double *f, int n, int nlor){
  //f contains list with n adjustable parameters:
  //f[ic]=e00, f[ic+1]=g1, f[ic+2]=Dw0,  f[ic+3]=wp1, etc
  //  int nlor = (n-ic-1)/3; //the number ot lorentzians
  int ic = n - 3*nlor - 1;
  double w = 1239.84193 / lamda;
  double dw = 0.0;
  complex e = cmplx(f[ic], 0.0);
  for(int i=0; i<nlor; i++){
    e = e + lrntz(w, dw, &f[i*3+ic+1]);
    dw += f[i*3+ic+2];
  }
  return e;
}

void prepare_lorentz(double *f, int n, int nlor){
  //change so that oscilator frequencies appear as a Dw
  //  int nlor = (n-ic-1)/3; 
  int ic = n - 3*nlor - 1;
  double dw = 0.0;
  for(int i=0; i<nlor; i++){
    f[i*3+ic+2] -= dw;
    dw += f[i*3+ic+2];
  }
}

void print_lorentz(double *f, int n, int nlor){
  //int nlor = (n-ic-1)/3;
  int ic = n - 3*nlor - 1;
  double dw = 0.0;
  cout <<"   -e_infinity: " << f[ic] << endl << 
    "   -lorentzians (g, w, f):" << endl;
  
  for(int i=0; i<nlor; i++){
    cout << "     lorentzian " <<  i+1 << ": " << 
      f[i*3+ic+1] << " " << f[i*3+ic+2]+dw << " " << f[i*3+ic+3] << endl ;
    dw += f[i*3+ic+2];
  }
  cout<<endl;
}

//END of LORENTZ model functions-------------------------------

complex cauchy(double lamda, double *f, int n, int ncau){
  //f contains list with n adjustable parameters:
  //f[0]=thickness, f[1]=B, f[2]=C, etc
  int ic = n - ncau;
  double indx = f[ic];
  for (int i=ic+1; i< n; i++) indx += f[i]/pow(lamda, 2.0*(i-1));   
  return cmplx(indx*indx, 0.0);
}
  
void print_cauchy(double *f, int n, int ncau){
  int ic = n - ncau;
  cout <<"   -cauchy parameters: ";
  for(int i=ic; i<n; i++) cout << f[i] << " "; 
  cout << endl<<endl;;
}


