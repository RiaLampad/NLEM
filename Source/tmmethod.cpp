//------------------------------------
//class "layer" for a dielectric layer
//functions for reflection calculation
//------------------------------------
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "complex.h"
#include "matrix.h"
#include "dielectric_models.h"
#include "tmmethod.h"
using namespace std;


void layer::setThickness(double t, int nwv){
  thickness = t;
  for(int i=0; i<nwv; i++){
    phase[i] = exp(((cmplx(0.0, 1.0)*kappa[i])*thickness));
    invphase[i] = inv(phase[i]);
  }
}
//------------------------------------------------------------------------
void layer::setIndex1(int i){
  index1[i] = index[i];
}
//------------------------------------------------------------------------
void layer::setIndex0(int i){
  index0[i] = index[i];
}
//------------------------------------------------------------------------
void layer::makeIndex(double lamda, complex z, int i){
  epsil[i] = z;
  index[i] = sqrt(z);
  kappa[i] = index[i] * (2.0 * acos(-1.0) / lamda);
  invindex[i] = inv(index[i]);
  phase[i] = exp( ((cmplx(0.0, 1.0)*kappa[i])*thickness) );
  invphase[i] = inv(phase[i]);
}
//------------------------------------------------------------------------
//create the interface matrix bewteen current layer and the next 
matrix interface(layer *a, int i){
  complex one = cmplx(1.0, 0.0);
  complex ratio = a[0].index[i] * a[1].invindex[i];
  complex cp = (one + ratio)*0.5;
  complex cm = (one - ratio)*0.5;
  matrix m(cp,cm,cm,cp);
  return m;
}
//------------------------------------------------------------------------
//create propagation matrix within current layer
matrix propagation(layer *a, int i){
  complex zero = cmplx(0.0, 0.0);
  matrix m(a[0].phase[i], zero, zero, a[0].invphase[i]);
  return m;
}
//------------------------------------------------------------------------
//calculate reflection from entire stack
//stack is an array of n layers
double reflection(layer *a, int n, int i){
  matrix m;
    for(int j=0; j < n-2; j++){
    m = interface(&a[j], i) * m;
    m = propagation(&a[j+1], i) * m;
  }
  m = interface(&a[n-2], i) * m;
  return abssq(m.c * inv(m.d));

  //   return abssq((a[0].index[i]-a[1].index[i])*inv(a[0].index[i]+a[1].index[i]));

}
//
//-------------------------------------------------------------------------------
void dataset::setupLayer(int j, double **x, int n){
  //j:   layer whose index we store
  //x:   array with wavelengths and index from input file
  //n:   number of points red from the file
  for(int i=0; i<nwv; i++){
    complex c = findindex(lamda[i], x, n);
    stack[j].makeIndex(lamda[i], c*c, i);
    stack[j].setIndex0(i);
  }
}
//
//----------------------------------------------------------------------
void dataset::setupReflectance(double **x, int n){
  for(int i=0; i<nwv; i++)
    Rmeasured[i] = findrefl(lamda[i], x, n);
}
//
//-----------------------------------------------------------------------
void dataset::setupMinimizer(double **x, int n){
  nparam = n;
  for(int i=0; i<n; i++){
    param[i] = x[0][i];
    lowerlim_param[i] = x[1][i];
    upperlim_param[i] = x[2][i];
  }
}
//
//-----------------------------------------------------------------------
void dataset::calcReflection(){
  for(int i=0; i<nwv; i++) Rfitted[i] = reflection(&stack[0], nlr, i);
  //  for(int i=0; i<nlr; i++) cout<<i<<" "<<stack[i].index[nwv-1].x<<endl;
}
//
//----------------------------------------------------------------------
double dataset::getError(){
  double error = 0.0;
    for(int i=0; i<nwv; i++) error += pow(Rmeasured[i] - Rfitted[i],2);
  //   /abs(Rmeasured[i]);
  
  return error;
}
//
//----------------------------------------------------------------------
void dataset::print_reflectance(){
  for(int i=0; i<nwv; i++) 
    cout << lamda[i]<<" "<<Rmeasured[i] <<" "<<Rfitted[i]<<endl;
}

void dataset::write_reflectance(){

  FILE *refl;
  refl = fopen("results_reflectance", "w");

  fprintf(refl, "%10s\t%10s\t%10s\n", "lambda[nm]", "R_FDTD", "R_NLEM");
  
  for(int i=0; i<nwv; i++)
    fprintf(refl, "%10.4f\t%10.4f\t%10.4f\n", lamda[i], Rmeasured[i], Rfitted[i]);

  fclose(refl);
}

//
//----------------------------------------------------------------------
void dataset::print_index(int j){
  for(int i=0; i<nwv; i++)
    cout << lamda[i]<<" "<<Rmeasured[i] <<" "<<Rfitted[i]<<" "<<
      real(stack[j].index[i]) <<" "<<imag(stack[j].index[i])<<endl;
}
//
//-----------------------------------------------------------------------
complex findindex(double lamda, double **x, int n){
  //find in arrays x the data point corresponding to lamda
  //return the corresponding index of refraction

  int imin = findposition(lamda, &x[0][0], n);
  double f = (lamda-x[0][imin])/(x[0][imin+1]-x[0][imin]);
  double yf = x[1][imin] + f * (x[1][imin+1]-x[1][imin]);
  double zf = x[2][imin] + f * (x[2][imin+1]-x[2][imin]);
  complex c = cmplx(yf,zf);
  return c;
}
//
//-----------------------------------------------------------------------
double findrefl(double lamda, double **x, int n){
  //find in arrays x the data point corresponding to lamda
  //return the corresponding reflection

  int imin = findposition(lamda, &x[0][0], n);

  double f = (lamda-x[0][imin])/(x[0][imin+1]-x[0][imin]);
  return x[1][imin] + f * (x[1][imin+1]-x[1][imin]);
}
//
//-------------------------------------------------------------
int findposition(double lamda, double *x, int n){
  //find in array x the data point corresponding 
  //to lamda, and return its position

  int imin = 0; int imax = n-1;

  while (imax > imin+1){
    int icent = (imin+imax)/2;
    if(lamda < x[icent]) imax = icent;
    else imin = icent;
  }

  return imin;
}
