//---------------------------------------------
//class "complex" functions for complex numbers
//---------------------------------------------
#include <cmath>
#include "complex.h"

complex::complex(){ x = y = 0;} 

complex::complex(double x1, double y1){  
  x = x1;
  y = y1;
}

void complex::set(double x1, double y1){ 
  x = x1;
  y = y1;
}

complex complex::operator+(complex v){
  complex u;
  u.x = x + v.x;
  u.y = y + v.y;
  return u;
}

complex complex::operator-(complex v){
  complex u;
  u.x = x - v.x;
  u.y = y - v.y;
  return u;
}

complex complex::operator*(complex v){
  complex u;
  u.x = x*v.x - y*v.y;
  u.y = x*v.y + y*v.x;
  return u;
}

complex complex::operator*(double m){
  complex u;
  u.x = x * m;
  u.y = y * m;
  return u;
}

complex cmplx(double x1, double y1){
  complex u;
  u.x = x1;
  u.y = y1;
  return u;
}

double abs(complex v){
  return sqrt(v.x*v.x + v.y*v.y);
}

double abssq(complex v){
  return v.x*v.x + v.y*v.y;
}

complex conjg(complex v){
  complex u;
  u.x = v.x;
  u.y = -v.y;
  return u;
}

double real(complex v){
  return v.x;
}

double imag(complex v){
  return v.y;
}

complex exp(complex v){
  complex u;
  double r = real(v);
  double i = imag(v);
  double e = exp(r);
  u.x = e * cos(i);
  u.y = e * sin(i);
  return u;
}

complex inv(complex v){
  complex u = conjg(v)*(1.0/abssq(v));
  return u;
}

complex sqrt(complex v){
  double r = abs(v);
  double x = real(v);
  complex u = cmplx(sqrt(0.5*(r+x)),sqrt(0.5*(r-x)));
  return u;
} 
