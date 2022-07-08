//-------------------------------------------------
//class "matrix" functions for 2x2 complex matrices
//-------------------------------------------------
#include <cmath>
#include "complex.h"
#include "matrix.h"

matrix::matrix(){
  a.set(1.0, 0.0); d = a;
  b.set(0.0 ,0.0); c = b;
}

matrix::matrix(complex a1, complex b1, complex c1, complex d1){
  a = a1;
  b = b1;
  c = c1;
  d = d1;
}

void matrix::set(complex a1, complex b1, complex c1, complex d1){
  a = a1;
  b = b1;
  c = c1;
  d = d1;
}

matrix matrix::operator+(matrix v){
  matrix u;
  u.a = a + v.a;
  u.b = b + v.b;
  u.c = c + v.c;
  u.d = d + v.d;
  return u;
}

matrix matrix::operator-(matrix v){
  matrix u;
  u.a = a - v.a;
  u.b = b - v.b;
  u.c = c - v.c;
  u.d = d - v.d;
  return u;
}

matrix matrix::operator*(matrix v){
  matrix u;
  u.a = (a*v.a) + (b*v.c);
  u.b = (a*v.b) + (b*v.d);
  u.c = (c*v.a) + (d*v.c);
  u.d = (c*v.b) + (d*v.d);
  return u;
}

matrix matrix::operator*(complex v){
  matrix u;
  u.a = a*v;
  u.b = b*v;
  u.c = c*v;
  u.d = d*v;
  return u;
}

matrix matrix::operator*(double v){
  matrix u;
  u.a = a*v;
  u.b = b*v;
  u.c = c*v;
  u.d = d*v;
  return u;
}

complex det(matrix v){
  return (v.a*v.d) - (v.b*v.c);
}

matrix inv(matrix v){
  matrix u;
  u.a=v.d;
  u.d=v.a;
  u.b.x=-v.b.x;
  u.b.y=-v.b.y;
  u.c.x=-v.c.x;
  u.c.y=-v.c.y;
  complex de = inv(det(v));
  u.a = u.a * de;
  u.b = u.b * de;
  u.c = u.c * de;
  u.d = u.d * de;
  return u;
}

