//----------------------------------------------
//class header "matrix" for 2x2 complex matrices
//----------------------------------------------


class matrix{
 public:
  complex a,b,c,d;
  matrix();
  matrix(complex,complex,complex,complex);
  void set(complex,complex,complex,complex);
  matrix operator+(matrix);
  matrix operator-(matrix);
  matrix operator*(matrix);
  matrix operator*(double);
  matrix operator*(complex);
};

complex det(matrix);
matrix inv(matrix);
