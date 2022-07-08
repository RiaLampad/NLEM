//------------------------------------------
//class header "complex" for complex numbers
//------------------------------------------
class complex{
 public:
  double x,y;
  complex();
  complex(double, double);
  void set(double, double);
  complex operator+(complex); 
  complex operator-(complex); 
  complex operator*(complex); 
  complex operator*(double); 
};

complex cmplx(double, double);
double abs(complex);
double abssq(complex);
complex conjg(complex);
double real(complex);
double imag(complex);
complex exp(complex);
complex inv(complex);
complex sqrt(complex);

