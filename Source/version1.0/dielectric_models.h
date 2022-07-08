//--------------------------------------------------
//functions for dielectric function models available
//--------------------------------------------------
complex lrntz(double, double, double *);
complex lorentz(double, double *, int, int);
void prepare_lorentz(double *, int, int);
void print_lorentz(double *, int, int);

complex cauchy(double, double *, int, int);
void print_cauchy(double *, int, int);

