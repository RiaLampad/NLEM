#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <cstdlib>
#include <iomanip>

using namespace std;

#define MXV 50                     //  Maximum number of variables.

class local_minimization {

private:
	int itnocs;

//  Machine constants.
	struct {
		double comacc;
		double bignum;
		double dwarf;
		double dwalog;
		double blim;
	} machine_constants;

//  Function details.
	struct {
		int nod;
		double (*funmin)(double [], int);
		int ixat[MXV];
		int marg[MXV];
		double xll[MXV];
		double xrl[MXV];
		double point[MXV];
		double val;
	} detail;

//  Vector/matrix operations.
	void vassign ( double [], double [], int );
	double vvmult ( double [], double [], int );
	void lmult ( double *, double *, double *, int, int );
	double llmult ( double *, double [], int, int );
	void ltmult ( double *, double [], double [], int, int );
	void llcomp ( double *, int , int );
	void llsolv ( double *, double [], double [], int, int );

//  Termination criteria.
	int igconv ( double [], double [], int, double, double );
	int ifconv ( double, double, double );
	int ixconv ( double [], double [], int, double );

//  Derivatives.
	void dermer ( double [], double, double [], int *, int * );
	void cderi ( double [], int, double, double, double *, int * );
	double acsq1 ( double [], int, double );

//  Development aids.
	void print_hessian ( const char [], double *, int, int );
	void print_vector ( const char [], double [], int );

//  Utility.
	double axp ( double );
	int freeva ( void );
	double acsq ( double [], int );
	void mdis ( double [], double, int, int, int, int, int );
	void mabort ( const char [], const char [] );

//  For the TRUST minimizer.
	void trusmv ( double *, double [], double [], int, double *, int, double, double, 
                    int, double, double, double, int *, int *, int *, int *, int, int [],
	              int, double [], double [], double [], double [], double [], int, double *, int [] );
	void ttdri ( double *, int, double [], double *, int, double, double *, double, int *,
                   double [], int, int [], double [], double [], double [], int * );
	void trupd ( double [], int, double, double, double *, int, double *, double, int *,
	             double *, double [], double [], int, int [], double [], int );
	void dogleg ( double, double [], double [], int, double [], double, double, double, 
                    double, double, double );
	void hpdd ( double *, double [], int, int, double, int * );
	void lldec ( double *, double [], int, int, double, double *, double );
	void jrot ( double *, int, double, double, int, int );
	void llupd ( double *, double *, double *, int, int );
	double fcnt ( double [], int );
	void trtm ( int, char * );

//  Transformations.
	void xtoy ( double [], int );
	void ytox ( double [], int );
	void gtrans ( double [], double [] );
	void ltrans ( double *, int, double [] );

//  Minimization routines.
	void trust ( int );

public:
	local_minimization ( double (*)(double [], int), int );

//  "Set" functions.
	void set_variable ( int, double );
	void set_lower_bound ( int, double );
	void set_upper_bound ( int, double );

//  "Get" functions.
	double get_minimizer ( int );
	double get_lower_value ( void );
	void get_gradient ( double [] );
	int get_function_counter ( void );

//  Minimization routines.
	void run ( int );

};

#define idx(i,j,nr) ((nr)*(j)+(i))
#define sign(a,b) ( (b)>=0 ? abs(a) : -abs(a) )

void granal ( int, double [], double [] );
