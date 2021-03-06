//---------------------------------------------------------
//classes and functions for fitting of optical reflectance
//based on the Fresnel equations formulated with the 
//trasfer matrix method. Normal incidence considered only. 
//---------------------------------------------------------
#include <cstring>
#include <cstdlib>

using namespace std;

#define NXV 50
#define nwvs 5000    //max number of wavelengths
#define nlrs 100     //max number of layers

//class "layer" for a dielectric layer
class layer{
 public:
  string filename;      //filename containing wavelength and index
  double thickness;       //layer thickness in nm
  int IfToFit;
  complex epsil[nwvs];
  complex index0[nwvs];
  complex index1[nwvs];
  complex index[nwvs];    //complex refractive index
  complex invindex[nwvs]; //inverse of refractive index
  complex kappa[nwvs];    //complex propagation constant index*w/c
  complex phase[nwvs];    //propagation phase aquired in layer
  complex invphase[nwvs]; //inverse propagation phase
  void setThickness(double, int);
  void makeIndex(double, complex, int);
  void setIndex0(int);
  void setIndex1(int);
};

//class "dataset" for the data to be fitted
class dataset{
public:
int Nruns, Rlrs;
int FitThickness, FitRoughness, FitPorocity, FitIndex;
int ToPrint, ToSmooth;
int NLorentz, NCauchy, NTaucLor; 
double lamda_min, lamda_max;
int nlr;                    //number of layers in the stack
int nwv;                    //number of wavelengths to be considered
int jfit1, jfit2;                   //layer to be fitted
int nparam;
string reffile;           //file with experimental reflectance data
string modelname;         //name of dielectric model to use
double lamda[nwvs];         //wavelengths to be considered
double Rmeasured[nwvs];     //measured reflectance
double Rfitted[nwvs];       //fitted reflectance
layer stack[nlrs];          //stack of dielectric layers to be fitted
double param1[NXV];
double param[NXV];
double lowerlim_param[NXV];
double upperlim_param[NXV];
double emt_p[3], emt_g[3], emt_d[3], emt_semi[3];
double weight_a[3], weight_b[3];
void setupLayer(int, double **, int);
void setupReflectance(double **, int);
void setupMinimizer(double **, int);
void calcReflection();
double getError();
void print_reflectance();
void write_reflectance();
void print_index(int);
};

//functions for reflection calculation
matrix interface(layer *, int);
matrix propagation(layer *, int);
double reflection(layer *, int, int);

//functions for the setting up process
complex findindex(double, double **, int);
double findrefl(double, double **, int);
int findposition(double, double *, int);

