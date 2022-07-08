//-------------------------------------------------------
//classes and functions for color analysis of reflectance
//-------------------------------------------------------
#include <cstring>
#include <cstdlib>

using namespace std;

#define cie_num 471    //CIE: number of wavelengths
#define cie_min 360    //CIE: minimum wavelength
#define cie_del 1      //CIE: delta wavelength

class color{
 public:
  double illuminationTemp;
  string file;
  double X,Y,Z,L,a,b,R,G,B,x,y,Tsr;
  double ColorMatching[3][cie_num]; //cie color matching functions
  double Illumination[cie_num];    //bb illumination of temperature t 
  void makeIllumination();
  void setupCIE(double **, int);
  void getColor(double *, double *, int);
  double Get_x();
  double Get_y();
  double Get_L();
  double Get_a();
  double Get_b();
  double Get_R();
  double Get_G();
  double Get_B();
  double Get_Tsr();
};







