//------------------------------------------------------
//class "color" for chromaticity analysis of reflection
//------------------------------------------------------
#include <cmath>
#include <iostream>
#include "color.h"
using namespace std;


void color::makeIllumination(){
  //calculate the black body spectrum at temperature t. 
  //Standard solar spectrum is at 6500K.
  //Illumination spectrum comes normalized with peak value 1.

  double lamda, c1, c2, bb;
  double bbmax = 0;
  double t = illuminationTemp; 

  for(int i=0; i<cie_num; i++){

    lamda = (cie_min + i * cie_del) * 1e-9;
    c1 = 3.74177152*1e-16;  //c1 = 2*pi*h*c^2
    c2 = 1.43877696*1e-2;   //c2 = h*c/k
    bb = c1/( pow(lamda,5) * (exp(c2/(lamda*t)) - 1.0) );

    bbmax = (bb > bbmax) ? bb : bbmax;
    Illumination[i] = bb;
  }

  for(int i=0; i<cie_num; i++) Illumination[i] *= 1.0/bbmax;

}

void color::getColor(double *lamda, double *refl, int n){
  
  double lam, ref, ill, dlam, f, X0, Y0, Z0, cmx,cmy,cmz;
  double x1,y1,z1, R0,G0,B0;
  int j;

  X=Y=Z=X0=Y0=Z0=Tsr=0;

  //integration to get the XYZ tristimulus values and Tsr
  for (int i=0; i<n-1; i++){
    
    lam = 0.5 * (lamda[i] + lamda[i+1]);
    ref = 0.5 * (refl[i] + refl[i+1]);
    dlam = abs(lamda[i] - lamda[i+1]);

    //get the Tsr
    Tsr += ref * dlam;

    //get the XYZ tristimulus values
    if(lam >= cie_min && lam <= cie_min+cie_num*cie_del){
      
      j = (int) ((lam-cie_min)/cie_del);      
      f = (lam - (cie_min+j*cie_del))/cie_del;

      cmx = ColorMatching[0][j] + f * (ColorMatching[0][j+1]-ColorMatching[0][j]);
      cmy = ColorMatching[1][j] + f * (ColorMatching[1][j+1]-ColorMatching[1][j]);
      cmz = ColorMatching[2][j] + f * (ColorMatching[2][j+1]-ColorMatching[2][j]);
      ill = Illumination[j]    + f * (Illumination[j+1]-Illumination[j]);

      //XYZ tristimulus values for the reflected light
      X  += cmx * ref * ill * dlam;
      Y  += cmy * ref * ill * dlam;
      Z  += cmz * ref * ill * dlam;
      //XYZ tristimulus values for the illuminating light
      X0 += cmx * ill * dlam;
      Y0 += cmy * ill * dlam;
      Z0 += cmz * ill * dlam;
    }
  }
  Tsr *=  2.0 / abs(lamda[0]+lamda[1]-lamda[n-1]-lamda[n-2]);
  
  //get the La*b* from the XYZ tristimulus values
  y1 = Y/Y0; 
  if(y1 > 0.008856) L = 116 * pow(y1, 1.0/3.0) - 16;
  else L = 903.3 * y1;
  
  x1= X/X0;
  if(x1 > 0.008856) a = 500.0 * (pow(x1, 1.0/3.0) - pow(y1, 1.0/3.0));
  else a = 3893.5 * (x1 - y1);
  
  z1 = Z/Z0;
  if(z1 > 0.008856) b = 200.0 * (pow(y1, 1.0/3.0) - pow(z1, 1.0/3.0));
  else b = 1557.4 * (y1 - z1);

  //get the RGB (Adobe1998) from the XYZ tristimulus values
  double xyz2rgb[3][3];
  xyz2rgb[0][0] =  2.0414;
  xyz2rgb[1][0] = -0.9693;
  xyz2rgb[2][0] =  0.0134;
  xyz2rgb[0][1] = -0.5649;
  xyz2rgb[1][1] =  1.8760;
  xyz2rgb[2][1] = -0.1184;
  xyz2rgb[0][2] = -0.3447;
  xyz2rgb[1][2] =  0.0416;
  xyz2rgb[2][2] =  1.0154;

  //unscaled RGB coordinates of the reflected light
  R = xyz2rgb[0][0]*X + xyz2rgb[0][1]*Y + xyz2rgb[0][2]*Z;
  G = xyz2rgb[1][0]*X + xyz2rgb[1][1]*Y + xyz2rgb[1][2]*Z;
  B = xyz2rgb[2][0]*X + xyz2rgb[2][1]*Y + xyz2rgb[2][2]*Z;
  
  //unscaled RGB coordinates of the incident light
  R0 = xyz2rgb[0][0]*X0 + xyz2rgb[0][1]*Y0 + xyz2rgb[0][2]*Z0;
  G0 = xyz2rgb[1][0]*X0 + xyz2rgb[1][1]*Y0 + xyz2rgb[1][2]*Z0;
  B0 = xyz2rgb[2][0]*X0 + xyz2rgb[2][1]*Y0 + xyz2rgb[2][2]*Z0;

  //scale RGB so that max of R0G0B0 is 255.
  double rgbmax = 0;
  rgbmax = (R0 > G0) ? R0 : G0;
  rgbmax = (B0 > rgbmax) ? B0 : rgbmax;
  R *= 255 / rgbmax;
  G *= 255 / rgbmax;
  B *= 255 / rgbmax;
  R = (R > 0) ? R : 0;
  G = (G > 0) ? G : 0;
  B = (B > 0) ? B : 0;

  //get the chromaticity coordinates x,y from XYZ
  x = X/(X+Y+Z);
  y = Y/(X+Y+Z);

}

void color::setupCIE(double **xx, int n){
  double lam;

  if(n != cie_num) cout<<"ERROR in reading CIE"<<endl;

  for (int i=0; i<n; i++){
    lam = cie_min + i*cie_del;
    if(lam != xx[0][i]) cout<<"ERROR in reading CIE "<<lam<<" "<< xx[0][i]<<endl;
    for(int j=0; j<3; j++) ColorMatching[j][i] = xx[j+1][i];
  }
}

double color::Get_x(){return x;}
double color::Get_y(){return y;}
double color::Get_L(){return L;}
double color::Get_a(){return a;}
double color::Get_b(){return b;}
double color::Get_R(){return R;}
double color::Get_G(){return G;}
double color::Get_B(){return B;}
double color::Get_Tsr(){return Tsr;}
