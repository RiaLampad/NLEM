#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cstdio>

#include "complex.h"
#include "matrix.h"
#include "dielectric_models.h"
#include "local_minimization/local_minimization.h"
#include "tmmethod.h"


using namespace std;

dataset dset;

double funmin(double*, int);
complex make_rough_eps(complex, complex, double, double);
complex mg_emt(complex, complex, double, double);
double get_rough_filling(double, int, int);
double get_p_emt_slice(double, double);
double get_p_emt_semi(double, double);
double get_slope(double ,double, int, double);
double get_weight(complex , double );
void readInput();
void readMaterial(int);
void readReflectance();
int smoothFile(double**, int, int);
int check_lamdas(double, double, double*, int);
void sortFile(double**, int, int);
void swap(double*, double*);
void set_minimization(local_minimization*);
void reset_minimization(local_minimization*);
void store_new_optimal(local_minimization*);
double max(double, double);
double min(double, double);


int main(){
  double error, errormin; int noc;
  double thick, dthick, thick_min, thick_max;
  
  

  readInput();
  
  for(int i=0; i<dset.nlr; i++) readMaterial(i);
  readReflectance();  

  
  if(dset.ToPrint==1){
    cout <<endl<< " ANALYSIS AND FITTING OF REFLECTANCE:"<<endl<<
      "Fitting from " << dset.lamda[0] <<" to " <<
      dset.lamda[dset.nwv-1] <<" nm using " << dset.nwv << " points" <<endl;
    

  }
  
  //---------------------------------------------------------------------------------------//
  
  local_minimization *a;
  a = new local_minimization(&funmin, dset.nparam);

  set_minimization(a);
  errormin=funmin(dset.param, dset.nparam);
  
  errormin=1000000.0;
  noc = 0;

  for(int irun=0; irun<dset.Nruns; irun++){
  
     a->run(10000);
     error = a->get_lower_value();
    
    if(error < errormin){
      errormin = error;
      store_new_optimal(a);
      cout << "error "<< errormin <<" in run "<< irun << endl;
    }

    noc  = a->get_function_counter();
    reset_minimization(a);
  }
  delete a;
  //-------------------------------------------

  if(dset.ToPrint==1){
    cout << "the final error is " << errormin << ", reached after " << noc <<
      " minimizer calls" << endl <<endl <<"Fit results:"<<endl;
    
    error = funmin(dset.param1,dset.nparam);
    if(dset.FitIndex == 1) dset.print_index(dset.jfit2);
    else dset.write_reflectance();
    cout << endl << "Fitted paramaters: "<<endl;
    
    int ic = 0;


    if(dset.FitRoughness == 1){
      cout << "Roughness height [nm]: " <<dset.param1[ic++] << endl;
      cout << "Roughness shape:       " <<dset.param1[ic++] << endl;

    }

  }


}


//-------------------------------------------------------

double funmin(double *param, int nparam){
  complex e, eh, e1;
  double L, height, thick, shape, weight, porous, f, w, slope;
  
  L = 32.0; //***//

  int j1 = dset.jfit1;
  int j2 = dset.jfit2;
  
  int ic = 0;
  
  if(dset.FitRoughness == 1){
    height = param[ic++];
    thick  = height / dset.Rlrs;
    shape  = param[ic++];
    for(int j=j1+1; j<j2; j++) dset.stack[j].setThickness(thick, dset.nwv);
    
  }
 
  
  if(dset.FitRoughness == 1){
    for(int i=0; i<dset.nwv; i++){
      eh = dset.stack[j1].index0[i] * dset.stack[j1].index0[i];
      if(dset.FitPorocity == 1) e =  dset.stack[j2].index1[i] * dset.stack[j2].index1[i];
      else e =  dset.stack[j2].index0[i] * dset.stack[j2].index0[i];
      
      for(int j=j1+1; j<j2; j++){
	f      = get_rough_filling(shape, j-j1-1, j2-j1-1);
        slope  = get_slope(shape,L, j-1, height);	
	weight = get_weight(e,  slope);
	e1 = make_rough_eps(e, eh, f,weight);
	dset.stack[j].makeIndex(dset.lamda[i], e1, i);
      }
    }
  }
   
  dset.calcReflection();
  
  return dset.getError(); 
}

//----------------------------------------------------------------

complex make_rough_eps(complex e, complex eh, double f, double weight){
  double n = sqrt(e).x;
  double p_slice = get_p_emt_slice(f, n);
  double p_semi  = get_p_emt_semi(f, n);
  complex e_slice = mg_emt(e, eh, f, p_slice);
  complex e_semi  = mg_emt(e, eh, f, p_semi);
  return e_slice * weight + e_semi * (1.0-weight);
}  
//-----------------------------------------------------------

//------------------------------------------------------------

complex mg_emt(complex e, complex eh, double f, double p){
    complex z1, z2;
    z1 = eh * ((e * (f*p + 1.0)) + (eh * (p - f*p)));
    z2 = inv((e * (1.0-f)) + (eh * (p+f)));
    return z1 * z2;
  }
//-----------------------------------------------------------
//-------------------------------------------------

double get_rough_filling(double shape, int  j, int jm){
  double y = (j + 0.5) / jm;
  return pow(y, 2.0/shape);
  }

//----------------------------------------------------------------

double get_p_emt_slice(double f, double n){
  double p = dset.emt_p[0] + dset.emt_p[1] * n + dset.emt_p[2] * n*n;
  double g = dset.emt_g[0] + dset.emt_g[1] * n + dset.emt_g[2] * n*n;
  double d = dset.emt_d[0] + dset.emt_d[1] * n + dset.emt_d[2] * n*n;
  return p + g * pow((pow(f,d) - 0.5), 2);
}
//----------------------------------------------------------------

double get_weight(complex e, double slope){
  double n = sqrt(e).x;
  double a = dset.weight_a[0] + dset.weight_a[1]*n;
  double b = dset.weight_b[0] + dset.weight_b[1]*n;
  return a + b * slope;
}

//------------------------------------------------------------------

double get_p_emt_semi(double f, double n){
  return dset.emt_semi[0] + dset.emt_semi[1] * n + dset.emt_semi[2] * n*n;
}

//-------------------------------------------------------------------
double get_slope(double shape,double L, int j, double height){
  double y2,y1,dy;
  double x2,x1,dx;
  double a;

  a = height/pow((L/2.),shape);
  y1 = j+0.5;
  y2 = j+1+0.5;
  x1 = pow( (y1/a),1/shape );
  x2 = pow( (y2/a),1/shape ); 
  dy = y2-y1;
  dx = x2-x1;

  return atan(dx/dy)*180/acos(-1.);
}


//-------------------------------------------------------------------

void readInput(){
  string ch, ch1, line, name;  
  int layer = 0;
  int ip = 0;
  dset.NLorentz = dset.NCauchy = dset.NTaucLor = 0;
  
  ifstream input("input");
  while(!input.eof()){
      input>>ch;
      
      if(ch.at(0) == '#'){
	getline(input, line);
      }
      else if(ch == "REFLECTANCE"){
	input >> dset.reffile;
	getline(input, line);
      }
      
      else if(ch == "MATERIAL"){
	input >> dset.stack[layer].filename;
	input >> dset.stack[layer].thickness;
	input >> dset.stack[layer++].IfToFit;
	getline(input, line);
      }

      else if(ch == "FIT_TYPE"){
	input >> ch1;
	if (ch1 == "ROUGHNESS") input >> dset.FitRoughness;
	getline(input, line);
      }

      else if(ch == "FIT_LIMITS"){
	input >> ch1;
	if (ch1 == "RUNS_NUM")  input >> dset.Nruns;
	else if (ch1 == "ROUGH_LRS") input >> dset.Rlrs;
	else if (ch1 == "LAMDA_MIN") input >> dset.lamda_min;
	else if (ch1 == "LAMDA_MAX") input >> dset.lamda_max;
	else if (ch1 == "LAMDA_NUM") input >> dset.nwv;
	getline(input, line);
      }


      else if(ch == "FIT_PARAM"){
	input >> ch1;	
	if (ch1 == "ROUGH"){
	  if(dset.FitRoughness == 1)
	    input >> dset.param[ip] >> dset.lowerlim_param[ip] >> dset.upperlim_param[ip++];
	}
      }

      else if(ch == "EMT_PARAM"){
	input >> ch1;
	if (ch1 == "PMAX")
	  input >> dset.emt_p[0] >> dset.emt_p[1] >> dset.emt_p[2];
	else if (ch1 == "GAMM")
	  input >> dset.emt_g[0] >> dset.emt_g[1] >> dset.emt_g[2];
	else if (ch1 == "DELT")
	  input >> dset.emt_d[0] >> dset.emt_d[1] >> dset.emt_d[2];
	else if (ch1 == "SEMI")
	  input >> dset.emt_semi[0] >> dset.emt_semi[1] >> dset.emt_semi[2];
	getline(input, line);
      }
      //__________________________________________________________________________//
      else if (ch == "WEIGHT_PARAM"){
	input >> ch1;
	if (ch1 == "BETA")
	  input >> dset.weight_b[0] >> dset.weight_b[1] >> dset.weight_b[2];
	else if (ch1 == "ALPHA")
	  input >> dset.weight_a[0] >> dset.weight_a[1] >> dset.weight_a[2];
	getline(input, line);
      }
      else if(ch == "OPTION"){
	input >> ch1;
	if (ch1 == "PRINT")
	  input >> dset.ToPrint;
	else if (ch1 == "SMOOTH")
	  input >> dset.ToSmooth;
	getline(input, line);
      }	
      else{
	cout << "ATTENTION: " << ch << " not proper input" << endl;
	getline(input, line);
      }
  }

  input.close();

  double omega1 = 1./dset.lamda_max;
  double omega2 = 1./dset.lamda_min;

  for(int i=0; i<dset.nwv; i++){

    double omega  = omega2 + (omega1-omega2)*i/(dset.nwv-1);
    dset.lamda[i] = 1./omega;
  }
  
  dset.nparam = ip;
  for(int ip=0; ip<dset.nparam; ip++) dset.param1[ip] = dset.param[ip];
  
  dset.nlr = layer;
  int isum = 0;
  for(int j=0; j<layer; j++){
    isum += dset.stack[j].IfToFit;
    if(dset.stack[j].IfToFit == 1) dset.jfit1 = j-1;
  }

  if(isum == 0) cout <<"WARNING!!!! No layer chosen for fit" << endl;
  if(isum > 1) cout  <<"WARNING!!!! Too many layers chosen for fit" << endl;

  dset.jfit2 = dset.jfit1 + 1;
  
  if(dset.FitRoughness == 1){

    for(int j=layer-1; j>dset.jfit1; j--){
	dset.stack[j+dset.Rlrs].filename  = dset.stack[j].filename;
	dset.stack[j+dset.Rlrs].thickness = dset.stack[j].thickness;
	dset.stack[j+dset.Rlrs].IfToFit   = dset.stack[j].IfToFit;   
    }

    dset.jfit2 = dset.jfit1 + dset.Rlrs + 1;

    for(int j=dset.jfit1+1; j<dset.jfit2; j++){
      dset.stack[j].filename  = dset.stack[dset.jfit2].filename;
      dset.stack[j].thickness = dset.stack[dset.jfit2].thickness * 0.0;
      dset.stack[j].IfToFit   = dset.stack[dset.jfit2].IfToFit; 
    }
    dset.nlr += dset.Rlrs; 
  }
    
}

//-------------------------------------------------------------------

void readMaterial(int j){
  string line; 

  if(j <= dset.jfit1 || j >= dset.jfit2)
    cout<<"Loading data for material "<<dset.stack[j].filename<<endl;
  
  double x[3][1000]; double *xx[3];
  for(int k=0; k<3; k++) xx[k] = &x[k][0];
    
  ifstream input(dset.stack[j].filename.c_str());

  int i=0;
  while(!input.eof() && i<1000) input >> x[0][i] >> x[1][i] >> x[2][i++]; getline(input, line);
  int n = i-1;
  
  if(i >= 1000) cout <<" ERROR!!! file "<<dset.stack[j].filename<<" too large"<<endl;
  
  input.close();
  sortFile(xx, n, 3);

  int lamdaOK = check_lamdas(dset.lamda_min, dset.lamda_max, &x[0][0], n);

  if(lamdaOK == 0){
    dset.setupLayer(j, xx, n);
  }
  else {
    cout<<"ERROR!!! "<<dset.stack[j].filename<<" does not cover requested wavelengths"<<endl;
    cout<<dset.stack[j].filename<<" is limited from "<<x[0][0]<<" to "<<x[0][n-1]<<endl;
    exit(1);
  }
  
}

//------------------------------------------------------------------

void readReflectance(){
  string line;
  
  cout<<"Loading refelctance data from file "<<dset.reffile<<endl;

  double x[2][5000]; double *xx[2];
  for(int k=0; k<2; k++) xx[k] = &x[k][0];
  
  ifstream input(dset.reffile.c_str());

  int i=0;
  while(!input.eof() && i<5000) input >> x[0][i] >> x[1][i++]; getline(input, line);
  int n = i-1;

  if(i >= 5000) cout <<" ERROR!!! file "<<dset.reffile<<" too large"<<endl;

  input.close();
  sortFile(xx, n, 2);
  
  int lamdaOK = check_lamdas(dset.lamda_min, dset.lamda_max, &x[0][0], n);
  if(lamdaOK == 0){
    if(dset.ToSmooth > 1) n = smoothFile(xx, n, 2);      
    dset.setupReflectance(xx, n);
  }
  else {
    cout<<"ERROR!!! "<<dset.reffile<<" does not cover requested wavelengths"<<endl;
    cout<<dset.reffile<<" is limited from "<<x[0][0]<<" to "<<x[0][n-1]<<endl;
    exit(1);
  }
}

//----------------------------------------------------------------

int smoothFile(double **x, int n, int kc){

  int i, j, k, ii, nn;
  double *sum;
  sum = new double [kc];

  ii = 0;
  nn = 0;
  int nnew = n / dset.ToSmooth;
  for (j=0; j<nnew; j++){
    
    for(k=0; k<kc; k++) sum[k] = 0;
    for(i=0; i<dset.ToSmooth; i++){    
      for(k=0; k<kc; k++) sum[k] += x[k][ii];
      ii++;
    }

    for(k=0; k<kc; k++) x[k][nn] = sum[k] / dset.ToSmooth;
    nn++;
  }
  delete[] sum;
  return nnew;
}

//---------------------------------------------------------------


//----------------------------------------------------------------------
int check_lamdas(double lmin, double lmax, double *x, int n){

  if( x[0] > lmin || x[n-1] < lmax)
    return 1;
  else
    return 0;
}

//-----------------------------------------------------------------------

void sortFile(double **x, int n, int kc){
  for(int j=n-1; j>0; j--)
    for(int i=0; i<j; i++)
      if(x[0][i] > x[0][i+1])
	for(int k=0; k<kc; k++) swap(&x[k][i],&x[k][i+1]);
}

//-----------------------------------------------------------------------

void swap(double *x, double *y){
  double z=*x; *x = *y; *y = z;
}

//------------------------------------------------------------------------

void set_minimization(local_minimization *a){
  for(int i=0; i<dset.nparam; i++){
    a->set_variable(i,dset.param[i]);
    a->set_lower_bound(i, dset.lowerlim_param[i]);
    a->set_upper_bound(i, dset.upperlim_param[i]);
  }
}

//----------------------------------------------------------------------

void reset_minimization(local_minimization *a){
  double x, y;
  for(int i=0; i<dset.nparam; i++){
    x = (rand()%1000)/1000.0;
    y = dset.lowerlim_param[i] + x * (dset.upperlim_param[i]-dset.lowerlim_param[i]);
    a->set_variable(i, y);
    dset.param[i] = y;
  }
}

//-----------------------------------------------------------------------

void store_new_optimal(local_minimization *a){
  for(int i=0; i<dset.nparam; i++)
    dset.param1[i] = a->get_minimizer(i);
}

//-----------------------------------------------------------------------

double max(double x, double y){
  return (x>y) ? x : y;
}

//-------------------------------------------------------------------

double min(double x, double y){
  return (x<y) ? x : y;
}

//------------------------------------------------------------------
