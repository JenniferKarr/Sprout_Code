/* File: ccodes.i */
/*  9/1/2011   Created (Hyosun)                                             */
/*  3/3/2012   Revisions are made for the Shakura-Sunyaev function (Hiro)   */
/*                                                                          */
/*                                       Hyosun Kim and Hiro Takami         */

%module ccodes
%include "carrays.i"
%array_functions(double, darr);

%{
void VectorRotate(double new_xyz[3], double xyz[3], double xyz0[3]);
void StokesMuller(double new_IQUV[4], double new_UVW[3], 
		  double IQUV[4], double UVW[3]);

double tau_per_AU_Shakura_Sunyaev(double x,double y,double z);

int NewPosition(double new_xyz[3], double xyz[3], double vel[3]);
int NewPosition4givenTau(double tau_goal,
                double new_xyz[3], double xyz[3], double vel[3]);
void Initialize(double IQUV[4], double xyz[3], double vel[3]);
void init_genrand(unsigned long s);
%}

%inline %{
  const double pi=3.1415926535897931;
  const double rad2deg=57.295779513082323;

  /* 2/22/2012: dR is added   */
  /* 3/3/2012 : dR is removed */
  //extern double g, albedo, Rmin, Rmax, Tdisk, tau, dR;
  extern double g, albedo, Rmin, Rmax, Tdisk, tau, dpath, PSI;
  extern double Shakura_Sunyaev_tau_0,alpha,beta,h_0,stellar_radius_in_AU; 
  extern double S11[181], S12[181], S33[181], S34[181], probability_func[201];

/*   double *new_a(int size) { */
/*     return (double (*)) malloc(size*sizeof(double)); */
/*   } */
/*   void free_a(double (*a)) { */
/*     free(a); */
/*   } */
/*   double set_a(double a[], int i, double val){ */
/*     a[i] = val; */
/*   } */
/*   double get_a(double a[], int i){ */
/*     return a[i]; */
/*   } */

  double (*new_matrix(int size))[91] {
    return (double (*)[91]) malloc(91*size*sizeof(double));
  }
  void set_matrix(double x[][91], long i, long j, double val) {
    x[i][j] = val;
  }
/*   double get_matrix(double x[301][91], long i, long j) { */
/*     return x[i][j]; */
/*   } */
  void free_matrix(double (*x)[91]) {
    free(x);
  }

/* double rand1(double min, double max) { */
/*   /\* random number generator between min and max *\/ */
/*   return min + (max-min)*( rand()/(RAND_MAX+1.0) ); */
/* } */
%}

void VectorRotate(double new_xyz[3], double xyz[3], double xyz0[3]);
void StokesMuller(double new_IQUV[4], double new_UVW[3], 
		  double IQUV[4], double UVW[3]);

double tau_per_AU_Shakura_Sunyaev(double x,double y,double z);

int NewPosition(double new_xyz[3], double xyz[3], double vel[3]);
int NewPosition4givenTau(double tau_goal,
                double new_xyz[3], double xyz[3], double vel[3]);
void Initialize(double IQUV[4], double xyz[3], double vel[3]);
void init_genrand(unsigned long s);

/* REFERENCES */
/* http://www.swig.org/Doc1.3/SWIG.html#SWIG_nn15                            */
/* http://www.swig.org/Doc1.3/Python.html#Python_nn19                        */
/* http://www.swig.org/Doc1.3/Library.html#Library_carrays                   */
/* http://www.cplusplus.com/doc/tutorial/program_structure                   */
/* http://www.exforsys.com/tutorials/c-language/c-structures-and-unions.html */
/* http://www.acm.uiuc.edu/webmonkeys/book/c_guide                           */
