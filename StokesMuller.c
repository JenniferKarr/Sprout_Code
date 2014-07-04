/* File: StokesMuller.c 


    Sprout is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sprout is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sprout.  If not, see <http://www.gnu.org/licenses/>.
 
   Copyright 2013, Michihiro Takami, Hyosun Kim, Jennifer Karr
   All rights reserved.   
   Contact: Jennifer Karr (jkarr@asiaa.sinica.edu.tw)


*/
#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"

const double pi, rad2deg;
double g, albedo, S11[181], S12[181], S33[181], S34[181], probability_func[201], PSI;
void VectorRotate(double new_xyz[3], double xyz[3], double xyz0[3]);

double scatterAngle_HG() {
/* WARNING : check if 0 < g < 1 */
/* WARNING : check if 0 < g < 1 */
/* WARNING : check if 0 < g < 1 */
/* WARNING : check if 0 < g < 1 */
/*                                                                       */
/*  DESCRIPTION                                                          */
/*     phase function P(theta;g)                                         */
/*       the Henyey-Greenstein phase function for 0.15 < g <= 1          */
/*       the Rayleigh scattering phase function for 0 <= g < 0.15        */
/*  INPUTS                                                               */
/*     (g as a global variable)                                          */
/*     (P as a random number)                                            */
/*  OUTPUT                                                               */
/*     cos(theta)                                                        */
/*                                                                       */
/*                                              2011-08-24  Hyosun Kim   */

  double g2,A,B,C,D,cosT,P;

  P=genrand_real2();
  if (g > 0.15) {/* Henyey-Greenstein phase function for 0.15 < g <= 1 */
    g2 = g*g;
    A = 1.0 + g2;
    B = 1.0 - g2;
    C = 1.0 - g + 2.0*g*P;
    D = B / C;
    cosT = ( A - D*D ) / ( 2.0 * g );
  } else {     /* Rayleigh scattering phase function for 0 <= g < 0.15 */
    B = 1.0 - 2.0*P;
    C = pow(sqrt(4.*B*B+1.)-2.*B , 1./3.);
    cosT = C - 1.0/C;
  }
  return cosT;
}

double scatterAngle_table_lookup(){
/*                                                                       */
/*  DESCRIPTION                                                          */
/*     phase function P(theta;g)                                         */
/*       the table-lookup method                                         */
/*  INPUTS                                                               */
/*     (probability_func[201])                                           */
/*     (P as a random number)                                            */
/*  OUTPUT                                                               */
/*     cos(theta)                                                        */
/*                                                                       */
/*                                              2011-09-17  Hiro Takami  */	
/*                                                                       */
/*     Note: check completed on 9/19/2011                                */

  double cosT,P,i,i_residual;
  int    i_int;
	
  P          = genrand_real2();
  i          = P*200; // number of elements - 1 for "probability func"
  i_int      = i;
  i_residual = i-i_int*1.;
  cosT       = probability_func[i_int]*(1.-i_residual)+probability_func[i_int+1]*i_residual;
	
  return cosT;
}


double calVelo(double new_vel[3],double vel[3],double cosPSI,double IQUV[4]) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*                                                                       */
/*  INPUTS                                                               */
/*     vel[3]                                                            */
/*     cosPSI                                                            */
/*     IQUV[4]                                                           */
/*     (S11[181], S12[181])                                              */
/*     (2 random numbers for PHI)                                        */
/*  OUTPUTS                                                              */
/*     new_vel[3]                                                        */
/*     PHI                                                               */
/*                                              2011-08-24  Hyosun Kim   */

  double CRIT, A,B,C, dx_norm, dy_norm, dz_norm;
  double PHI, sinPHI, cosPHI;
  //double PSI, sinPSI;
  double sinPSI;
  int intPSI;

  dx_norm=vel[0];  dy_norm=vel[1];  dz_norm=vel[2];

  PHI = genrand_real2();
  PHI = pi * (PHI + (PHI>0.5)*0.5);
/*   if (PHI < 0.5) { */
/*     PHI = pi * PHI; */
/*   } else { */
/*     PHI = pi * (PHI + 0.5); */
/*   } */

  PSI    = acos(cosPSI);   if (PSI > pi) PSI = 0.0;
  intPSI = (int) (PSI * rad2deg);
  CRIT   = (IQUV[1]*cos(2.*PHI)+IQUV[2]*sin(2.*PHI))/IQUV[0];
  CRIT   = (1.+S12[intPSI]/S11[intPSI]*CRIT)*0.5;

  if (genrand_real2() > CRIT) PHI = PHI + pi*0.5;

  sinPSI = sqrt(1.-cosPSI*cosPSI);
  sinPHI = sin(PHI);
  cosPHI = cos(PHI);

  A = sinPSI * cosPHI;
  B = sinPSI * sinPHI;
  C = sqrt(1.-dz_norm*dz_norm);

  if (abs(dz_norm) == 1) {
    new_vel[0] = A;
    new_vel[1] = B;
    new_vel[2] = cosPSI*dz_norm;
  } else {
    new_vel[0] = ( A*dz_norm*dx_norm - B*dy_norm ) / C + cosPSI*dx_norm;
    new_vel[1] = ( A*dz_norm*dy_norm + B*dx_norm ) / C + cosPSI*dy_norm;
    new_vel[2] =  -A*C + cosPSI*dz_norm;
  }
  return PHI;
}


void StokesMuller(double new_IQUV[4], double new_vel[3], 
		  double IQUV[4], double vel[3]) {
/* WARNING : check what do W and W2 mean */
/* WARNING : check what do W and W2 mean */
/* WARNING : check what do W and W2 mean */
/* WARNING : check what do W and W2 mean */
/*                                                                       */
/*  DESCRIPTION                                                          */
/*     RPO * XMUL * ROP                                                  */
/*     = | 1   0    0   0 | * | s11 s12 0    0  | * | 1    0     0   0 | */
/*       | 0  cosW sinW 0 |   | s12 s11 0    0  |   | 0  cosW2 sinW2 0 | */
/*       | 0 -sinW cosW 0 |   |  0   0 s33 -s34 |   | 0 -sinW2 cosW2 0 | */
/*       | 0   0    0   1 |   |  0   0 s34  s33 |   | 0    0     0   1 | */
/*                                                                       */
/*     IQUV_out = (RPO * XMUL * ROP) * IQUV_in                           */
/*     IQUV_out *= I_in / I_out * albedo                                 */
/*                                                                       */
/*  INPUTS                                                               */
/*     IQUV[4]                                                           */
/*     vel[3]                                                            */
/*     (S11[181], S12[181], S33[181], S34[181], albedo)                  */
/*  OUTPUTS                                                              */
/*     new_IQUV[4]                                                       */
/*     new_vel[3]                                                        */
/*                                                                       */
/*                                              2011-08-30  Hyosun Kim   */

  int i, ANGDIFF;
  double W, W2, cosW,sinW,cosW2,sinW2, cosPSI, vel_modified[3];
  double I1,Q1,U1,V1, I2,Q2,U2,V2, s11,s12,s33,s34, scale; 

  /* cosPSI */
  /*  The probability function is replaced on 9/19/2011  */
  //cosPSI = scatterAngle_HG();
  cosPSI = scatterAngle_table_lookup();

  /* W2 (= 2 PHI) */
  W2 = 2. * calVelo(new_vel, vel, cosPSI, IQUV);  /* get new_vel */

  /* W (= 2 THETA) */
  VectorRotate(vel_modified, vel, new_vel);  /* get vel_modified */
  scale = sqrt(pow(vel_modified[1],2) + pow(vel_modified[2],2));
  if (scale == 0.) {
    W = 0.;
  } else {
    W = 2. * acos(-vel_modified[1] / scale);
  }
  if (W > 2.*pi) W = 0.;
  if (vel_modified[2] < 0.) W = -W;

  /* cosW, sinW, cosW2, sinW2 */
  cosW = cos(W);   cosW2 = cos(W2);
  sinW = sin(W);   sinW2 = sin(W2);

  /* s11, s12, s33, s34 */
  ANGDIFF=(int) (acos(cosPSI)*rad2deg + 0.5);  if (ANGDIFF>180) ANGDIFF=180;
  s11=S11[ANGDIFF];  s12=S12[ANGDIFF];  s33=S33[ANGDIFF];  s34=S34[ANGDIFF];

  /* new_IQUV  via.  matrix calculation */
  I1 = IQUV[0];
  Q1 = IQUV[1] * cosW2 + IQUV[2] * sinW2;
  U1 =-IQUV[1] * sinW2 + IQUV[2] * cosW2;
  V1 = IQUV[3];

  I2 = s11 * I1 + s12 * Q1;
  Q2 = s12 * I1 + s11 * Q1;
  U2 = s33 * U1 - s34 * V1;
  V2 = s34 * U1 + s33 * V1;

  new_IQUV[0] = I2;
  new_IQUV[1] = Q2 * cosW + U2 * sinW;
  new_IQUV[2] =-Q2 * sinW + U2 * cosW;
  new_IQUV[3] = V2;

  scale = I1 / I2 * albedo;
  for (i=0; i<4; i++) {
    new_IQUV[i] *= scale;
  }
}
