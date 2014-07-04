/* File: Initialize.c 


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
#include <math.h>
#include "mt19937ar.h"

const double pi;
double Tdisk;

void Initialize(double IQUV[4], double xyz[3], double vel[3]) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*     initialize (I,Q,U,V), (x,y,z), (dx_norm, dy_norm, dz_norm)        */
/*     for each photon                                                   */
/*  INPUTS                                                               */
/*     (Tdisk)                                                           */
/*     (2 random numbers for direction vector vel[3])                    */
/*  OUTPUTS                                                              */
/*     IQUV[4]                                                           */
/*     xyz[3]                                                            */
/*     vel[3]                                                            */
/*                                                                       */
/*                                              2011-09-01  Hyosun Kim   */
  double velz, len, ang;
  
  IQUV[0]=1.;  IQUV[1]=0.;  IQUV[2]=0.;  IQUV[3]=0.;  
  xyz[0]=0.;  xyz[1]=0.;  xyz[2]=0.;

  velz= (2.*genrand_real2()-1.) * sin(Tdisk);
  ang = (2.*genrand_real2()-1.) * pi;
  len = sqrt( 1. - velz*velz );
  vel[0] = len * cos(ang);
  vel[1] = len * sin(ang);
  vel[2] = velz;
}
