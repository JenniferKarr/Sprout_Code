/* File: NewPosition.c 



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
#define NEVERHAPPEN 1
#undef NEVERHAPPEN

#ifdef NEVERHAPPEN
#include <stdio.h>
#endif

#include <stdlib.h>
#include <math.h>
#include "mt19937ar.h"

const double rad2deg;
double Rmin, Rmax;       //  2/22/2012: their unit has been changed to AU
//double dR;             //  2/22/2012: added for the revision this time
                         //  3/03/2012: removed when we directly introduce the SS function
                         //              without an array.
double tau;
double dpath=0.1;        //  3/12/2012: step for integration tau (AU). Originally in NewPosition4givenTau.

double Shakura_Sunyaev_tau_0,alpha,beta,h_0,stellar_radius_in_AU; 
                         // 3/3/2012: added for the SS function

double tau_per_AU_Shakura_Sunyaev(double x,double y,double z) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*                                                                       */
/*  3/3/2012                                                             */
/*  I am moving the original table look-up alrogirm here. (Hiro)         */
/*                                                                       */
/*  INPUTS                                                               */
/*     x, y, z                                                           */
/*  OUTPUT                                                               */
/*     None (the value tau per AU is returned.)                          */
/*                                                                       */

  double r_in_stellar_radius,h;

  r_in_stellar_radius = sqrt(x*x + y*y) / stellar_radius_in_AU;
  h  = h_0 * pow(r_in_stellar_radius, beta);

  return Shakura_Sunyaev_tau_0 / pow(r_in_stellar_radius,alpha) * exp (-0.5*pow(z/h,2));
  //return Shakura_Sunyaev_tau_0 / pow(r_in_stellar_radius,alpha) * exp (-0.5*z*z/(h*h));

}

int NewPosition(double new_xyz[3], double xyz[3], double vel[3]) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*                                                                       */
/*     9/19/2011                                                         */
/*     I split up Sunny's NewPosition (see description below) into two:  */
/*     i.e., giving tau_goal using a random number here, and do the      */
/*     remaining calculations in the other subroutine, namely,           */
/*     "NewPosition4givenTau".                                           */
/*                                                                       */
/*     --- (Sunny's description                                          */
/*     Get a new position new_xyz[3] from the current position xyz[3],   */
/*     given unit displacement vector vel[3], by calculating optical     */
/*     depth (in situ) of the model nebula extended from Rmin to Rmax    */
/*       Destiny 1 : photon_will_be_out of the entire nebula             */
/*       Destiny 2 : achieved, tau -> tau_goal in the photon path        */
/*                                                                       */
/*      distance2border = X = distance from xyz[3] to a border at r0     */
/*      xyz[3] = (x,y,z)                   = position vector             */
/*      vel[3] = (dx_norm,dy_norm,dz_norm) = unit velocity vector        */
/*      -->  | vec{r} + X vec{dr_norm} |^2 = r0^2                        */
/*      -->  X^2 + 2 B X + C = 0   where  B = vec{r} * vec{dr_norm}      */
/*                                        C = r^2 - r0^2                 */
/*      -->  X = -B +/- sqrt( B^2 - C )                                  */
/*                                                                       */
/*  INPUTS                                                               */
/*     xyz[3]                                                            */
/*     vel[3]                                                            */
/*     disk[][91]                                                        */
/*     (Rmin, Rmax)                                                      */
/*  OUTPUT                                                               */
/*     new_xyz[3]                                                        */
/*     photon_will_be_out                                                */
/*                                                                       */
/*                                              2011-09-01  Hyosun Kim   */
/*                                                                       */
/*  ---                                                                  */
/*                                                                       */
/*  Revision by Hiro                                                     */
/*                                                                       */
/*  9/19/2011                                                            */
/*  Now the correct tau is provided when the photon goes out from the    */
/*   nebula. (The calculation for the last step was skipped in Hyoson's  */
/*   original code.)                                                     */
/*  I also have temporarily set dpath constant (0.3).                    */
/*                                                                       */
/*  2/22/2012                                                            */
/*  I have changed the unit of the coordinate to AU.                     */
/*                                                                       */
/*  3/3/2012                                                             */
/*  The array for opacity distribution ("disk") is removed.              */
/*  We alternatively use the Shakura & Sunyaev function.                 */
/*                                                                       */

  int photon_goes_out;
  double tau_goal;

  tau_goal = -log(genrand_real2());
  photon_goes_out=NewPosition4givenTau(tau_goal,new_xyz,xyz,vel);

  return photon_goes_out;
}


int NewPosition4givenTau(double tau_goal,
                double new_xyz[3], double xyz[3], double vel[3]) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*                                                                       */
/*                                                                       */
/*     Get a new position new_xyz[3] from the current position xyz[3],   */
/*     given unit displacement vector vel[3], by calculating optical     */
/*     depth (in situ) of the model nebula extended from Rmin to Rmax    */
/*       Destiny 1 : photon_will_be_out of the entire nebula             */
/*       Destiny 2 : achieved, tau -> tau_goal in the photon path        */
/*                                                                       */
/*      distance2border = X = distance from xyz[3] to a border at r0     */
/*      xyz[3] = (x,y,z)                   = position vector             */
/*      vel[3] = (dx_norm,dy_norm,dz_norm) = unit velocity vector        */
/*      -->  | vec{r} + X vec{dr_norm} |^2 = r0^2                        */
/*      -->  X^2 + 2 B X + C = 0   where  B = vec{r} * vec{dr_norm}      */
/*                                        C = r^2 - r0^2                 */
/*      -->  X = -B +/- sqrt( B^2 - C )                                  */
/*                                                                       */
/*  INPUTS                                                               */
/*     xyz[3]                                                            */
/*     vel[3]                                                            */
/*     disk[][91]                                                        */
/*     (Rmin, Rmax)                                                      */
/*  OUTPUT                                                               */
/*     new_xyz[3]                                                        */
/*     photon_will_be_out                                                */
/*                                                                       */
/*                                              2011-09-01  Hyosun Kim   */
/*  3/12/2012                                                            */
/*  "dpath" is now a global parameter. See above for more details of     */
/*  revision with the Shakura & Sunyaev function.                        */
/*                                                                       */

  int photon_heads_for_cavity=0;
  int photon_enters_in_cavity=0;
  int photon_goes_out        =0;

  double Rmin2, Rmax2;
  double dx_norm, dy_norm, dz_norm, x, y, z, x_new, y_new, z_new;
  //double dpath, dpath_real, path_length_old=0., path_length_new;
  double dpath_real, path_length_old=0., path_length_new;
  //double dtau, dummy, r2, B, C, D, sqrtD;
  double dtau, r2, B, C, D, sqrtD;
  double distance2border, sol_small, sol_large, destination;
  /* revised for the SS function */
  //int index_r_old, index_t_old, index_r_new, index_t_new;
  double tau_per_AU,tau_per_AU_new;

  tau   = 0.;
  Rmin2 = Rmin * Rmin;
  Rmax2 = Rmax * Rmax;

/***** x, y, z, r2, dx_norm, dy_norm, dz_norm *****/
  dx_norm=vel[0]; dy_norm=vel[1]; dz_norm=vel[2];
  x=xyz[0]; y=xyz[1]; z=xyz[2];   r2=x*x+y*y+z*z;

#ifdef NEVERHAPPEN
  if (r2 > Rmax2) {printf("OUT? (r2=%f)>(Rmax2=%f) ?\n", r2, Rmax2);}
#endif
  
/***** photon_heads_for_cavity, distance2border, destination, sol_large *****/
  /* r=0 is the only case that photon starts from cavity (ref. MC.py) */
  if (r2 == 0.) { /* photon within cavity : modify to start from Rmin */
    /* photon freely moves through inner cavity (x,y,z,r2 modified) */
    x += Rmin * dx_norm;
    y += Rmin * dy_norm;
    z += Rmin * dz_norm;
    r2 = Rmin2;
#ifdef NEVERHAPPEN
    r2 = x*x + y*y + z*z;
    if (abs(r2-Rmin2) > 0.00001) {
      printf("AT INNER BORDER? : (r2=%f) = (Rmin2=%f) ?\n", r2, Rmin2);
    } else {
      r2 = Rmin2;
    }
#endif
    distance2border = Rmax - Rmin;
    destination = distance2border;
    /* photon_heads_for_cavity = 0; (default) */
  } else { /* photon starting between Rmin and Rmax */
    B = x*dx_norm + y*dy_norm + z*dz_norm;
    
    /*** roots for inner boundary (r0 = Rmin) ***/
    C = r2 - Rmin2;
    D = B*B - C;
    if (D > 0.) { /* inner cavity is in photon path (regardless direction) */
      /* D=0 is not this special case because photon just passes the border */
      sqrtD = sqrt(D);
      sol_small = -sqrtD-B;
      sol_large =  sqrtD-B;
      
      /* possibilty of sols=(-,+) has already been excluded */
      /* by shifting photon within the inner cavity to Rmin */
#ifdef NEVERHAPPEN
      if (sol_large * sol_small < 0.){
	printf("BE NOT IN CAVITY! %f * %f > 0 ? \n", sol_small, sol_large);
      }
#endif

      if (sol_small > 0.) { /* sols=(+,+); photon heads for cavity */
	photon_heads_for_cavity = 1;
	distance2border = sol_small; /* closer one from present position */
      } /* else {photon_heads_for_cavity = 0;} (default) */
    } /* else {photon_heads_for_cavity = 0;} (default) */
    
    /*** roots for outer boundary (r0 = Rmax) ***/
    C = r2 - Rmax2;
    D = B*B - C;
    destination = sqrt(D)-B;
    
#ifdef NEVERHAPPEN
    if (D < 0.) {printf("FAR AWAY FROM NEBULA?\n");}
    if (destination < 0.) {printf("(-,-)?\n");}
    if ((-sqrt(D)-B) > 0.) {printf("(+,+)?\n");}
#endif

    if (photon_heads_for_cavity == 0) {
      distance2border = destination;
    }
  } /* end if-else */
  
/***** index_r_old,  index_t_old *****/
/* WARNING : int, +0.5 : maybe need finer grids than 1 deg */
/* WARNING : int, +0.5 : maybe need finer grids than 1 deg */
/* WARNING : int, +0.5 : maybe need finer grids than 1 deg */
/* WARNING : int, +0.5 : maybe need finer grids than 1 deg */
/* WARNING : see also the same calculation in while loop */

///* 2/22/2012: as the unit of r,r2, z is now AU, we should divide them by dR. */
////  index_r_old = (int) (sqrt(r2) + 0.5);
//  index_r_old = (int) (sqrt(r2) / dR + 0.5);
//
//  index_t_old = (int) (atan2(abs(z),sqrt(r2-z*z))*rad2deg + 0.5); 
  tau_per_AU = tau_per_AU_Shakura_Sunyaev(x,y,z);

/***** while loop *****/


  while (1) {

    /*** path_length_new,  dpath_real,  photon_enters_in_cavity ***/
    path_length_new = path_length_old + dpath;
    if (path_length_new > distance2border){ /* photon crosses border */
      if (photon_heads_for_cavity == 1) { /* border is inner cavity */
	photon_enters_in_cavity = 1;
      } else {      /* border is outer boundary of the model nebula */
        photon_goes_out=1;
      }
      path_length_new = distance2border;

    } /* else { photon_enters_in_cavity = 0; } (default) */
    dpath_real = path_length_new - path_length_old;

    /*** index_r_new,  index_t_new ***/
    x_new = x + path_length_new * dx_norm;
    y_new = y + path_length_new * dy_norm;
    z_new = z + path_length_new * dz_norm;

//    dummy = x_new*x_new + y_new*y_new;
//    /* 2/22/2012: as the unit of r,r2, z is now AU, we should divide them by dR. */
////    index_r_new = (int) (sqrt(dummy+z_new*z_new) + 0.5);
//    index_r_new = (int) (sqrt(dummy+z_new*z_new)/dR + 0.5);
//
//    index_t_new = (int) (atan2(abs(z_new),sqrt(dummy))*rad2deg + 0.5); 

    /*** dtau ***/
/* WARNING : should use better way */
/* WARNING : should use better way */
/* WARNING : should use better way */
/* WARNING : should use better way */
//    dtau = (disk[index_r_old][index_t_old]+disk[index_r_new][index_t_new])*0.5
//      * dpath_real;
    tau_per_AU_new=tau_per_AU_Shakura_Sunyaev(x_new,y_new,z_new);
    dtau = (tau_per_AU + tau_per_AU_new) *0.5 * dpath_real;

    /*** achieved goal?  cross inner cavity? ***/
    // revised on 9/19/2011, Hiro Takami
    if (photon_goes_out == 1){
      break;
    } else if (tau+dtau >= tau_goal) { /* achieved */
      path_length_new = path_length_old + (tau_goal-tau) / dtau * dpath_real;
      x_new = x + path_length_new * dx_norm;
      y_new = y + path_length_new * dy_norm;
      z_new = z + path_length_new * dz_norm;
      break;
    } else if (photon_enters_in_cavity == 1) { /* enter in the inner cavity */
      path_length_new = sol_large; /* jump to the opposite border of cavity */
      distance2border = destination;
      photon_heads_for_cavity = 0;
      photon_enters_in_cavity = 0;
    } else { /* just going to the outer boundary */
/* WARNING : why need this ?*/
/* WARNING : why need this ?*/
/* WARNING : why need this ?*/
/* WARNING : why need this ?*/
//      if (dtau > 0.1) {
//	dpath *= 0.5;
//      } else {
//	dpath *= 2.0;
//      }
//      if (dpath > 2.) {
//	dpath = 2.;
//      }
    }

    /*** updates ***/
    tau += dtau;

    /* revised on 3/3/2012 */
    //index_r_old     = index_r_new;
    //index_t_old     = index_t_new;
    tau_per_AU = tau_per_AU_new;

    path_length_old = path_length_new;

  } /* end while */

  /* return */
  new_xyz[0]=x_new;  new_xyz[1]=y_new;  new_xyz[2]=z_new;

  // revised on 9/19/2011, Hiro Takami
  //return 0; /* photon_will_be_out = 0 */
  return photon_goes_out;
}

