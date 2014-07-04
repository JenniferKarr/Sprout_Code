/* File: VectorRotate.c 


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

void VectorRotate(double new_xyz[3], double xyz[3], double xyz0[3]) {
/*                                                                       */
/*  DESCRIPTION                                                          */
/*     Rotate the vector (x,y,z) in the sense that the                   */
/*     other (normal) vector (x0,y0,z0) goes to (1,0,0),                 */
/*     given x0^2+y0^2+z0^2=1.                                           */
/*       |       x0        y0    z0 | | x0 |   | 1 |                     */
/*       |    -sinT      cosT     0 | | y0 | = | 0 |                     */
/*       | -dz*cosT  -dz*sinT  sing | | z0 |   | 0 |                     */
/*                                                                       */
/*       where  sing = sqrt(x0^2+y0^2) = sqrt(1-z0^2)                    */
/*              cosT = x0/sing                                           */
/*              sinT = y0/sing                                           */
/*  INPUTS                                                               */
/*     xyz[3]                                                            */
/*     xyz0[3]                                                           */
/*  OUTPUT                                                               */
/*     new_xyz[3]                                                        */
/*                                                                       */
/*                                              2011-07-31  Hyosun Kim   */

  double x,y,z,x0,y0,z0;
  double theta,cosT,sinT,sing;
  
  x=xyz[0]; y =xyz[1]; z=xyz[2]; x0=xyz0[0]; y0=xyz0[1]; z0=xyz0[2];

  theta = atan2(y0,x0); /* if z0=1 (x0=y0=0), theta=0 */
  cosT = cos(theta);
  sinT = sin(theta);
  sing = sqrt(1.-z0*z0);

  new_xyz[0] = x*x0  + y*y0 + z*z0;
  new_xyz[1] = y*cosT - x*sinT;
  new_xyz[2] = z*sing - (x*cosT + y*sinT)*z0;
}
