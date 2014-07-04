#!/usr/bin/env python

"""


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

  History of revision:

  9/17/2011    array of probability function for scattering angle is added.
  9/20/2011    the particles Wood+02, Model 1 are added.
               ice coating (10%, 50% relative to core radius) is added
  4/3/2011     BHCOAT is replaced as the  wiscomb code.
               The ice coat is removed.
               Instead, size distribution (Cotera x 3), (Cotera x 10) are added.
  10/19/2012   Test with the grid for array_for_probability_function_of_cospsi
               of 2000
  10/29/2012   The probablity function for the models with amorphous was not saved
               correctly because of a typo. This is now corrected.
  11/5/2012    Now the calculations for dust properties are made for 0.1-deg. steps.
               I decide to still save S11, S12, S33, S34 with 1-deg. step, and
               the probability function with 200 steps. The simulations using the saved
               probability function reproduce well the phase function for
               KMH, C01, and C01x15 dusts at 1.65 um. (run test_probability_function.py)
               So I decide to use this  for the RY Tau paper.
  12/13/2012   Calculations are added with amorphous carbon by Jaeger et al. (1998).
               (See also the revision note for Wood02_wiscomb.py and physical_constants.py.)
               Now we can save comments in the output file.
               I also tided up the code as we now calculate a fairly large number of cases.
"""

from pylab import *
import pickle

#################---------------------------------------------------------------
### Constants ###---------------------------------------------------------------
#################---------------------------------------------------------------

dust_library_filename='dust_library.pickle'

d_angle = 0.1

n_bin_angles = int(180/d_angle)+1

###########################-----------------------------------------------------
### Classes and Modules ###-----------------------------------------------------
###########################-----------------------------------------------------

class Load_dust_data:
  def __init__(self,dust,wavelength_in_micron):
    self.dust                 = dust
    self.wavelength_in_micron = wavelength_in_micron

    opened=open(dust_library_filename)
    dust_library=pickle.load(opened)
    opened.close()

    self.S11=dust_library[dust][wavelength_in_micron]['S11']
    self.S12=dust_library[dust][wavelength_in_micron]['S12']
    self.S33=dust_library[dust][wavelength_in_micron]['S33']
    self.S34=dust_library[dust][wavelength_in_micron]['S34']

    self.probability_func=dust_library[dust][wavelength_in_micron]['probability function for scattering']


    self.kappa_ext=dust_library[dust][wavelength_in_micron]['kappa_ext (cm2 g-1)']
    self.kappa_sca=dust_library[dust][wavelength_in_micron]['kappa_sca (cm2 g-1)']

    self.g     =dust_library[dust][wavelength_in_micron]['g']
    self.albedo=dust_library[dust][wavelength_in_micron]['albedo']

    self.dust_to_H_mass_ratio=dust_library[dust][wavelength_in_micron]['dust to H mass ratio']

#----------------------------------------------------------------------------------
theta        = array(range(n_bin_angles))*d_angle
sint         = sin(radians(theta))
cost         = cos(radians(theta))
cost_reverse = list(cost)
cost_reverse.reverse()
cost_reverse = array(cost_reverse)

def array_for_probability_function_of_cospsi(S11,steps=200):
  """
    Provide the probability function for the scattering angle (cos psi) as an array.
    Originally developed in test_new_scattering_ag_generator3.py.
    The 1-D array with 201 elements (default) is provided. The element i shows
    to cos_psi corresponding to a random number (0-1) P*200.

    input parameter:  S11
    output parameter: integrated probability function to give cos (psi)
                      (for a random number 0-1, with an interval 0.005 for default)

                                                9/17/2011  Hiro Takami
  """

  S11_reverse=list(S11)
  S11_reverse.reverse()
  S11_reverse=array(S11_reverse)

  interp=interp1d(cost_reverse,S11_reverse)

  integ  = [0.]
  cospsi = (array(range(200001))-100000)*0.00001

  for i in range(len(cospsi)-1):
    d_integ=(interp(cospsi[i])+interp(cospsi[i+1]))*0.5
    integ.append(integ[-1]+d_integ)

  integ=array(integ)
  integ/=integ[-1]

  integ_interp=interp1d(integ,cospsi)

  return integ_interp(array(range(steps+1))*1./steps)

#----------------------------------------------------------------------------------
def shrink_array(ar):
  """
    This subroutine is made to shrink the S11/S12/S33/S34 arrays
    to 1-deg. step (thereby 181 elements).
  """
  n           = len(ar)
  step        = int(1/d_angle)

  final_array = []

  for i in range(0,n,step):
    final_array.append(ar[i])

  return array(final_array)

#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

### Begin (to save the data) ###

if __name__ == "__main__":

  from scipy.interpolate import interp1d
  import KMH94_wiscomb  as K
  import Wood02_wiscomb as W

  dust_list=['KMH94, R_V=3.1',
               'Cotera+01, amorphous',
               'Cotera+01 x 15, amorphous',
               'Cotera+01, graphite',
               'Cotera+01 x 15, graphite',

               'Cotera+01, amorphous (J98,400C)',
               'Cotera+01, amorphous (J98,600C)',
               'Cotera+01, amorphous (J98,800C)',
               'Cotera+01, amorphous (J98,1000C)',
               'Cotera+01 x 15, amorphous (J98,400C)',
               'Cotera+01 x 15, amorphous (J98,600C)',
               'Cotera+01 x 15, amorphous (J98,800C)',
               'Cotera+01 x 15, amorphous (J98,1000C)'

               #'Wood+02 Model 1, amorphous',
               #'Wood+02 Model 1, graphite',
               ]

  wl_list=[0.55,1.25,1.65,2.2]

  ### Initialize the parameters

  dust_library={}
  params={}
  for dust in dust_list:
    dust_library[dust]={}
    params[dust]      ={}
    for wl in wl_list:
      dust_library[dust][wl]={}
      params[dust][wl]      ={}

  ### get and write parameters

  d_save=open(dust_library_filename,'w')

  for wl in (0.55,1.25,1.65,2.2):

    print
    print "### Wavength = %.2f um ###" % wl
    print

    ### no ice ###
    params['KMH94, R_V=3.1'][wl]      = K.Set_particle_params_KMH94(wavelength_in_micron=wl)

    print 'Cotera+01, amorphous'
    params['Cotera+01, amorphous'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01, amorphous (J98,400C)'
    params['Cotera+01, amorphous (J98,400C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,400C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01, amorphous (J98,600C)'
    params['Cotera+01, amorphous (J98,600C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,600C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01, amorphous (J98,800C)'
    params['Cotera+01, amorphous (J98,800C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,800C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01, amorphous (J98,1000C)'
    params['Cotera+01, amorphous (J98,1000C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,1000C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01, graphite'
    params['Cotera+01, graphite'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='graphite',
                    a_min=1e-3,
                    particle_size_log10_step=0.003)

    print 'Cotera+01 x 15, amorphous'
    params['Cotera+01 x 15, amorphous'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous',
                    a_min=1e-3,
                    particle_size_log10_step=0.003,
                    scaling_in_particle_number=0.000296296, 
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_size=15.)

    print 'Cotera+01 x 15, amorphous (J98,400C)'
    params['Cotera+01 x 15, amorphous (J98,400C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,400C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003,
                    scaling_in_particle_number=0.000296296, 
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_size=15.)

    print 'Cotera+01 x 15, amorphous (J98,600C)'
    params['Cotera+01 x 15, amorphous (J98,600C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,600C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003,
                    scaling_in_particle_number=0.000296296, 
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_size=15.)

    print 'Cotera+01 x 15, amorphous (J98,800C)'
    params['Cotera+01 x 15, amorphous (J98,800C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,800C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003,
                    scaling_in_particle_number=0.000296296, 
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_size=15.)

    print 'Cotera+01 x 15, amorphous (J98,1000C)'
    params['Cotera+01 x 15, amorphous (J98,1000C)'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='amorphous (J98,1000C)',
                    a_min=1e-3,
                    particle_size_log10_step=0.003,
                    scaling_in_particle_number=0.000296296, 
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_size=15.)

    print 'Cotera+01 x 15, graphite'
    params['Cotera+01 x 15, graphite'][wl] = W.Set_particle_params_Wood02(
                    wavelength_in_micron=wl,
                    dust_model='Cotera',
                    carbon_dust='graphite',
                    a_min=1e-3,
                    particle_size_log10_step=0.000296296,
                                              # do not forget to change for
                                              # different scaling in size
                    scaling_in_particle_number=1e-3,
                    scaling_in_size=15.)

    #print 'Wood+02 Model 1, amorphous'
    #params['Wood+02 Model 1, amorphous'][wl]   = W.Set_particle_params_Wood02(
    #                wavelength_in_micron=wl,
    #                dust_model='Model 1',
    #                carbon_dust='amorphous',
    #                a_min=1e-4,
    #                particle_size_log10_step=0.003)

    #print 'Wood+02 Model 1, graphite'
    #params['Wood+02 Model 1, graphite'][wl]   = W.Set_particle_params_Wood02(
    #                wavelength_in_micron=wl,
    #                dust_model='Model 1',
    #                carbon_dust='graphite',
    #                a_min=1e-3,
    #                particle_size_log10_step=0.003)

    #--
    for dust in dust_list:

      dust_library[dust][wl]['S11']=shrink_array(params[dust][wl].S11)
      dust_library[dust][wl]['S12']=shrink_array(params[dust][wl].S12)
      dust_library[dust][wl]['S33']=shrink_array(params[dust][wl].S33)
      dust_library[dust][wl]['S34']=shrink_array(params[dust][wl].S34)
      dust_library[dust][wl]['kappa_ext (cm2 g-1)']=params[dust][wl].kappa_ext
      dust_library[dust][wl]['kappa_sca (cm2 g-1)']=params[dust][wl].kappa_sca
      dust_library[dust][wl]['g']                  =params[dust][wl].g
      dust_library[dust][wl]['albedo']             =params[dust][wl].albedo
      dust_library[dust][wl]['dust to H mass ratio'] = params[dust][wl].dust_to_H_mass_ratio

      pf = array_for_probability_function_of_cospsi(params[dust][wl].S11)
      dust_library[dust][wl]['probability function for scattering'] = pf

      print "kappa_ext (%s)     = %.2e cm2 g-1" % (dust,params[dust][wl].kappa_ext)


  # add comments
  dust_library['Comments']=raw_input("Add comments:\n")

  # save the dust parameters
  pickle.dump(dust_library,d_save)

  d_save.close()

