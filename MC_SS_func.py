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

  ### Note: This code is adapted for use with the  Sprout GUI.           ###
  ###       See the revision history for 4/23/2013 about what I did. ###

  The original code is MC.py (2011-09-02  Hyosun Kim)
  This consists of the following module & class:

    carray (module): 
    Used to store the numpy array to c-type array.
    Used in the class below and not externally called.

    run_MC_simulation (class):
    Called for main part of Mote-Carlo scattering simulations.

    Shakura_Sunyaev (class):
    Set parameters for the Shakura and Sunyaev disk

  ### History of Revision ###

  9/15/2011  ... Created: 
                the density distribution is replaced as the Shakura & Sunyaev disk.
  9/19/2011  ... The new probability function is inserted.
  2/21/2012  ... Comments are added.
  3/3/2012   ... revised for the new C code with the Shakura & Sunyaev function.
  4/3/2012   ... Save_load_KMH94_Cotera01_Wood02 is replaced as Save_load_KMH94_Cotera01_Wood02_wiscomb
  4/4/2012   ... bugs for normalizing the flux are fixed.
  4/?/2012   ... the image for the intensity-weighted last scattering angle is now provided.
  4/11/2012  ... IGRIDX and IGRIDY are unified as IGRID.
  10/3/2012  ... tidied just a little
  10/18/2012 ... one more iteration is added to the main loop in 
                 MC_simulations.run_simulation to make MAXSCA accurate
                 (removed on 10/19/2012 for the revision described below)
  10/19/2012 ... The photons without scattering are not stacked to the detector anymore.
                 The photon flux is now analytically calculated and saved as
                   MC_simulations.photons_from_the_star (1-D array for 11 detectors).

                 The optical thickness to the star is also saved as
                   MC_simulations.tau_to_the_star.
                 Now the images are nomalized to the incident number of photons for each
                 direction. If you would like to normalize the image to the stellar flux,
                 you have to multiply exp(tau) to correct the extinction toward the star.
  10/25/2012 ... The parameter for colorbar() is saved as self.colorbar_params.
  11/5/2012  ... The classes and library for dust grains are updated.
  12/14/2012 ... self.dpath is defined in the class MC_simulation.

  4/23/2013  ... The module "get_2D_density_distribution" is copied from Shakura_Sunyaev_func.py
                 Jennifer has. This is slightly different from my original version, as she
                 waht to save the density distribution calculated here.

  ### bugs to be fixed ###
  1/8/2013   ... The optical thickness is not correctly scaled with a combination of tau and specific angle
                 (tau=1 for theta=25 or 30 deg.),
                 if beta=1.3, h_50AU=25 AU, alpha < 0.5.

"""

from pylab import *
import ccodes as D
from Save_load_dust_parameters import Load_dust_data
#from particle import dust_parameters        # for power-law distribution
import matplotlib.font_manager

#################--------------------------------------------------------------------------------------
### Constants ###--------------------------------------------------------------------------------------
#################--------------------------------------------------------------------------------------

D.init_genrand(0)                           # initialize the random number generator

# from Shakura_Sunyaev_func.py

AU      = 149597870.7e3 # m
R_solar = 0.0046491     # AU
M_solar = 1.98892e33    # g

legend_fontsize=matplotlib.font_manager.FontProperties(size=10)

###########################----------------------------------------------------------------------------
### Classes and Modules ###----------------------------------------------------------------------------
###########################----------------------------------------------------------------------------

def carray(parr):
  """
    Used to store the numpy array to c-type array. More speficically, used to store
    the following to the C code:-
     - elements for the Muller matrix      (dust.S11,dust.S12,dust.S33,dust.S34)
     - probability function for scattering (dust.probability_func)

  """
  parr = array(parr)
  carr = D.new_darr(parr.size)
  for i in xrange(parr.size): D.darr_setitem(carr, i, parr[i])
  return carr

###----------------------------------------------------------------------------------------------------
class MC_simulation:
  """
    Execute Monte-Carlo scattering simulations. This consists of the following
    two modules:-

    __init__                    : set the parameters to run the simulation code.
    run_simulation              : the main body of the simulation code  
    get_photon_onto_detector    : as described. Called from run_simulation. Not extermally called.

    The program should be used, e.g., as follows:

      params = MC_simulation(
          dust,wave,disk,rho_0,disk_mass,alpha,beta,h_0,rstar,
          rmin,rmax,kappa,tau_crit,tmax,
          ntot, incl, igrid, rang, scat, dpath,
          h_for_arbitrary_r_in_AU,r_for_h_in_AU,
          rho_for_arbitrary_r,r_for_rho_in_AU,
          tau_for_arbitrary_theta,theta_for_arbitrary_tau
          )

      params.run_simulation()

    Revision:

    10/8/2012
    I have changed the order of the input parameters below. I should talk to Jennifer.

  """
  def __init__(self,dust_model,wavelength,              # Input parameres for dust. See
                                                        # Save_load_KMH94_Cotera01_Wood02_wiscomb.py for details

                    disk_model,                         # Input parameres for the disk. See Shakura_Sunyaev_func.py
                    alpha,beta,                         # for details
                    Rmin_in_AU,Rmax_in_AU, 

                    NTOT,                               # Total number of photons
                    IGRID,                              # Pixel number for the final image
                    INCL,                               # Inclination for the specific veiwing angle
                    IRANGE,                             # Range for the angle for the above inclination angle
                    dpath,                              # step for integrating tau for radiative transfer (AU)

                    ### the parameters with default values ###

                    # (scaling the density of the disk)

                    rho_0                     =  1.,    # g cm-3
                    disk_mass                 = -1.,    # solar mass
                    rho_for_arbitrary_r       = -1.,    # If set, the value is used to scale rho_0
                    r_for_rho_in_AU           =  1.,    # in unit of AU
                    tau_for_arbitrary_theta   = -1.,    # If set, the value is used to scale rho_0
                    theta_for_arbitrary_tau   = 20.,    # deg.

                    # (scale height of the disk)

                    h_0                       =  0.017, # in unit of stellar radius
                    h_for_arbitrary_r_in_AU   = -1.,    # If set, the value is used to scale h_0
                    r_for_h_in_AU             =  1.,    # in unit of AU

                    # (for the other model parameters)
                    
                    kappa_ext                 = -1.,    # As described. Given by the dust model
                    R_star_in_solar_radius    =  1.2,

                    # (for simulations)

                    half_image_size_in_AU     = -1,     # if default, it will match the outer radius of the disk

                    MAXSCA                    = 10,     # Maximum number of scattering

                    tau_critical              = 0.001,  # Critical optical thickness for radial direction.
                                                        # The photons are ejected from the star at
                                                        # angles for which the optical thickness larger than
                                                        # tau_critical.

                    Theta_max                 = -1,     # The critical angle from the midplane.
                                                        # The photons are ejected from the star within
                                                        # this angle. If -1 is set, it should be
                                                        # tau_critical.

                    seed_for_rn_generator     = 0,      # seed of the random number genetator
                                                        # (added in Oct 2012)

                    show_2D_density_distribution =  False

               ):


    #print "%-25s" % ("dust_model"),dust_model
    #print "%-25s" % ("wavelength"),wavelength
    #print "%-25s" % ("disk_model"),disk_model
    #print "%-25s" % ("rho_0"),rho_0
    #print "%-25s" % ("disk_mass"),disk_mass
    #print "%-25s" % ("alpha"),alpha
    #print "%-25s" % ("beta"),beta
    #print "%-25s" % ("h_0"),h_0
    #print "%-25s" % ("R_star_in_solar_radius"),R_star_in_solar_radius
    #print "%-25s" % ("Rmin_in_AU"),Rmin_in_AU
    #print "%-25s" % ("Rmax_in_AU"),Rmax_in_AU
    #print "%-25s" % ("kappa_ext"),kappa_ext # RBIN is removed added on 3/12/2012 
    #print "%-25s" % ("tau_critical"),tau_critical
    #print "%-25s" % ("Theta_max"),Theta_max
    #print "%-25s" % ("NTOT"),NTOT
    #print "%-25s" % ("INCL"),INCL
    #print "%-25s" % ("IGRID"),IGRID
    #print "%-25s" % ("IRANGE"),IRANGE
    #print "%-25s" % ("MAXSCA"),MAXSCA
    #
    #                ### the parameters below are newly added on 3/12/2012 ###
    #
    #print "%-25s" % ("dpath"),dpath
    #print "%-25s" % ("h_for_arbitrary_r_in_AU"),h_for_arbitrary_r_in_AU
    #print "%-25s" % ("r_for_h_in_AU"),r_for_h_in_AU
    #print "%-25s" % ("rho_for_arbitrary_r"),rho_for_arbitrary_r
    #print "%-25s" % ("r_for_rho_in_AU"),r_for_rho_in_AU
    #print "%-25s" % ("tau_for_arbitrary_theta"),tau_for_arbitrary_theta
    #print "%-25s" % ("theta_for_arbitrary_tau"),theta_for_arbitrary_tau


    self.dust_model=dust_model
    self.NTOT=int(NTOT)
    self.IGRID=IGRID
    self.INCL=INCL
    self.IRANGE=IRANGE
    self.MAXSCA=MAXSCA
    self.tau_critical=tau_critical
    self.Rmax_in_AU=Rmax_in_AU
    # added in 12/14/2012
    self.dpath=dpath

    # added in Oct 2012
    if half_image_size_in_AU > 0:
      self.half_image_size_in_AU=half_image_size_in_AU
    else:
      self.half_image_size_in_AU=Rmax_in_AU

    # Commented out in Oct 2012

    # self.Theta_max=Theta_max


    ###########################################################################
    #  set physical parameters
    ###########################################################################

    ###  dust parameters ###
    #print "preparing for dust parameters..."
    #dust = dust_parameters(wavelength        = wavelength,
    #                       amax              = max_dust,
    #                       amin              = min_dust,
    #                       amin_core         = min_core,
    #                       power_law         = power_law,
    #                       thickness_of_coat = thickness_of_coat,
    #                       number_of_bins    = number_of_bins_for_dust_size,
    #                       core_number_fractions   = core_number_fractions,
    #                       core_refractive_indexes = core_refractive_indexes,
    #                       coat_refractive_index   = coat_refractive_index)
    #dust=Load_dust_data(dust='Cotera+01, amorphous',wavelength_in_micron=1.65)
    self.dust=Load_dust_data(dust_model,wavelength)

    D.cvar.albedo = self.dust.albedo
    D.cvar.S11              = carray(self.dust.S11)
    D.cvar.S12              = carray(self.dust.S12)
    D.cvar.S33              = carray(self.dust.S33)
    D.cvar.S34              = carray(self.dust.S34)
    D.cvar.probability_func = carray(self.dust.probability_func)

    # We do not use the HG function for scattering anymore,
    # so it should not be necessary unseed_for_rn_generatorless we revise the code (Oct 2012)

    # D.cvar.g = self.dust.g 

    ###  disk parameters ###
    #print "preparing for disk parameters..."

    if kappa_ext <= 0:
      kappa_ext=self.dust.kappa_ext

    self.SS_disk  = Shakura_Sunyaev(
                    parameters = disk_model,
                    rho_0      = rho_0,
                    disk_mass  = disk_mass,
                    rho_for_arbitrary_r       = rho_for_arbitrary_r,
                    r_for_rho_in_AU           = r_for_rho_in_AU,
                    tau_for_arbitrary_theta   = tau_for_arbitrary_theta,
                    theta_for_arbitrary_tau   = theta_for_arbitrary_tau,
                    h_0                       = h_0,
                    h_for_arbitrary_r_in_AU   = h_for_arbitrary_r_in_AU,
                    r_for_h_in_AU             = r_for_h_in_AU,
                    alpha                     = alpha,
                    beta                      = beta,
                    Rmin_in_AU                = Rmin_in_AU,
                    Rmax_in_AU                = Rmax_in_AU,
                    kappa_ext                 = kappa_ext,
                    R_star_in_solar_radius    = R_star_in_solar_radius,
                    show_2D_density_distribution =  show_2D_density_distribution,
                    dpath                        =  dpath,     # step for integrating tau for radiative transfer (AU)


                    tau_critical              = tau_critical,
                    Theta_max                 = Theta_max      # maximum angle from the midplane to eject the photons

                    )

    if disk_mass <= 0:
      self.SS_disk.get_total_mass_analytically()

    # Initialize the random number generator. This has been moved here in Oct 2012.
    D.init_genrand(seed_for_rn_generator)                           

  #----------------------------------------------------------------
  def run_simulation(self):
    """
      The part below is now split for the class to be able to use to get lines for tau=0.5,1,2.

        3/13/2012 Hiro Takami
    """

    ###########################################################################
    #  initialization of detectors 
    ###########################################################################

    self.number_of_photons_killed = 0
    self.NDIFS=zeros(11) # number of scattering
    self.NPHOT=zeros(11) # number of photons onto each detector
    self.image_I=zeros((11,self.IGRID,self.IGRID)) # Stokes parameter I @ each pixel
    self.image_Q=zeros((11,self.IGRID,self.IGRID)) # Stokes parameter Q @ each pixel
    self.image_U=zeros((11,self.IGRID,self.IGRID)) # Stokes parameter U @ each pixel
    self.image_V=zeros((11,self.IGRID,self.IGRID)) # Stokes parameter V @ each pixel

    self.image_average_number_of_scattering      =zeros((11,self.IGRID,self.IGRID)) # as described
    self.image_averaged_angle_for_last_scattering=zeros((11,self.IGRID,self.IGRID)) # as described

    INCLMIN = cos(radians(min(90., self.INCL+self.IRANGE*0.5))) # inclination <= 90 deg
    INCLMAX = cos(radians(max( 0., self.INCL-self.IRANGE*0.5))) # inclination >=  0 deg

    ###########################################################################
    #  main part of the simulation
    ###########################################################################

    #print "starting the scattering simulation"

    newIQUV=D.new_darr(4);   newXYZ=D.new_darr(3);   newUVW=D.new_darr(3)
    IQUV=D.new_darr(4)   ;   XYZ=D.new_darr(3)   ;   UVW=D.new_darr(3)

    # Those lines were in "set physical parameters", but I have moved here as
    # the new algorism needs the definition for the above parameters.
    # 
    #   3/13/2012 Hiro Takami
    #
    # This part has been moved to Shakura_Sunyaev. (Oct 2012)
    #
    #if self.Theta_max <= 0:
    #  Theta_max=self.SS_disk.get_Theta_max(self.tau_critical,IQUV,XYZ,UVW,newXYZ)
    #else:
    #  Theta_max  = self.Theta_max
    #D.cvar.Tdisk = radians(Theta_max)
    #-

    for I_PHOTON in xrange(self.NTOT):
      D.Initialize(IQUV, XYZ, UVW)

      # Algorism before the 1st scattering.
      # The photons without scattering are not stacked onto the detector anymore
      # This part is pulled put from the loop below on 10/19/2012.

      photon_is_out_of_nebula = D.NewPosition(newXYZ, XYZ, UVW)

      if photon_is_out_of_nebula:
        continue

      for i in xrange(3):
        D.darr_setitem(XYZ, i, D.darr_getitem(newXYZ, i))

      # Algorism for and after the 1st scattering.
      # In Sunny's algorism, she calculated the next position first,
      # then determine the new photon vector and Stokes parameters.
      # I have swapped the order.
      #
      # 10/19/2012
      
      for N_scatter in range(1,self.MAXSCA+1):

        D.StokesMuller(newIQUV, newUVW, IQUV, UVW)        
        for i in xrange(4):
          D.darr_setitem(IQUV, i, D.darr_getitem(newIQUV, i))
        for i in xrange(3):
          D.darr_setitem(UVW, i, D.darr_getitem(newUVW, i))

        photon_is_out_of_nebula = D.NewPosition(newXYZ, XYZ, UVW)
        if photon_is_out_of_nebula:
          self.get_photon_onto_detector(IQUV, XYZ, UVW, N_scatter,
                                        self.IGRID,INCLMIN,INCLMAX)

                                        # self.Rmax_in_AU is removed from the input parameter.
                                        # The module now uses self.half_image_size_in_AU. (Oct 2012)
          break
        for i in xrange(3):
          D.darr_setitem(XYZ, i, D.darr_getitem(newXYZ, i))

      else: # if the number of scattering excees MAXSCA
        self.number_of_photons_killed += 1

    D.delete_darr(IQUV)   ;  D.delete_darr(XYZ)   ;  D.delete_darr(UVW)
    D.delete_darr(newIQUV);  D.delete_darr(newXYZ);  D.delete_darr(newUVW)

    self.image_average_number_of_scattering       /= self.image_I
    self.image_averaged_angle_for_last_scattering *= 180./pi/self.image_I

    # now IQUV is normalized to the total number of photon.
    # Note that the this does not provide accurate normalizeation to the stellar photon:
    # we have to include the effect of extinction. The optical thickness to the star
    # for each viewing angle is now separately calculated and saved.
    scaling_factor                    = self.NTOT / sin(radians(self.SS_disk.Theta_max)) / 10
    scaling_factor_for_specific_angle = self.NTOT / sin(radians(self.SS_disk.Theta_max)) * (INCLMAX-INCLMIN)

    self.image_I[:10,:,:] /= scaling_factor
    self.image_Q[:10,:,:] /= scaling_factor
    self.image_U[:10,:,:] /= scaling_factor
    self.image_V[:10,:,:] /= scaling_factor
    self.image_I[10,:,:]  /= scaling_factor_for_specific_angle
    self.image_Q[10,:,:]  /= scaling_factor_for_specific_angle
    self.image_U[10,:,:]  /= scaling_factor_for_specific_angle
    self.image_V[10,:,:]  /= scaling_factor_for_specific_angle

    # Calcuate and save the optical thickness to the star.

    self.tau_to_the_star      =[]
    for theta in (180/pi*arccos(0.95-arange(10)*0.1)):
      tau=self.SS_disk.get_tau_for_theta_equal_to(90-theta)  # "90-" is added, as theta in the array
                                                             # is the viewing angle, while what module
                                                             # wants is the angle from the midplane.
      self.tau_to_the_star.append(tau)

    tau=self.SS_disk.get_tau_for_theta_equal_to(90-self.INCL)
    self.tau_to_the_star.append(tau)

  #----------------------------------------------------------------
  def get_photon_onto_detector(self, IQUV, XYZ, UVW, N_scatter,
                               IGRID,INCLMIN,INCLMAX):
    """
      (Comments added on 10/1/2012)
      Now we have the following detectors

        self.image_I
        self.image_Q
        self.image_U
        self.image_V
        self.image_average_number_of_scattering
        self.image_averaged_angle_for_last_scattering[CAPT,IMAP,JMAP]

      The image size is determined using Rmax.
      > Now it is determined by self.half_image_size_in_AU.
        So Rmax is removed from the input parameter,
        and that in the remaining is replaced as self.half_image_size_in_AU.
        (Oct 2012)
    """

    W = D.darr_getitem(UVW, 2)
    if (W < 0.):
      D.darr_setitem(XYZ, 1, -D.darr_getitem(XYZ, 1))
      D.darr_setitem(XYZ, 2, -D.darr_getitem(XYZ, 2)) 
      D.darr_setitem(UVW, 1, -D.darr_getitem(UVW, 1))
      W=-W;                 D.darr_setitem(UVW, 2, W)

    photon_IQUV = array([D.darr_getitem(IQUV, i) for i in xrange(4)])
    photon_IQUV[2:] *= -1

    XYZ_modified=D.new_darr(3);    D.VectorRotate(XYZ_modified, XYZ, UVW)
    YPRIM = D.darr_getitem(XYZ_modified, 1)
    ZPRIM = D.darr_getitem(XYZ_modified, 2)
    D.delete_darr(XYZ_modified)

    IMAP=int((YPRIM+self.half_image_size_in_AU)/(2.*self.half_image_size_in_AU)*IGRID);  IMAP=min(IGRID-1, max(0, IMAP))
    JMAP=int((ZPRIM+self.half_image_size_in_AU)/(2.*self.half_image_size_in_AU)*IGRID);  JMAP=min(IGRID-1, max(0, JMAP))

    CAPT = min(9, int( (1-W)*10 ))
    self.NPHOT[CAPT] += 1
    self.NDIFS[CAPT] += N_scatter
    self.image_I[CAPT,IMAP,JMAP] += photon_IQUV[0]
    self.image_Q[CAPT,IMAP,JMAP] += photon_IQUV[1]
    self.image_U[CAPT,IMAP,JMAP] += photon_IQUV[2]
    self.image_V[CAPT,IMAP,JMAP] += photon_IQUV[3]
    self.image_average_number_of_scattering[CAPT,IMAP,JMAP] += photon_IQUV[0]*N_scatter
    self.image_averaged_angle_for_last_scattering[CAPT,IMAP,JMAP] += photon_IQUV[0]*D.cvar.PSI

    if W >= INCLMIN and W < INCLMAX:
      self.NPHOT[10] += 1
      self.NDIFS[10] += N_scatter
      self.image_I[10,IMAP,JMAP] += photon_IQUV[0]
      self.image_Q[10,IMAP,JMAP] += photon_IQUV[1]
      self.image_U[10,IMAP,JMAP] += photon_IQUV[2]
      self.image_V[10,IMAP,JMAP] += photon_IQUV[3]
      self.image_average_number_of_scattering[10,IMAP,JMAP] += photon_IQUV[0]*N_scatter
      self.image_averaged_angle_for_last_scattering[10,IMAP,JMAP] += photon_IQUV[0]*D.cvar.PSI

###----------------------------------------------------------------------------------------------------

class Shakura_Sunyaev:
  """

    Prepare for density distribution based on the Sakura & Sunyaev equation.
    Originally in Shakura_Sunyaev_func.py but copied here to share
    the c codes (ccode/D) for some calculations such as the distance from the star
    at a certain opacity.

    The input parameters are a sort of duplicated to determine the density
    distribution of the disk. The density distribution is determined with
    the following priorities.

    (1)
    If 'Cotera+01, 1.65 um' or 'Wood+02, 1.65 um' is set as the input parameter
    "parameters":
    > The class sets the parameters ignoring all the other input
      parameters but RBIN.

    (2)
    If the input parameter "tau_for_arbitrary_theta" is set:
    > The density distribution is scaled using the input parameters
      "tau_for_arbitrary_theta" and "theta_for_arbitrary_tau", ignoring
      the input parameters "rho_0", "disk_mass", "rho_for_arbitrary_r" and
      "r_for_rho_in_AU".

    (3)
    If the input parameter "disk_mass" is set:
    > The density distribution is scaled using this parameter,
      ignoring the input parameters "rho_0", "rho_for_arbitrary_r" and
      "r_for_rho_in_AU".

    (4)
    If the input parameter "rho_for_arbitrary_r" is set:
    > The density distribution is scaled using the input parameters
      "rho_for_arbitrary_r" and "r_for_rho_in_AU", ignoring the input
      parameters "rho_0". (The input parameter "R_star" is still used,
      but it does not affect the resultant scale height.)

    (5)
    If the input parameter "h_for_arbitrary_r_in_AU" is set:
    > The scale height is determined based on two input parameters
      "h_for_arbitrary_r_in_AU" and "r_for_h_in_AU", ignoring the input
      parameters "h_0". (The input parameter "R_star" is still used,
      but it does not affect the resultant scale height.)


    ### Input parameters ###
      parameters                ... disk parameters, which should be either
                                    'Cotera+01, 1.65 um','Wood+02, 1.65 um',
                                    or anything else for manual input
                                                        (default =  '')
      rho_0                     ... density at the stellar radius (g cm-3) 
                                                        (default =  1.)
      disk_mass                 ... in unit of M_solar.
                                    If -1 is set the disk mass is determined
                                    from the other parameters. Otherwise, the
                                    density disbturbtion in the entire disk
                                    is scaled based on the total mass defined.
                                                        (default = -1)
      alpha                     ... power index         (default =  2.25)
      beta                      ... power index         (default =  1.25)
      h_0                       ... scale height at the stellar radius, in unit of stellar radius
                                                        (default =  0.017)
      R_star_in_solar_radius    ... as described        (default =  1.2)
      Rmin_in_AU                ... as described. The class will provide an error message
                                      if minus.         (default = -1)
      Rmax_in_AU                ... as described. The class will provide an error message
                                      if minus.         (default = -1)
      kappa_ext                 ... as described, in unit of cm2 g-1.
                                                        (default = -1)
                                    You must specify if you do not set the specific disk model
                                    (either 'Cotera+01, 1.65 um' or 'Wood+02, 1.65 um' for
                                    "parameters").
                                    If the specific model is set, you may leave it as the default
                                    value, and that in the paper will be used.

      h_for_arbitrary_r_in_AU   ... If set, the value is used to scale h_0.
                                    Ignored if -1.      (default = -1)
      r_for_h_in_AU             ... radius for the scale height defined above.
                                    In unit of AU.      (default = 1.)
      rho_for_arbitrary_r       ... If set, the value is used to scale rho_0 ().
                                    Ignored if -1.      (default = -1)
      r_for_rho_in_AU           ... radius for density defined above.
                                    In unit of AU.      (default = 1.)
      tau_for_arbitrary_theta   ... If set, the entire density distribution is scaled
                                    using this optical thickness integrated over r.
                                    Ignored if -1.      (default = -1)
      theta_for_arbitrary_tau   ... theta for tau_for_arbitrary_theta
                                                        (default = 20)

      == the parameter below is not used for simulations anymore, but I still keep it ==
      == just in case for the analysis of the disk itself. (4/11/2012)                ==

      RBIN                      ... as described        (default = 200)


    ### Output parameters ###

      self.rho_0                ... Same as the input parameters.
      self.RBIN                 ...
      self.Rmin_in_AU           ... 
      self.Rmax_in_AU           ... 
      self.alpha                ... 
      self.beta                 ... 
      self.kappa_ext            ... 

      self.R_star,self.h_0      ... Same as the input parameters but in AU

      self.tau_0                ... tau per AU at the same position as rho_0

      self.total_mass_analytic  ... total mass of the disk determined analytically.
                                    It should provide a (sligtly) larger value
                                    than the disk limited in the sphere R_max, as the value is
                                    determined by intergrating over the whole z.

      == the parameters below are not used for simulations anymore, but I still keep it ==
      == just in case for the analysis of the disk itself. (4/11/2012)                  ==

      self.dR                   ... Step for binning (AU)


      self.density_distribution ... array for density. The size is (self.RBIN+1,91).
                                    Set when self.get_density_distribution is executed.

      self.opacity_distribution ... array for opacity. The size is (self.RBIN+1,91).
                                    Set when self.opacity_density_distribution is executed.

      self.op_r                 ... optical thickness integrated over r as a function of theta.
                                    Set if self.get_opacity_as_a_function_of_theta is executed.

      self.total_mass_numerical ... total mass of the disk determined numerically.
                                    This does not work well for some disk we use, as the grid
                                    for the angle is too coarse. So self.total_mass_analytic
                                    is alternatively used elsewhere.

    ### Modules ###

    self.set_params_Cotera01_H()    ... Set parameters for Cotera+01.
                                        Kappa_ext is for Cotera+01, amorphous, 1.61 um
                                        This also makes the following parameters:-
                                           self.density_distribution
                                           self.mass_for_each_cell
                                           self.total_mass_analytic

    self.set_params_Wood02_H()      ... Set parameters for Wood+02.
                                        Kappa_ext is for Cotera+01, amorphous, 1.61 um
                                        This also makes the following parameters:-
                                           self.density_distribution
                                           self.mass_for_each_cell
                                           self.total_mass_analytic

    self.h(r)                       ... Give the scale hight for arbitrary r (in AU)
                                        based on the Shakura & Sunyaev equation.

    self.rho(r,z)                   ... Give the density for arbitrary r, z (in AU)
                                        based on the Shakura & Sunyaev equation.

    self.get_total_mass_analytically()     ... Get total mass analytically. This makes:
                                           self.total_mass_analytic

    self.scaling_density_and_mass_with_disk_mass(disk_mass)
                                    ... as described. The disk mass should be in unit of
                                        solar mass. This makes the following parameters
                                        if not created:-
                                           self.density_distribution
                                           self.mass_for_each_cell
                                           self.total_mass_analytic
                                        
    self.show_2D_density_distribution()
                                    ... As described.

    self.print_params()             ... Print parameters as txt format.

    get_Theta_max                   ... get Theta_max from tau_critical. See below for details of parameters.
    get_lines_for_tau_from_star     ... as described.

  """

  # Ihe order of the input parameters has changed. It should not affect the other part
  # of the code. (10/1/2012)

  def __init__(self,
                parameters                   = '',
                                              # 'Cotera+01, 1.65 um',
                                              # 'Wood+02, 1.65 um',
                                              #  or anything else for manual input

                ### for scaling the density ###
                rho_0                        =   1.,
                disk_mass                    =  -1,
                rho_for_arbitrary_r          =  -1.,   # Added on 9/27/2011. If set,
                                                       # the value is used to scale rho_0
                r_for_rho_in_AU              =   1.,   # in unit of AU
                tau_for_arbitrary_theta      =  -1.,   # change the default value to -1 at some point,
                                                       #   and this scaling will not work as default.
                theta_for_arbitrary_tau      =  20.,   # degree

                ### for scale height ###
                h_0                          =  0.017, # in unit of stellar radius
                h_for_arbitrary_r_in_AU      =  -1.,   # Added on 9/27/2011. If set,
                                                       # the value is used to scale h_0
                r_for_h_in_AU                =   1.,   # in unit of AU

                ### other disk parameters ###
                alpha                        =  2.25,
                beta                         =  1.25,
                Rmin_in_AU                   =  -1,
                Rmax_in_AU                   =  -1,
                kappa_ext                    =  -1.,
                R_star_in_solar_radius       =  1.2,

                ### model parameters ###
                RBIN                         =  200,   # grid number to show the density profile
                dpath                        =  1,     # step for integrating tau for radiative transfer (AU)
                show_2D_density_distribution =  False,

                tau_critical                 =  0.001,
                Theta_max                    =  0      # maximum angle from the midplane to eject the photons
                                                       # from the star.
               ):
    """
      The above parameters will be set in this module.
      If the disk mass is not set, it will be calculated from these parameters.
      Otherwise, rho_0 is scaled based on the disk mass provided.
    """
    self.rho_0        = rho_0
    self.RBIN         = RBIN

    if parameters == 'Cotera+01, 1.65 um':
      self.Rmin_in_AU = Rmin_in_AU  # the values will be the same as the paper if these are -1
      self.Rmax_in_AU = Rmax_in_AU
      self.kappa_ext  = kappa_ext

      self.set_params_Cotera01_H()

    elif parameters == 'Wood+02, 1.65 um':
      self.Rmin_in_AU = Rmin_in_AU  # the values will be the same as the paper if these are -1
      self.Rmax_in_AU = Rmax_in_AU
      self.kappa_ext  = kappa_ext

      self.set_params_Wood02_H()

    else:
      if Rmin_in_AU <= 0:
        print " ### Error!! ###"
        print "Set a positive value for Rmin_in_AU (inner boundary of the disk)."
        return
      if Rmax_in_AU <= 0:
        print " ### Error!! ###"
        print "Set a positive value for Rmax_in_AU (outer boundary of the disk)."
        return
      if kappa_ext  <= 0:
        print " ### Error!! ###"
        print "Set a positive value for kappa_ext."
        return
      self.alpha      = alpha
      self.beta       = beta
      self.R_star     = R_star_in_solar_radius * R_solar # AU
      self.h_0        = h_0 *  self.R_star               # AU
      self.Rmin_in_AU = Rmin_in_AU                       # AU
      self.Rmax_in_AU = Rmax_in_AU                       # AU
      self.dR         = Rmax_in_AU*1./RBIN               # AU
      self.kappa_ext  = kappa_ext

      # scaling h_0 if the h_for_arbitrary_r is set #
      if h_for_arbitrary_r_in_AU > 0.:
        self.h_0 = h_for_arbitrary_r_in_AU / self.h(r_for_h_in_AU) * self.h(self.R_star)

      # scaling rho_0 if the rho_for_arbitrary_r is set #
      if rho_for_arbitrary_r > 0.:
        self.rho_0 = rho_for_arbitrary_r / self.rho(r_for_rho_in_AU,0) * self.rho(self.R_star,0)

    if disk_mass > 0:
      # scaling rho_0 if the disk mass is given
      self.scaling_density_and_mass_with_disk_mass(disk_mass)
    else:
      self.get_total_mass_analytically()

    self.tau_0=self.rho_0*AU*1e2*self.kappa_ext

    # put parameters to the C code
    D.cvar.alpha                 = self.alpha
    D.cvar.beta                  = self.beta
    D.cvar.h_0                   = self.h_0
    D.cvar.stellar_radius_in_AU  = self.R_star
    D.cvar.Shakura_Sunyaev_tau_0 = self.tau_0
    D.cvar.Rmin                  = Rmin_in_AU
    D.cvar.Rmax                  = Rmax_in_AU
    D.cvar.dpath                 = dpath   

    # scaling the density distribution if tau(theta) is given
    if tau_for_arbitrary_theta >= 0:
      tau_now     = self.get_tau_for_theta_equal_to(theta_for_arbitrary_tau)
      self.rho_0 *= tau_for_arbitrary_theta/tau_now
      self.tau_0  =self.rho_0*AU*1e2*self.kappa_ext
      D.cvar.Shakura_Sunyaev_tau_0 = self.tau_0
      self.get_total_mass_analytically()

      ###
      # try to iterate once again as the above calculation is not always accurate
      # (see the bug report in the beginning of the file).
      # Unfortunately, it does not work well.
      ####
      #tau_now     = self.get_tau_for_theta_equal_to(theta_for_arbitrary_tau)
      #self.rho_0 *= tau_for_arbitrary_theta/tau_now
      #self.tau_0  =self.rho_0*AU*1e2*self.kappa_ext
      #D.cvar.Shakura_Sunyaev_tau_0 = self.tau_0
      #self.get_total_mass_analytically()
      #
      #tau_now     = self.get_tau_for_theta_equal_to(theta_for_arbitrary_tau)
      #print tau_now,tau_now,tau_now,tau_for_arbitrary_theta

    # The lines below have been moved from MC_SS in Oct 2012.
    # We do not always need Theta_max when we use this class
    # (in particular if we do not run the MC simulations),
    # but it makes MC_SS look simpler.
 
    if Theta_max <= 0:
      self.Theta_max=self.get_Theta_max(tau_critical)
    else:
      self.Theta_max=Theta_max
    D.cvar.Tdisk = radians(self.Theta_max)
    
    #-

    if show_2D_density_distribution:
      self.show_2D_density_distribution()

    """
      ### Copied from Jennifer's old Shakura_Sunyaev on 4/13/2013. ###
      ### Copied from Jennifer's old Shakura_Sunyaev on 4/13/2013. ###
      ### Copied from Jennifer's old Shakura_Sunyaev on 4/13/2013. ###

      As described. No new parameters are created.

      ### Input/Output parameters ###
      None. Input parameters are all set in __init__, etc.

      ### Note ###
      You should execute the figure command (and also subplot command if you want)
      before executing the program.

      While self.density_distribution is provided with r-theta coordinate, this require
      distribution in x-y coordinate. Such distribution is calculated here, thus
      we do not need to run any modules (e.g., self.get_density_distribution,
      self.get_opacity_distribution before executing this program.
    """

  #def get_2D_density_distribution(self):

    density_distribution=zeros((self.RBIN+1,self.RBIN+1),float)
    RBIN2=self.RBIN*self.RBIN

    for i_rr in range(self.RBIN+1):
      for i_zz in range(self.RBIN+1):
        if i_rr > 0:
          if (i_rr**2+i_zz**2) <= RBIN2:
            density_distribution[i_zz,i_rr]=self.rho(i_rr*self.dR,i_zz*self.dR)
          else:
            density_distribution[i_zz,i_rr]=NaN

    self.twod_density_distribution=log10(density_distribution)

  #---
  def set_params_Cotera01_H(self):
    """
      The same parameter set as Cotera+01 is provided.
      The default values are also saved for rho_0, alpha, beta and h_0 for scaling.

      We can change the following parameters, if set in the beginning:-

        Rmin_in_AU (inner boundary of the disk)
        Rmax_in_AU (outer boundary of the disk)
        kappa_ext

      If not set, the default values in the paper are used.

    """

    #print "Setting parameters for the Cotera+01 disk..."

    self.alpha          = 2.36667
    self.beta           = 58./45
    self.R_star         = 2.5   * R_solar         # AU
    self.h_0            = 0.011 * self.R_star     # AU
    if self.Rmin_in_AU <= 0:
      self.Rmin_in_AU   = 6     * self.R_star     # AU
    if self.Rmax_in_AU <= 0:
      self.Rmax_in_AU   = 200                     # AU
    if self.kappa_ext  <= 0:
      self.kappa_ext    = 87.9                    # Cotera+01, amorphous, 1.61 um
                                                  # Calculated value using our code is 101.3 for 1.65 um
    disk_mass           = 6.7e-4                  # solar mass

    self.dR        = self.Rmax_in_AU/self.RBIN    # AU


    self.scaling_density_and_mass_with_disk_mass(disk_mass)

    return

  #---
  def set_params_Wood02_H(self):
    """
      The same parameter set as Wood+02 is provided.
      The default values are also saved for rho_0, alpha, beta and h_0 for scaling.

      We can change the following parameters, if set in the beginning:-

        Rmin_in_AU (inner boundary of the disk)
        Rmax_in_AU (outer boundary of the disk)
        kappa_ext

      If not set, the default values in the paper are used.

    """

    #print "Setting parameters for the Wood+02 disk..."

    self.alpha      = 2.25
    self.beta       = 1.25
    self.R_star     = 1.2   * R_solar           # AU
    self.h_0        = 0.017 * self.R_star       # AU
    if self.Rmin_in_AU <= 0:
      self.Rmin_in_AU = 6   * self.R_star       # AU
    if self.Rmax_in_AU <= 0:
      self.Rmax_in_AU = 200                     # AU
    if self.kappa_ext  <= 0:
      self.kappa_ext  = 87.9                    # Cotera+01, amorphous, 1.61 um
                                                # Calculated value using our code is 101.3 for 1.65 um
    #disk_mass       = 3.5e-4                   # solar mass
    self.rho_0       = 1.248e-7                # manually calculated
    self.dR         = self.Rmax_in_AU/self.RBIN # AU

    self.get_total_mass_analytically()           # maybe slightly lower than the original paper as
                                                # we set R_max in the (r,theta,phi) coordinate.

    return

  #---
  def rho(self,r_in_AU,z_in_AU):
    """
      Density provided by Shakura & Sunyaev equation in Cotera+01 and Wood+02.

      ### Input parameters ###
      r_in_AU,z_in_AU

      ### Output parameters ###
      density in g cm-3

    """

    r_in_stellar_radius = r_in_AU/self.R_star
    value = self.rho_0 / r_in_stellar_radius**(self.alpha) * exp (-0.5*(z_in_AU/self.h(r_in_AU))**2)

    return value

  #---
  def h(self,r_in_AU):
    """
      Scale height based on the Shakura & Sunyaev equation.

      ### Input parameter ###
      r_in_AU

      ### Output parameter ###
      h (scale hight) in AU

    """

    r_in_stellar_radius = r_in_AU/self.R_star
    return self.h_0 * (r_in_stellar_radius)**self.beta

  #---
  def get_total_mass_numerically(self):
    """
      The old version with an array has removed. I may add a new version using the algorism
      for radiative transfer.
    """

    return

  #---
  def get_total_mass_analytically(self):
    R_star = self.R_star     * AU * 100
    h_0    = self.h_0        * AU * 100
    R_min  = self.Rmin_in_AU * AU * 100
    R_max  = self.Rmax_in_AU * AU * 100

    const  = (2*pi)**1.5 * self.rho_0 * h_0 / M_solar
    power  = self.beta-self.alpha+2
    self.total_mass_analytic = const * R_star ** (self.alpha-self.beta) / power * (R_max ** power - R_min ** power)

  #---
  def scaling_density_and_mass_with_disk_mass(self,disk_mass):
    """
      This revises the following parameters to make them consistent with the disk mass provided:-

        self.rho_0                : density at the inner edge
        self.total_mass_analytic  : total mass

    """

    self.get_total_mass_analytically()  # get the disk mass for given rho_0

    scaling     = disk_mass*1./self.total_mass_analytic
    self.rho_0 *= scaling
    self.get_total_mass_analytically()  # just to reset the value self.total_mass_analytic

    return

  #---
  def show_2D_density_distribution(self,
                                   show_tau_lines = True,
                                   linewidth      = 3,
                                   scaling        = 1.
                                   ):
    """
      As described.

      (11/12/2012)
      A scaling factor is added for density.
      This is because, for the paper, I decide to show the distribution for dust mass,
      not the total mass of the disk.

      ### Input parameters ###
      show_tau_lines, linewidth ... as described
      scaling                   ... the scaling factor for density

      ### Output parameters ###
      None. 

      ### Note ###
      You should execute the figure command (and also subplot command if you want)
      before executing the program.

      While self.density_distribution is provided with r-theta coordinate, this require
      distribution in x-y coordinate. Such distribution is calculated here, thus
      we do not need to run any modules (e.g., self.get_density_distribution,
      self.get_opacity_distribution before executing this program.
    """

    #self.print_params()

    density_distribution=zeros((self.RBIN+1,self.RBIN+1),float)
    RBIN2=self.RBIN*self.RBIN

    for i_rr in range(self.RBIN+1):
      for i_zz in range(self.RBIN+1):
        if i_rr > 0:
          if (i_rr**2+i_zz**2) <= RBIN2:
            density_distribution[i_zz,i_rr]=self.rho(i_rr*self.dR,i_zz*self.dR)
          else:
            density_distribution[i_zz,i_rr]=NaN

    log_d=log10(density_distribution*scaling)

    imshow(log_d,origin='lower',interpolation='nearest',
           extent=[0,self.Rmax_in_AU,0,self.Rmax_in_AU],
           vmin=log_d[0,self.RBIN]-1,
           vmax=log_d[0,self.RBIN/10],)

    xlabel('r (AU)')
    ylabel('z (AU)')
    title('log rho (g cm-3)')
    self.colorbar_params=colorbar() # added on 10/25/2012

    if show_tau_lines:
      r_array={}
      z_array={}

      for tau in (0.5,1.,2.):
        key="tau = %3.1f" % tau
        theta_array,r_array[key],z_array[key]=self.get_lines_for_tau_from_star(tau,d_theta=0.1)
  
      plot(r_array['tau = 1.0'],z_array['tau = 1.0'],'w-' ,linewidth=linewidth)
      plot(r_array['tau = 0.5'],z_array['tau = 0.5'],'w--',linewidth=linewidth)
      plot(r_array['tau = 2.0'],z_array['tau = 2.0'],'w:' ,linewidth=linewidth)

      axis([0,self.Rmax_in_AU,0,self.Rmax_in_AU])

    show()

  #---
  def print_params(self):
    """
      As described.
    """
    print "rho_0                  : ",self.rho_0
    if self.has_param('total_mass_analytic'):
      print "total_mass (analytic)  : ",self.total_mass_analytic
    print "alpha                  : ",self.alpha
    print "beta                   : ",self.beta
    print "h_0                    : ",self.h_0
    print "R_star                 : ",self.R_star
    print "Rmin (AU)              : ",self.Rmin_in_AU
    print "Rmax (AU)              : ",self.Rmax_in_AU
    print "kappa_ext              : ",self.kappa_ext
    print "RBIN                   : ",self.RBIN

  #----------------------------------------------------------------
  def get_tau_for_theta_equal_to(self,theta,tau_max = 1e10):
    """
      As described. I use a new algorism to derive Theta_max without an array, using
      the same radiative transfere algorism as the main code.

      Added on 4/13/2012, Hiro Takami

      ### Input paramters ###
        theta   ... as desribe (in deg)
        tau_max ... maximum value for tau. The radiative transfer calculations stops
                    at this value. (Default: 1e10)

      ### Output parameters ###
        tau (as a returned parameter)

    """

    XYZ=D.new_darr(3) ; UVW=D.new_darr(3) ; newXYZ=D.new_darr(3);

    D.darr_setitem(XYZ,0,0.)
    D.darr_setitem(XYZ,1,0.)
    D.darr_setitem(XYZ,2,0.)
    D.darr_setitem(UVW,0,cos(radians(theta)))
    D.darr_setitem(UVW,1,0.)
    D.darr_setitem(UVW,2,sin(radians(theta)))

    photon_is_out_of_nebula = D.NewPosition4givenTau(tau_max,newXYZ, XYZ, UVW)

    D.delete_darr(XYZ) ; D.delete_darr(UVW) ; D.delete_darr(newXYZ)

    return D.cvar.tau

  #----------------------------------------------------------------
  def get_Theta_max(self,tau_crit,d_theta=1):
    """
      As described. I use a new algorism to derive Theta_max without an array, using
      the same radiative transfere algorism as the main code.

      Added on 3/13/2012, Hiro Takami
      Moved to Shakura_Sunyaev on 4/27/2012

      ### Input paramters ###
        tau_crit ... critical value for tau
        d_theta  ... step for determining Theta_max (in deg.)

      ### Output parameters ###
        Theta_max (as a returned parameter)

    """

    theta_deg               = 0
    photon_is_out_of_nebula = 0

    XYZ=D.new_darr(3) ; UVW=D.new_darr(3) ; newXYZ=D.new_darr(3);

    while photon_is_out_of_nebula == 0:
      cost=cos(radians(theta_deg))
      sint=sin(radians(theta_deg))

      D.darr_setitem(XYZ,0,0.)
      D.darr_setitem(XYZ,1,0.)
      D.darr_setitem(XYZ,2,0.)
      D.darr_setitem(UVW,0,cost)
      D.darr_setitem(UVW,1,0.)
      D.darr_setitem(UVW,2,sint)

      photon_is_out_of_nebula = D.NewPosition4givenTau(tau_crit,newXYZ, XYZ, UVW)

      theta_deg += d_theta

    D.delete_darr(XYZ) ; D.delete_darr(UVW) ; D.delete_darr(newXYZ)

    return theta_deg

  #----------------------------------------------------------------
  def get_lines_for_tau_from_star(self,tau,d_theta=0.1):
    """
      As described. 

      Comments added on 4/13/2012, Hiro Takami
      Moved to Shakura_Sunyaev on 4/27/2012

      ### Input paramters ###
        tau     ... as described
        d_theta ... step for angle (in deg.)

      ### Output parameters ###
        arrays for theta, r and z to be able to draw lines

    """

    XYZ=D.new_darr(3) ; UVW=D.new_darr(3) ; newXYZ=D.new_darr(3);

    theta_deg               = 0
    theta_deg_array         = []
    r_array                 = []
    z_array                 = []
    photon_is_out_of_nebula = 0

    while photon_is_out_of_nebula == 0:
      cost=cos(radians(theta_deg))
      sint=sin(radians(theta_deg))

      D.darr_setitem(XYZ,0,0.)
      D.darr_setitem(XYZ,1,0.)
      D.darr_setitem(XYZ,2,0.)
      D.darr_setitem(UVW,0,cost)
      D.darr_setitem(UVW,1,0.)
      D.darr_setitem(UVW,2,sint)

      photon_is_out_of_nebula = D.NewPosition4givenTau(tau, newXYZ, XYZ, UVW)
      r = D.darr_getitem(newXYZ, 0)
      z = D.darr_getitem(newXYZ, 2)

      theta_deg_array.append(theta_deg)
      r_array.append(r)
      z_array.append(z)

      theta_deg += d_theta

    D.delete_darr(XYZ) ; D.delete_darr(UVW) ; D.delete_darr(newXYZ)

    return array(theta_deg_array),array(r_array),array(z_array)

  #----------------------------------------------------------------
  def has_param(self,param):
    return dir(self).count(param) > 0


#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------

if __name__ == "__main__":

  HH30_auto=Shakura_Sunyaev(parameters='Wood+02, 1.65 um',Rmin_in_AU=1)

