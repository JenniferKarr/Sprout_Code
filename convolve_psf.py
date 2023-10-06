#!/usr/bin/env python

from pylab import *
from astropy.io import fits as pyfits
import congrid
from  scipy.signal import *


### Module ###----------------------------------------------------------------------------

def convolve(data,psf,conflag,data_pixel_binning=1):

  """

  This file is part of Sprout.

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


   -----

   This is a convolution routine for the data. There are two options for the 
   PSF, currently only the first is used. 

   Version 1.0 (release) Last Updated May 25 2013


   """

  ##DO CONVOLUTION USING SCIPY 
  if(conflag==1):

    data_convolve=convolve2d(data,psf,mode='same',boundary='symm')

  ##DO CONVOLUTION USING HIRO'S METHOD 
  if (conflag==0):

    nx_data,ny_data     = data.shape
    nx_data_scaled      = nx_data*data_pixel_binning
    ny_data_scaled      = ny_data*data_pixel_binning
    nx_data_half        = int(nx_data*0.5)
    ny_data_half        = int(ny_data*0.5)
    nx_data_scaled_half = int(nx_data_scaled*0.5)
    ny_data_scaled_half = int(ny_data_scaled*0.5)
    
    nx_psf,ny_psf       = psf.shape
    nx_psf_half         = int(nx_psf*0.5)
    ny_psf_half         = int(ny_psf*0.5)
    
    im_psf=zeros((nx_data_scaled,ny_data_scaled),dtype=float)
    im_psf[:nx_psf_half+1,:ny_psf_half+1]=copy(psf[-(nx_psf_half+1):,-(ny_psf_half+1):])
    im_psf[-nx_psf_half: ,-ny_psf_half:] =copy(psf[ :nx_psf_half    ,:ny_psf_half])
    im_psf[:nx_psf_half+1,-ny_psf_half:] =copy(psf[-(nx_psf_half+1):,:ny_psf_half])
    im_psf[-nx_psf_half: ,:ny_psf_half+1]=copy(psf[ :nx_psf_half    ,-(ny_psf_half+1):])
    
    im_obj_fft=fft2(data)
    im_psf_fft=fft2(im_psf)
    
    im_psf_fft2=zeros((nx_data,ny_data),complex)
    im_psf_fft2[:nx_data_half,:ny_data_half]  =copy(im_psf_fft[:nx_data_half,:ny_data_half])
    im_psf_fft2[-nx_data_half:,-ny_data_half:]=copy(im_psf_fft[-nx_data_half:,-ny_data_half:])
    im_psf_fft2[:nx_data_half,-ny_data_half:] =copy(im_psf_fft[:nx_data_half,-ny_data_half:])
    im_psf_fft2[-nx_data_half:,:ny_data_half] =copy(im_psf_fft[-nx_data_half:,:ny_data_half])
    
    data_convolve=real(ifft2(im_obj_fft*im_psf_fft2))



  return data_convolve
