#!/usr/bin/env python

from Tkinter import * 
import tkMessageBox 
import tkFileDialog 
from tkFont import Font 
import pickle
import numpy
#from post_program import convert_fits
import Pmw
import pylab
import os

class view_directories(Frame): 
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

   This is the routine for a quick look at the parametesr for sets of models
   created by the program run_models.py 

   Version 1.0 (release) Last Updated May 25 2013


"""

  def __init__(self,master=None): 
    Frame.__init__(self,master) 

    self.loaded=0
    fm=Frame(self, relief=RAISED, bd=1) 
    mb2=Button(fm,text='Quit',command=self.quit).pack(side=LEFT,anchor=W)
    mb2=Button(fm,text='Choose Directory',command=self.load_dir).pack(side=LEFT,anchor=W)
    mb2=Button(fm,text='Save to File',command=self.save_text).pack(side=LEFT,anchor=W)
    fm.pack(side=TOP,anchor=W,fill=X) 

    fm=Frame(self, relief=RAISED, bd=1)
    self.text=Text(fm,height=30,width=100)
    self.text.pack()
    fm.pack(side=TOP,anchor=W,fill=X)

  def load_dir(self):

      if(self.loaded == 1):
          self.text.delete(1.0, END)
          #self.text.insert(END, text)
          #self.text.mark_set(INSERT, 1.0)

      ##GET DIRECTORY NAME
      dirname=tkFileDialog.askdirectory()


      ##IF DIRECTORY EXISTS
      if os.path.exists(dirname): 
          self.loaded=1
          f = open(dirname+'/run_info.pkl','r')
          myObject=pickle.load(f)
          f.close()

          ##SET VALUES
          h_array=myObject[0]
          rho_array=myObject[1]
          beta_array=myObject[2]
          mod=myObject[3]
          npix=int(mod['igrid'][2])
          dist=float(mod['dist'][2])
          resamp=float(mod['resamp'][2])
          angle=float(mod['incl'][2])
          dust=myObject[4]
          ss=myObject[5]
          angles=myObject[6]
      
          options=myObject[7]
          vars=myObject[8]


          self.text.insert(INSERT,"Directory = "+dirname+"\n\n")

          self.text.insert(INSERT,str(npix)+" by "+str(npix)+" pix model at "+str(dist)+" pc, with a pixel resample of "+str(resamp)+"\n")
          self.text.insert(INSERT,"Data Pixel Scale= "+str(mod['pixscale'][2])+" arcsec\n")

          self.text.insert(INSERT,"Running "+str(len(h_array)*len(beta_array)*len(rho_array))+" models total.\n\n")

          self.text.insert(INSERT,"Model Parameters\n")
          self.text.insert(INSERT,"  "+str(mod['ntot'][2])+"e6 photons\n")
          self.text.insert(INSERT,"  Specific Viewing Angle = "+str(mod['incl'][2])+" degrees \n\n")

          self.text.insert(INSERT,"Dust Parameters\n")
          self.text.insert(INSERT,"  Reference Wavelength = "+str(dust['wave'][2])+" microns \n")
          self.text.insert(INSERT,"  "+dust['dust_model'][2]+" dust model\n")

          self.text.insert(INSERT,"\n")
          self.text.insert(INSERT,"Disk Model = "+str(ss['disk_model'][2])+"\n")
          self.text.insert(INSERT,"  Minimum Radius = "+str(ss['minr'][2])+" AU\n")
          self.text.insert(INSERT,"  Maximum Radius = "+str(ss['maxr'][2])+" AU\n")
          self.text.insert(INSERT,"  Stellar Radius = "+str(ss['r'][2])+" Rsol\n\n")


          
          ##BETA INFO
          self.text.insert(INSERT,"Beta values \n")

          if(len(beta_array) <= 2):
            self.text.insert(INSERT,"   Beta varies from "+str(numpy.min(beta_array))+" to "+str(numpy.max(beta_array))+" in "+str(len(beta_array))+" steps.\n")
          elif(len(beta_array) > 2 and vars['beta'][4]==0):
            self.text.insert(INSERT,"   Beta varies from "+str(numpy.min(beta_array))+" to "+str(numpy.max(beta_array))+" in "+str(len(beta_array))+" linear steps.\n")
          else:
            self.text.insert(INSERT,"   Beta varies from "+str(numpy.min(beta_array))+" to "+str(numpy.max(beta_array))+" in "+str(len(beta_array))+" logarithmic steps.\n")
          self.text.insert(INSERT,"\n")

          self.text.insert(INSERT,"h values \n")

          if(options['fixh']==1):
            self.text.insert(INSERT,"   h is fixed at R*. Units of h are R*\n")
          else:
            self.text.insert(INSERT,"   h is fixed at R = "+str(options['h0_val'])+" AU. Units of h are AU.\n")


          if(len(h_array) <= 2):
            self.text.insert(INSERT,"   h0 varies from "+str(numpy.min(h_array))+" to "+str(numpy.max(h_array))+" in "+str(len(h_array))+" steps.\n")
          elif(len(h_array) > 2 and vars['h0'][4]==0):
            self.text.insert(INSERT,"   h0 varies from "+str(numpy.min(h_array))+" to "+str(numpy.max(h_array))+" in "+str(len(h_array))+" linear steps.\n")
          else:
            self.text.insert(INSERT,"   h0 varies from "+str(numpy.min(h_array))+" to "+str(numpy.max(h_array))+" in "+str(len(h_array))+" logarithmic steps.\n")

          self.text.insert(INSERT,"\n")


          self.text.insert(INSERT,"Rho values \n")
          if(options['fixrho']==1):
            self.text.insert(INSERT,"   rho is fixed at R*."+"\n")
          elif(options['fixrho']==2):
            self.text.insert(INSERT,"   rho is fixed at R = "+str(options['rho_val'])+" AU."+"\n")
          elif(options['fixrho']==3):
            self.text.insert(INSERT,"   rho is scaled to satisfy tau = "+str(options['tau_val'])+" at theta = "+str(options['ang_val'])+"\n")
          elif(options['fixrho']==4):
            self.text.insert(INSERT,"   rho is scaled to satisfy mdisk = "+str(options['mass_val'])+"\n")

          if(options['fixrho']==1 or options['fixrho']==2):
            if(len(rho_array) <= 2):
              self.text.insert(INSERT,"   rho0 varies from "+str(numpy.min(rho_array))+" to "+str(numpy.max(rho_array))+" in "+str(len(rho_array))+" steps.\n")
            elif(len(rho_array) > 2 and vars['rho'][4]==0):
              self.text.insert(INSERT,"   rho0 varies from "+str(numpy.min(rho_array))+" to "+str(numpy.max(rho_array))+" in "+str(len(rho_array))+" linear steps.\n")
            else:
              self.text.insert(INSERT,"   rho0 varies from "+str(numpy.min(rho_array))+" to "+str(numpy.max(rho_array))+" in "+str(len(rho_array))+" logarithmic steps.\n")

          self.text.insert(INSERT,"\n")


  def save_text(self):
      save_file=tkFileDialog.asksaveasfilename()
      #if os.file.exists(save_file): 

      contents=self.text.get(1.0,END)
      f = open(save_file, "w")
      f.write(contents.rstrip())
      f.close
      


a=view_directories() 
a.pack() 
a.mainloop() 
