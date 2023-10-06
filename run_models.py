#!/usr/bin/env python

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

   This is a GUI to generate a set of scattered light images for circumstellar 
   disks. The results can be viewed using the program view_models.py

   Version 1.0 (release) Last Updated May 25 2013

"""

import matplotlib 
matplotlib.use("TkAgg")

from pylab import * 
from Tkinter import * 
import tkMessageBox 
import tkFileDialog 
from tkFont import Font 
import pickle
import numpy
import Pmw
import os
import time

from MC_SS_func import MC_simulation
import codecs
from textfile_window import *

##-----------------------------------------------------------------------

class counter_window(Frame): 
  """

  A Counter frame.  Shows the number of the model being run, the total number of models, 
   and the estimated time to completion (after the first model has run).
  
  To call:  counterwin=counter_window(n_models)  #where N_models is an integer.
            counterwin.pack(side=TOP,anchor=W,fill=X)
            counterwin.update_idletasks()
            counterwin.update()

  To update: counterwin.numdone['text']=str(num)  #where num is an integer and the current model number
             counterwin.timeleft['text']=str(time) #where time is the estimated completion time in minutes. 
             counterwin.update_idletasks()
             counterwin.update()

  """

  def __init__(self,tot,master=None): 
    Frame.__init__(self,master) 

    ##CREATE WINDOW AND INITIALIZE VARIABLES

    top=Toplevel(self) 
    top.title('Progress') 

    ##CREATE THE FRAME

    fm=Frame(top) 
    lab=Label(fm,text="Running Program").pack(side=LEFT,anchor=W) 
    fm.pack()
    fm=Frame(top) 
    lab=Label(fm,text="Model ").pack(side=LEFT,anchor=W) 

    self.numdone=Label(fm,text="")
    self.numdone.pack(side=LEFT,anchor=W) 
    lab=Label(fm,text=" of ").pack(side=LEFT,anchor=W) 
    lab=Label(fm,text=str(tot),width=4).pack(side=LEFT,anchor=W) 
    fm.pack()
    fm=Frame(top) 
    lab=Label(fm,text="Time to Completion ~").pack(side=LEFT,anchor=W) 
    self.timeleft=Label(fm,text="")
    self.timeleft.pack(side=LEFT,anchor=W) 
    lab=Label(fm,text=" minutes.").pack(side=LEFT,anchor=W) 
    fm.pack()

    ##SET TOTAL VALUES AND INITIALIZE TIME

    self.numdone['text']="0"
    self.timeleft['text']="--"

##---------------------------------------------------------------------

class MC_parameter_window(Frame): 
  """
  
  The main GUI window class. 

"""


  def __init__(self,master=None): 
    Frame.__init__(self,master) 

     ##INITTIALIZE VARIABLES. 

     ##GENERAL VARIABLES

    self.master.title('Parameters for Scattering Code') 
    
    self.savedat=IntVar()  ##FLAG FOR SAVING IN TEXT FORMAT 

     ##VARIABLES FOR SIMULATION

    self.mod={"ntot":["Number of Photons","/1e6",1,StringVar(),0],  #TOTAL NUMBER OF PHOTONS
              "incl":["User Defined Viewing Angle","degrees",35,StringVar(),0],  #USER SPECIFIED VIEWING ANGLE
              "irange":["Range for User Defined Viewing Angle","degrees",10,StringVar(),0], #RANGE FOR VIEWING ANGLE
              "igrid":["Grid","",99,StringVar(),0], #NUMBER OF PIXELS IN EACH AXIS (CALCLULATED)
              "maxsca":["Maximum Number of Scatterings","",3,StringVar(),0], #MAXIMUM NUMBER OF SCATTERS
              "dist":["Distance","pc",144,StringVar(),0], #DISTANCE OF SOURCE
              "resamp":["Pixel Sampling","",2,StringVar(),0], #RESAMPLING FACTOR. 2=PIXELS TWICE THE SIZE
              "dpath":["Step Size for Photon Path","AU",0.5,StringVar(),0], #STEP SIZE FOR INTEGRATION
              "pixscale":["Pixel Scale","arcsec",0.00948,StringVar(),0], #PIXEL SCALE OF OBSERVATIONS BEFORE RESAMPLING
              "pixpad":["Image Padding","percent",10.,StringVar(),0], #SIZE OF BORDER TO ADD TO THE FINAL IMAGE
              "logfile":["Name of Log File","","log.txt",StringVar(),0]} ##LOGFILE
              

    ##MODEL VARIABLES TO SAVE TO FILE. THIS IS NECESARY BECAUSE THE STRINGVAR() ETC 
    ##OF SELF.MOD CANNOT BE STORED IN A FILE. OTHERWISE, THIS IS THE SAME AS SELF.MOD 

    self.mod_store={"ntot":["Number of Photons","/1e6",0],
                    "incl":["User Defined Viewing Angle","degrees",0],
                    "irange":["Range for User Defined Viewing Angle","degrees",0],
                    "igrid":["Grid","",0],
                    "maxsca":["Maximum Number of Scatterings","",0],
                    "dist":["Distance","pc",0],
                    "resamp":["Pixel Sampling","",0],
                    "dpath":["Step Size for Photon Path","AU",0],
                    "pixscale":["Pixel Scale","arcsec",0],
                    "pixpad":["Image Padding","percent",0],
                    "logfile":["Log file","",0]}


    ##DEFINE WHICH VARIABLES ARE IN THE ADVANCED (CONST) CATEGORY OR ARE LIKELY TO BE CHANGED BY 
    ##AN AVERAGE USER (CHANGE)

    self.mod_const=["maxsca","dpath","irange","pixpad"]     ## LIST OF VARIABLES THAT USER PROBABLY WON'T CHANGE (ADVANCED)
    self.mod_change=["ntot","incl","dist","resamp","pixscale","logfile"]  ##LIST OF VARIABLES THAT USER MAY CHANGE (REGULAR WINDOW)
    self.mod_int=["maxsca"]
    self.mod_float=["dpath","irange","pixpad","ntot","incl","dist","resamp","pixscale"]

     ##VARIABLES FOR DISK SS MODEL
     
    self.ss={"disk_model":["Disk Model","","User Specified",StringVar(),0],  #WHICH DISK MODEL (ARBITRARY OR PRESET)
             "tau_critical":[u"\u03c4_crit","",0.001,StringVar(),0], #TAU CRIT VALUE ??
             "tmax":[u"\u03d1_max","",-1,StringVar(),0], #T MAX VALUE ??
             "minr":["Disk Minimum Radius","AU",5,StringVar(),0], #MINIMUM RADIUS OF DISK
             "maxr":["Disk Maximum Radius","AU",100,StringVar(),0], #MAXIMUM RADIUS OF DISK
             "r":["Stellar radius","R_sol",1.2,StringVar(),0]} #STELLAR RADIUS
        
    ##DISK VARIABLES TO SAVE TO FILE, AS ABOVE

    self.ss_store={"disk_model":["Disk Model","",0],
                   "tau_critical":["Tau Critical","",0],
                   "tmax":["T Max","",0],
                   "minr":["Minimum Radius","AU",0],
                   "maxr":["Maximum Radius","AU",0],
                   "r":["Stellar radius","Rsol",0]}
        
    ##VARIABELS FOR THE MODEL GRID. SAME FORMAT AS ABOVE VARIABLES, BUT WITH FOUR VARIABLES; MIN VALUE, MAX VALUE, NUMBER
    ##OF STEPS, AND LOG/LINEAR FLAG.

    self.ss_var={"h0":["Disk Scale Height h_o",0.01,0.02,3,0,
                       StringVar(""),StringVar(""),StringVar(""),IntVar(),
                       [],[],[]], #SCALE HEIGHT VARIABLES
                 "beta":["Disk Flaring Power "+u"\u03b2",1.4,1.6,3,0,StringVar(""),StringVar(""),StringVar(""),IntVar(),
                         [],[],[]], #FLARING VARIABLES
                 "rho":["Initial Density "+u"\u03c1 (g cm-3)",10**(-7),10**(-5),3,0,StringVar(""),StringVar(""),StringVar(""),IntVar(),
                        [],[],[]]} #DENSITY VARIABLES

    ##GRID VARIABLES TO SAVE TO FILE, AS ABOVE

    self.ss_var_store={"h0":["Disk Scale Height H",0,0,0,0],
                       "beta":["Disk Flaring Power",0,0,0,0],
                       "rho":["Initial Density +"u"u03c1(g cm-3)",0,0,0,0]}

    ##DEFINE WHICH VARIABLES ARE IN THE ADVANCED (CONST) CATEGORY OR ARE LIKELY TO BE CHANGED BY 
    ##AN AVERAGE USER (CHANGE)

    self.ss_const=["tmax","tau_critical"] ## ADVANCED PARAMTETERS
    self.ss_change=['minr','maxr','r']  ## MAIN WINDOW PARAMETERS
    self.ss_float=["tau_critical","tmax",'minr','maxr','r']
    #self.ss_float=['minr','maxr','r']

    ##DICTIONARY OF VARIABLES FOR THE OPTIONS FOR RHO AND H

    self.ss_options={'fixh':[1,IntVar()], #FLAG FOR FORMAT OF HO (AT R* OR X AU)
                     'fixrho':[1,IntVar()], #FLAG FOR FORMAT OF HO (AT R* OR X AU, OR TAU, OR MASS)
                     'h0_val':["",StringVar()], #RADIUS TO FIX H0 IF FIXH=2
                     'rho_val':["",StringVar()], #RADIUS TO FIX RHO IF FIXRHO=2
                     'tau_val':["",StringVar()], #TAU VAULUE OF FIXRHO=3
                     'ang_val':["",StringVar()], #ANGLE VALUE IF FIXRHO=3
                     'mass_val':["",StringVar()]} #DISK MASS IF RIXRHO=4

    ##DICTIONARY TO SAVE TO FILE, AS ABOVE

    self.ss_options_store={'fixh':0,
                     'fixrho':0,
                     'h0_val':0,
                     'rho_val':0,
                     'tau_val':0,
                     'ang_val':0,
                     'mass_val':0}

    #self.which_disk=StringVar()  ## WHICH DISK MODEL TO USE

    self.ss_entry=[]   ## USED FOR MASKING OUT NON RELEVENT ENTRY FIELDS
    self.ss_entry1=[] 

     ##VARIABLES FOR DUST

    ##VARIABLES FOR DUST GRAIN MODEL

    self.dust={"kappa_ext":[u"\u03ba_ext","/g",0,StringVar(),0],  ##OPACITY TO OVERRIDE MODEL
               "dust_model":["Dust Model","","KMH94, R_V=3.1",StringVar(),0], ##DUST MODEL
               "wave":["Wavelength","microns",1.65,StringVar(),0]} ##WAVELENGTH OF MODEL
 
    ##DUST VARIABLES FOR SAVING. 

    self.dust_store={"wave":["Wavelength","microns",0],
                     "dust_model":["Dust Model","",0],
                     "kappa_ext":[u"\u03ba_ext","/g",0]}

    self.dust_const=['kappa_ext'] ##ADVANCED PARAMTERS
    self.dust_change=[] ##REGULAR PARAMETERS (WAVELENGTH ADN DUST_MODEL ARE DROP-DOWN MODELS)

   ##---------------------------------------------------------------------
   ##SET UP THE GUI

    ##TOP BAR OF MENU

    fm=Frame(self, relief=RAISED, bd=1) 
    mb2=Button(fm,text='Quit',command=self.quit).pack(side=LEFT,anchor=W) 
    mb2=Button(fm,text='Help',command=self.help_all).pack(side=LEFT,anchor=W)
    fm.pack(side=TOP,anchor=W,fill=X) 

    ##SAVE/LOAD/RUN OPTIONS ON TOP MENU

    fm=Frame(self, relief=RAISED, bd=1) 
    mb2=Button(fm,text='Run Simulations',command=self.start_loop).pack(side=LEFT,anchor=W)
    mb2=Button(fm,text='Save Parameters',command=self.save_params).pack(side=LEFT,anchor=W)
    mb2=Button(fm,text='Load Parameters',command=self.load_params).pack(side=LEFT,anchor=W)
    mb2=Button(fm,text='Check Input',command=self.check_params).pack(side=LEFT,anchor=W)
    mb2=Checkbutton(fm,text='Save Results as Text',var=self.savedat).pack(side=LEFT,anchor=W)
    fm.pack(side=TOP,anchor=W,fill=X) 


    ##USES PMW NOTEBOOK TO SPLIT THE PARAMETERS INTO MULTIPLE TABS

    notebook=Pmw.NoteBook(self)  ##INITIALIZE THE NOTEBOOK

     ##---------------------------------------------------------------------

    ## FIRST PAGE; GENERAL PARAMETERS FOR THE MODEL

    page1 = notebook.add('Model Parameters')
    notebook.tab('Model Parameters').focus_set()

    ##TOP BAR

    ##ADD A BIT OF SPACE FOR AESTHETIC REASONS

    model=Frame(page1) 
    title=Label(model,text=" ").pack(side=LEFT,anchor=W) 
    title=Label(model,text=" ").pack(side=LEFT,anchor=W) 
    model.pack(side=TOP,anchor=W,fill=X) 

    ##CLEAR/RESET OPTIONS
    model=Frame(page1) 
    title1=Label(model,text="   Model Parameters   ",background="grey").pack(side=LEFT,anchor=W) 
    mb2=Button(model,text='Clear Values',command=self.clear_model).pack(side=LEFT,anchor=W)
    mb2=Button(model,text='Reset Values to Defaults',command=self.reset_model).pack(side=LEFT,anchor=W)
    model.pack(side=TOP,anchor=W,fill=X) 
   

     ##CYCLE THROUGH THE COMMON PARAMETERS AND SET UP ENTRIES

    for k in self.mod_change:
      l1=Frame(page1)
      self.label=Label(l1,text=self.mod[k][0],width=30).pack(side=LEFT,anchor=W)
      self.entry=Entry(l1,textvariable=self.mod[k][3],width=20).pack(side=LEFT,anchor=W)
      self.label=Label(l1,text=self.mod[k][1],width=20).pack(side=LEFT,anchor=W)
      l1.pack(side=TOP,anchor=W,fill=X) 
      self.mod[k][3].set(self.mod[k][2])

    ##SET TO THE DEFAULT VALUES
    self.reset_model()

     ##---------------------------------------------------------------------

    ## SECOND PAGE;  PARAMETERS FOR THE DUST


    page2 = notebook.add('Dust Parameters')

    ##TOP BAR

    ##A BIT OF PADDING FOR AESTHETIC REASONS

    dust=Frame(page2) 
    title=Label(dust,text=" ").pack(side=LEFT,anchor=W) 
    title=Label(dust,text=" ").pack(side=LEFT,anchor=W) 
    dust.pack(side=TOP,anchor=W,fill=X) 

    ##RESET PARAMETERS
    dust=Frame(page2) 
    title=Label(dust,text="   Physical Properites of Dust   ",background="grey").pack(side=LEFT,anchor=W) 
    mb2=Button(dust,text='Reset Values to Defaults',command=self.reset_dust).pack(side=LEFT,anchor=W)
    dust.pack(side=TOP,anchor=W,fill=X) 

    ##DROP-DOWN MENU FOR PRE-SET MODELS
    l1=Frame(page2)
    label=Label(l1,text="Dust Model  ").pack(side=LEFT,anchor=W)
    opt=OptionMenu(l1,self.dust['dust_model'][3],
                   'KMH94, R_V=3.1',
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
                   'Cotera+01 x 15, amorphous (J98,1000C)')




    opt.pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 


    ##DROP DOWN MENU FOR WAVELENGTH

    l1=Frame(page2)
    label=Label(l1,text="Wavelength  ").pack(side=LEFT,anchor=W)
    opt=OptionMenu(l1,self.dust['wave'][3],
                   '0.55',
                   '1.25',
                   '1.65',
                   '2.2')
    opt.pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 
      
    ##FILL IN THE DEFAULT VALUES

    self.reset_dust()

     ##---------------------------------------------------------------------

    ## THIRD PAGE;  PARAMETERS FOR THE DENSITY DISTRIBUTION

    page3 = notebook.add('Disk Parameters')

    ##A BIT OF PADDING FOR AESTHETIC REASONS

    geom_window=Frame(page3)
    title=Label(geom_window,text=" ").pack(side=LEFT,anchor=W) 
    title=Label(geom_window,text=" ").pack(side=LEFT,anchor=W) 
    geom_window.pack(side=TOP,anchor=W,fill=X) 
    geom_window=Frame(page3)

    ##RESET/CLEAR OPTIONS

    title=Label(geom_window,text="   Disk Model   ",background="grey").pack(side=LEFT,anchor=W) 
    mb2=Button(geom_window,text='Clear Values',command=self.clear_geom).pack(side=LEFT,anchor=W)
    mb2=Button(geom_window,text='Reset Values to Defaults',command=self.reset_geom).pack(side=LEFT,anchor=W)
    geom_window.pack(side=TOP,anchor=W,fill=X) 

    ## SET UP THE ENTRIES FOR THE COMMON VARIABLES FOR THE DISK

    for k in self.ss_change:
      l1=Frame(page3)
      label=Label(l1,text=self.ss[k][0],width=30).pack(side=LEFT,anchor=W)
      entry=Entry(l1,textvariable=self.ss[k][3],width=20)
      entry.pack(side=LEFT,anchor=W)
      self.ss_entry.append(entry)
      label=Label(l1,text=self.ss[k][1],width=20).pack(side=LEFT,anchor=W)
      l1.pack(side=TOP,anchor=W,fill=X) 
      self.ss[k][3].set(self.ss[k][2])


    ## ENTRIES FOR THE VALUES FOR THE GRID OF MODELS

    l1=Frame(page3)
    label=Label(l1,text="Grid Values",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 
    
    ##SET UP THE WINDOWS FOR MIN/MAX/STEPS/LOG

    for k in self.ss_var:
        l1=Frame(page3)
        self.label=Label(l1,text=self.ss_var[k][0],width=30).pack(side=LEFT,anchor=W)
        for i in range(3):
          entry=Entry(l1,textvariable=self.ss_var[k][i+5],width=10)
          entry.pack(side=LEFT,anchor=W)
          self.ss_var[k][9+i].append(entry)
        entry=Checkbutton(l1,text="Log Steps",var=self.ss_var[k][8])
        entry.pack(side=LEFT,anchor=W)          
        l1.pack(side=TOP,anchor=W,fill=X) 

    ##OPTIONS FOR H_0

    l1=Frame(page3)
    label=Label(l1,text="Options for h_o",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 

    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Fix h_o at R* (h in R*)",variable=self.ss_options['fixh'][1],value=1,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 
    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Fix h_o at Arbitrary Radius (h in AU)",variable=self.ss_options['fixh'][1],value=2,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    self.harb=Entry(l1,textvariable=self.ss_options["h0_val"][1],width=10)
    self.harb.pack(side=LEFT,anchor=W)  
    self.label=Label(l1,text="AU").pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 

    ##OPTIONS FOR RHO

    l1=Frame(page3)
    label=Label(l1,text="Options for "+u"\u03c1 ",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 

    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Fix "+u"\u03c1"+"_o at R*",variable=self.ss_options['fixrho'][1],value=1,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 

    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Fix "+u"\u03c1"+"_o at Arbitrary Radius",variable=self.ss_options['fixrho'][1],value=2,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    self.rhoarb=Entry(l1,textvariable=self.ss_options["rho_val"][1],width=10)
    self.rhoarb.pack(side=LEFT,anchor=W)  
    self.label=Label(l1,text="AU").pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 

    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Scale "+u"\u03c1"+" for "u"\u03c4=",variable=self.ss_options['fixrho'][1],value=3,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    self.tauarb=Entry(l1,textvariable=self.ss_options["tau_val"][1],width=10)
    self.tauarb.pack(side=LEFT,anchor=W)  
    self.label=Label(l1,text="(optical depth) at "+u"\u03d1=").pack(side=LEFT,anchor=W)  
    self.angarb=Entry(l1,textvariable=self.ss_options["ang_val"][1],width=10)
    self.angarb.pack(side=LEFT,anchor=W)  
    self.label=Label(l1,text="degrees (over-rides density grid option)").pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 

    l1=Frame(page3)
    self.check=Radiobutton(l1,text="Fix Disk Mass at ",variable=self.ss_options['fixrho'][1],value=4,command=self.check_disk_state1).pack(side=LEFT,anchor=W)  
    self.massarb=Entry(l1,textvariable=self.ss_options["mass_val"][1],width=10)
    self.massarb.pack(side=LEFT,anchor=W)  
    self.label=Label(l1,text="solar masses (over-rides density grid option)").pack(side=LEFT,anchor=W)  
    l1.pack(side=TOP,anchor=W,fill=X) 


    ##FILL IN THE DEFAULT VALUES
    self.reset_geom()

     ##---------------------------------------------------------------------

    ##ADVANCED PARAMETERS FOR USERS

    page4 = notebook.add('Advanced Parameters')


    ##ADVANCED MODEL PARAMETERS

    l1=Frame(page4)
    title=Label(l1,text="").pack(side=LEFT,anchor=W) 
    title=Label(l1,text="").pack(side=LEFT,anchor=W) 
    label=Label(l1,text="Advanced Model Parameters",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 

    for k in self.ss_const:
      l1=Frame(page4)
      self.label=Label(l1,text=self.ss[k][0],width=30).pack(side=LEFT,anchor=W)
      entry=Entry(l1,textvariable=self.ss[k][3],width=20)
      entry.pack(side=LEFT,anchor=W)
      self.ss_entry.append(entry)
      self.label=Label(l1,text=self.ss[k][1],width=20).pack(side=LEFT,anchor=W)
      l1.pack(side=TOP,anchor=W,fill=X) 
      self.ss[k][3].set(self.ss[k][2])


    ##ADVANCED MODEL PARAMETESR

    l1=Frame(page4)
    title=Label(l1,text="").pack(side=LEFT,anchor=W) 
    title=Label(l1,text="").pack(side=LEFT,anchor=W) 
    label=Label(l1,text="Advanced Model Parameters  ",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 

    for k in self.mod_const:
      l1=Frame(page4)
      self.label=Label(l1,text=self.mod[k][0],width=30).pack(side=LEFT,anchor=W)
      self.entry=Entry(l1,textvariable=self.mod[k][3],width=20).pack(side=LEFT,anchor=W)
      self.label=Label(l1,text=self.mod[k][1],width=20).pack(side=LEFT,anchor=W)
      l1.pack(side=TOP,anchor=W,fill=X) 
      self.mod[k][3].set(self.mod[k][2])

   ##ADVANCED DUST PARAMETESR

    l1=Frame(page4)
    label=Label(l1,text="Advanced Dust Parameters",background="grey").pack(side=LEFT,anchor=W)
    l1.pack(side=TOP,anchor=W,fill=X) 

    for k in self.dust_const:
      l1=Frame(page4)
      self.label=Label(l1,text=self.dust[k][0],width=30).pack(side=LEFT,anchor=W)
      self.entry=Entry(l1,textvariable=self.dust[k][3],width=20).pack(side=LEFT,anchor=W)
      self.label=Label(l1,text=self.dust[k][1],width=20).pack(side=LEFT,anchor=W)
      l1.pack(side=TOP,anchor=W,fill=X) 
      self.dust[k][3].set(self.dust[k][2])


    ##NO NEED TO FILL IN DEFAULTS - THIS WILL BE DONE BY THE PREVIOUS CALLS

     ##---------------------------------------------------------------------

    ##FINALIZE AND PACK THE NOTEBOOK
    notebook.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
    notebook.setnaturalsize()

    ##CHECKE THE PARAMETERS CHOSEN, AND GREY OUT UN-AVAILABLE OPTIONS

    self.check_disk_state(1)

     ##---------------------------------------------------------------------

  def check_disk_state1(self):
    
    """

    A simple wrapper for check_disk_state, so you can call it from the menu
    or internally. 

    """

    
    ##SIMPLE WRAPPER FOR CHECK_DISK_STATE, SO YOU CAN CALL IT FROM A FUNCTION
    ##CALL OR NOT. 

    z=self.check_disk_state(1)

     ##---------------------------------------------------------------------

  def check_disk_state(self,value):

    """
    
    Check out the paramaters chosen and adjust the interface - greying out options,
    changing units, etc. 


    """

    ##FIRST, RESET ALL ENTRIES TO ENABLED

    for k in self.ss_var:
      for i in range(3):
        self.ss_var[k][i+9][0].config(state=NORMAL)

    self.harb.config(state=NORMAL)
    self.rhoarb.config(state=NORMAL)
    self.angarb.config(state=NORMAL)
    self.tauarb.config(state=NORMAL)
    self.massarb.config(state=NORMAL)

    ##DISABLE UNUSED WINDOWS FOR HO OPTIONS

    if(self.ss_options['fixrho'][1].get()>=3):
      for i in range(3):
        self.ss_var['rho'][i+9][0].config(state='disabled')

    if(self.ss_options['fixh'][1].get()==1):
      self.harb.config(state='disabled')
    else:
      self.harb.config(state=NORMAL)
    
    ##DISABLE UNUSED WINDOWS FOR RHO OPTIONS

    if(self.ss_options['fixrho'][1].get()==1):
      self.rhoarb.config(state='disabled')
      self.tauarb.config(state='disabled')
      self.angarb.config(state='disabled')
      self.massarb.config(state='disabled')
    if(self.ss_options['fixrho'][1].get()==2):
      self.tauarb.config(state='disabled')
      self.angarb.config(state='disabled')
      self.massarb.config(state='disabled')
    if(self.ss_options['fixrho'][1].get()==3):
      self.rhoarb.config(state='disabled')
      self.massarb.config(state='disabled')
    if(self.ss_options['fixrho'][1].get()==4):
      self.rhoarb.config(state='disabled')
      self.tauarb.config(state='disabled')
      self.angarb.config(state='disabled')

     ##---------------------------------------------------------------------

  def quit(self):

    """
    Quit.

    """

    self.master.destroy()  # quit from the window 

  ##---------------------------------------------------------------------

  def help_all(self): 


    """

    Create a help window, and display the file. 

    """

    ##DISPLAY THE HELP FILE WINDOW

    self.helpwin=textfile_window("Help File","help_run.txt")
    self.helpwin.pack(side=TOP,anchor=W,fill=X)
    self.helpwin.update()

  ##--------------------------------------------------------------------

  def clear_model(self): 

    """

    Clear the Model window entries. 
    
    """

     ##CLEAR THE MODEL PARAMETER WINDOWS

    for k in self.mod_change:
      self.mod[k][3].set("")
    for k in self.mod_const:
      self.mod[k][3].set("")

  ##---------------------------------------------------------------------

  def clear_geom(self): 
 
    """

    Clear the Disk window entries. 
    
    """

    ##CLEAR THE GEOMETRY VALUES UPDATE THE SS_SET FLAG

    for k in self.ss_const:
      self.ss[k][3].set("")
    for k in self.ss_change:
      self.ss[k][3].set("")
    for k in self.ss_var:
      for i in range(4):
        self.ss_var[k][i+5].set("")
    for k in self.ss_options:
      self.ss_options[k][1].set(self.ss_options[k][0])

 ##---------------------------------------------
   
  def reset_dust(self): 

    """

    Reset the Dust parameter windows to default values. 

    """

     ##RESET THE DUST PARAMETER WINDOWS TO DEFAULT VALUES
       
    for k in self.dust:
      self.dust[k][3].set(self.dust[k][2])

  ##---------------------------------------------------------------------

  def reset_model(self): 

    """

    Reset the Model parameter windows to default values. 

    """

     ##RESET THE MODEL PARAMETER WINDOWS TO DEFAULT VALUES

    for k in self.mod_change:
      self.mod[k][3].set(self.mod[k][2])
    for k in self.mod_const:
      self.mod[k][3].set(self.mod[k][2])

  ##---------------------------------------------------------------------
      
  def reset_geom(self): 

    """

    Reset the Disk parameter windows to default values. 

    """

    ##RESET TO DEFAULT VALUES UPDATE THE SS_SET FLAG

    for k in self.ss_const:
      self.ss[k][3].set(self.ss[k][2])
    for k in self.ss_change:
      self.ss[k][3].set(self.ss[k][2])

    for k in self.ss_var:
      for i in range(4):
        self.ss_var[k][i+5].set(self.ss_var[k][i+1])

    for k in self.ss_options:
      self.ss_options[k][1].set(self.ss_options[k][0])

    self.check_disk_state(1)

 ##---------------------------------------------
  
  def check_params(self):

    """
    
    Check the input parameters for validity. Will check for 
    non-numeric input, positive values, etc, and display 
     a warning message or confirmation. If correct, values
     will be stored. 

    """

    ##INITIALIZE VARIABLES
    
    self.goodinput=0

    ##WARNINGS AND FLAGS FOR DUST, GEOMETRY AND MODEL PARAMETERS

    warning=""
 
    modnum=0
    geomnum=0

   ##CHECKS THE DUST PARAMETERS FOR VALIDITY. 

   ##CHECK FOR NUMERIC INPUT

    for k in self.dust_const:
      try:
        x=float(self.dust[k][3].get())
      except:
        warning=warning+"Non Numeric Input in Dust Parameters!\n"

   ##CHECKS THE MODEL PARAMETERS FOR VALIDITY

   ##CHECK FOR NUMERIC INPUT

    for k in self.mod_float:
      try:
        x=float(self.mod[k][3].get())
      except:
        warning=warning+" "+self.mod[k][0]+" "+" must be a floating point number.\n"
        modnum=1
    for k in self.mod_int:
      try:
        x=int(self.mod[k][3].get())
        if((x-float(self.mod[k][3].get()) > 0.00001)):
          warning=warning+" "+self.mod[k][0]+" "+" must be an integer.\n"
      except:
        warning=warning+" "+self.mod[k][0]+" "+" must be an integer.\n"
        modnum=1
        
    if modnum != 1 :
   
   ##IF THE VALUES ARE NUMBERIC, CHECK FOR PHYSICALLY PLAUSIBLE NUMBERS

      if float(self.mod["ntot"][3].get()) <= 0:
        warning=warning+" "+self.mod["ntot"][0]+" must be > 0!\n"
      if float(self.mod["resamp"][3].get()) <= 0:
        warning=warning+" "+self.mod["resamp"][0]+" must be > 0!\n"
      if float(self.mod["dist"][3].get()) <= 0:
        warning=warning+" "+self.mod["dist"][0]+" must be > 0!\n"
      if float(self.mod["incl"][3].get()) < 0 or  float(self.mod["incl"][3].get()) > 90:
        warning=warning+" "+self.mod["incl"][0]+" must be >=0 and <=90 degrees!\n"
      if float(self.mod["irange"][3].get()) < 0 or  float(self.mod["irange"][3].get()) > 90:
        warning=warning+" "+self.mod["irange"][0]+" must be  >=0 and <=90 degrees!\n"
      if float(self.mod["maxsca"][3].get()) <= 0:
        warning=warning+" "+self.mod["maxsca"][0]+" must be > 0!\n"
      if float(self.mod["dpath"][3].get()) <= 0:
        warning=warning+" "+self.mod["dpath"][0]+" must be >0 \n"


   ##CHECKS THE DISK PARAMETERS FOR VALIDITY

   ##CHECK MAIN DISK PARAMETESR

    for k in self.ss_float:
      try:
        x=float(self.ss[k][3].get())
      except:
        warning=warning+" "+self.ss[k][0]+" must be a floating point number.\n"
        geomnum=1

    try:
      if(float(self.ss['tmax'][3].get()) < 0 and float(self.ss['tau_critical'][3].get()) < 0 ):
        warning=warning+"Either "+self.ss['tmax'][0]+" or "+self.ss['taucrit'][0]+" must be set.\n"
    except:
      print 
      

    ##NOW THE GRID VARIABLES - CHECK FOR INTEGERS AND MIN<MAX

    for k in self.ss_var:
      try:
        x=float(self.ss_var[k][5].get())
        if(x < 0):
          warning=warning+"Minimum value for "+self.ss_var[k][0]+" must be > 0\n"
          geomnum=1
      except:
        warning=warning+"Minimum value for "+self.ss_var[k][0]+" must be a floating point number.\n"
        geommnum=1
      try:
        x=float(self.ss_var[k][6].get())
        if(x < 0):
          warning=warning+"Maximum value for "+self.ss_var[k][0]+" must be > 0\n"
          geomnum=1
      except:
        warning=warning+"Maximum value for "+self.ss_var[k][0]+" must be a floating point number.\n"
        geomnum=1
      try:
        x=int(self.ss_var[k][7].get())
        if((x-float(self.ss_var[k][7].get()) > 0.001)):
          warning=warning+"Number of steps for "+self.ss_var[k][0]+" must be an integer.\n"
          geomnum=1
      except:
        warning=warning+"Number of steps for "+self.ss_var[k][0]+" must be an integer.\n"
        geomnum=1

      try:
        if(float(self.ss_var[k][6].get())-float(self.ss_var[k][5].get()) < 0):
          warning=warning+"Maximum value of "+self.ss_var[k][0]+" must be greater than the minimum value.\n"
      except:
        warning=warning

    ##NOW CHECK OPTIONS FOR

    if(self.ss_options['fixh'][1].get()==2):
      try:
        z=float(self.ss_options['h0_val'][1].get())
        if(z < 0):
          warning=warning+"Distance for ho must be a positive floating point number.\n"
      except:
        warning=warning+"Distance for ho must be a positive floating point number.\n"
        
    ##NOW CHECK OPTIONS FOR RHO

    if(self.ss_options['fixrho'][1].get()==2):
      try:
        z=float(self.ss_options['rho_val'][1].get())
        if(z < 0):
          warning=warning+"Distance for "+u"\u03c1 must be a positive floating point number.\n"
      except:
        warning=warning+"Distance for "+u"\u03c1 must be a positive floating point number.\n"

    if(self.ss_options['fixrho'][1].get()==3):
      try:
        z1=float(self.ss_options['tau_val'][1].get())
        if(z1 <=0):
          warning=warning+"Value for "+u"\u03c4 must be a floating point number.\n"
      except:
        warning=warning+"Value for "+u"\u03c4 must be a floating point number.\n"

      try:
        z2=float(self.ss_options['ang_val'][1].get())
        if(z1 <0 or z1 > 90):
          warning=warning+"Value for "+u"\u03d1 must be > 0 and < 90\n"
      except:
        warning=warning+"Value for "+u"\u03d1 must be a floating point number.\n"

    if(self.ss_options['fixrho'][1].get()==4):
      try:
        z1=float(self.ss_options['mass_val'][1].get())
        if(z1 < 0):
          warning=warning+"Value for disk mass must be > 0\n"
      except:
        warning=warning+"Value for disk mass must be a floating point number.\n"

    ##MAKE SURE ONE OF TMAX OR TAUCRIT ARE SET
  


    ##DISPLAY THE RESULTS

    if(warning != "" ):
      tkMessageBox.showwarning("Warning", warning)
    else:
      tkMessageBox.showinfo("","Input Parameters Okay!")
      self.goodinput=1


    ##NOW STORE THE PARAMETERS IN THE SAVE VARAIBLES (_STORE)



    ##FOR THE DUST
    for k in self.dust:
      self.dust_store[k][2]=self.dust[k][3].get()

    ##FOR THE MODEL
    for k in self.mod_const:
      self.mod_store[k][2]=self.mod[k][3].get()
    for k in self.mod_change:
      self.mod_store[k][2]=self.mod[k][3].get()

    ##FOR THE DISK GEOMETRY

    for k in self.ss_const:
      self.ss_store[k][2]=self.ss[k][3].get()
    for k in self.ss_change:
      self.ss_store[k][2]=self.ss[k][3].get()
    #self.ss_store['disk_model'][2]=self.ss['disk_model'][3].get() 
    self.ss_store['disk_model'][2]='User Specified'

    for k in self.ss_var:
      for i in range(4):
        self.ss_var_store[k][i+1]=self.ss_var[k][i+5].get()

    for k in self.ss_options:
      self.ss_options_store[k]=self.ss_options[k][1].get()

    print "zzz",self.ss_var_store['rho'][1]
    print "zzz",self.ss_var_store['rho'][2]
    print "zzz",self.ss_var_store['rho'][3]
    print "zzz",self.ss_var_store['rho'][4]

   ##---------------------------------------------

  def save_params(self):

    """

    Save parameters to an external file. 

    """


    ##CHECK IF THE PARAMETERS ARE VALID

    self.check_params()

    ##IF SO, ASSEMBLE THE DATA STRUCTURE. tHIS CONTAINS THE DICTIONARIES FOR THE MODEL, DUST AND DISK, 
    ##DISK GRID, AND DISK OPTIONS

    if self.goodinput == 1:

      object=[self.mod_store,self.dust_store,self.ss_store,self.ss_var_store,self.ss_options_store]

      ##QUERY FOR THE FILE NAME AND SAVE

      save_window=tkFileDialog.asksaveasfilename()
      open_file=open(save_window,'w')
      pickle.dump(object,open_file)
      open_file.close
        
      tkMessageBox.showinfo("","Successfully Saved!")

  ##---------------------------------------------------------------------

  def load_params(self):

    """

    Load parameters from a save file. 

    """

   ##LOAD MODEL PARAMETERS SAVED IN SAVE_PARAMS

    ##PROMPT USER FOR THE SAVE FILE

    try:
      save_window=tkFileDialog.askopenfilename()
      open_file=open(save_window,'r')
      object=pickle.load(open_file)
      open_file.close
    except:
      tkMessageBox.showwarning("", "This does not appear to be the correct type of file!")
      return

    ##PUT THE SAVE FILE INTO THE RIGHT PARAMETERS

    self.mod_store=object[0]
    self.dust_store=object[1]
    self.ss_store=object[2]
    self.ss_var_store=object[3]
    self.ss_options_store=object[4]

    ##AND FILL IN THE ENTRY BOXES

    ##FOR THE MODEL
    for k in self.mod_const:
      self.mod[k][3].set(self.mod_store[k][2])
    for k in self.mod_change:
      self.mod[k][3].set(self.mod_store[k][2])

    ##FOR THE DUST
    for k in self.dust:
      self.dust[k][3].set(self.dust_store[k][2])

    ##FOR THE DISK 
    for k in self.ss_const:
      self.ss[k][3].set(self.ss_store[k][2])
    for k in self.ss_change:
      self.ss[k][3].set(self.ss_store[k][2])

    ##FOR THE DISK GRID
    for k in self.ss_var_store:
        for i in range(4):
          self.ss_var[k][i+5].set(self.ss_var_store[k][i+1])

    ##FOR THE DISK OPTIONS

    for k in self.ss_options_store:
      self.ss_options[k][1].set(self.ss_options_store[k])

    ##AND GREY OUT BOXES AS APPROPRIATE

    a=self.check_disk_state(1)
    

   ##---------------------------------------------

  def start_loop(self):

    """

    This is the workhorse routine, which sets up the variables, runs
    the code, saves the output, and writes the log file.  

    """

    ##CHECK FOR APPROPRIATE INPUT, AND IF IT'S OKAY, START

    self.check_params()

    if(self.goodinput==1):

    ##SET UP VARIABLES FOR LOOPS AND CREATE THE ARRAY OF VALUES, IN LOG OR 
    ##LINEAR SPACE
 
      ##H0

      h0_start=float(self.ss_var_store['h0'][1])
      h0_end=float(self.ss_var_store['h0'][2])
      h0_nstep=int(self.ss_var_store['h0'][3])

      if(self.ss_var['h0'][8].get() == 0):
        h_array=numpy.linspace(h0_start,h0_end,num=h0_nstep,endpoint="True")
      if(self.ss_var['h0'][8].get() == 1):
        h_array=10**numpy.linspace(log10(h0_start),log10(h0_end),h0_nstep,endpoint="True")

      ##IF WE'RE VARYING RHO
      if(self.ss_options['fixrho'][1].get()<3):

        rho_start=float(self.ss_var_store['rho'][1])
        rho_end=float(self.ss_var_store['rho'][2])
        rho_nstep=int(self.ss_var_store['rho'][3])

        if(self.ss_var['rho'][8].get() == 0):
          rho_array=numpy.linspace(rho_start,rho_end,num=rho_nstep,endpoint="True")
        if(self.ss_var['rho'][8].get() == 1):
          rho_array=10**numpy.linspace(log10(rho_start),log10(rho_end),rho_nstep,endpoint="True")
      else:

        rho_array=[1]

      ##BETA

      beta_start=float(self.ss_var_store['beta'][1])
      beta_end=float(self.ss_var_store['beta'][2])
      beta_nstep=int(self.ss_var_store['beta'][3])
      beta_step=(beta_end-beta_start)/beta_nstep
      if(self.ss_var['beta'][8].get() == 0):
        beta_array=numpy.linspace(beta_start,beta_end,num=beta_nstep,endpoint="True")
      if(self.ss_var['beta'][8].get() == 1):
        beta_array=10**numpy.linspace(log10(beta_start),log10(beta_end),beta_nstep,endpoint="True")

    ##SET THE VARIABLES FOR THE MODEL

      pixsize=float(self.mod_store['pixscale'][2])
      dist=float(self.mod_store['dist'][2])
      resamp=float(self.mod_store['resamp'][2])
      maxsca=int(self.mod_store['maxsca'][2])
      dpath=float(self.mod_store['dpath'][2])
      ntot=int(float(self.mod_store['ntot'][2])*1e6)
      incl=float(self.mod_store['incl'][2])
      irange=float(self.mod_store['irange'][2])

    ##SET THE VARIABLES FOR THE DISK

      rstar=float(self.ss_store['r'][2])
      rmin=float(self.ss_store['minr'][2])
      maxr=float(self.ss_store['maxr'][2])


      if(float(self.ss_store['tmax'][2]) > 0):
        tmax=float(self.ss_store['tmax'][2])
        tau_critical=-1
      else:
        tau_critical=float(self.ss_store['tau_critical'][2])
        tmax=-1
      
    ##SET OPTIONS FOR H
      if(self.ss_options['fixh'][1].get()==1):
        h_for_arbitrary_r_in_AU   = -1.
        r_for_h_in_AU             =  1.
      else:
        r_for_h_in_AU             =  float(self.ss_options['h0_val'][1].get())
     
      rho_for_arbitrary_r       = -1.
      r_for_rho_in_AU           = -1.
      tau_for_arbitrary_theta   = -1.
      theta_for_arbitrary_tau   = -1.
      disk_mass                 = -1.

    ##SET OPTIONS FOR R

      #if(self.ss_options['fixrho'][1].get()==1):
      #  print ""

      if (self.ss_options['fixrho'][1].get()==2):
        r_for_rho_in_AU           =  float(self.ss_options['rho_val'][1].get())


      elif (self.ss_options['fixrho'][1].get()==3):
        theta_for_arbitrary_tau   = float(self.ss_options['ang_val'][1].get())
        tau_for_arbitrary_theta   = float(self.ss_options['tau_val'][1].get())

      elif (self.ss_options['fixrho'][1].get()==4):
        disk_mass=float(self.ss_options['mass_val'][1].get())
    
      disk_model=self.ss_store['disk_model'][2]

    ##SET THE VARIABLES FOR THE DUST

      dust_model=self.dust_store['dust_model'][2]
      wavelength=float(self.dust_store['wave'][2])
      kappa_ext=float(self.dust_store['kappa_ext'][2])

      ##THE ARRAY OF ANGLES

      self.angles=[12.9,31.4,41.3,49.4,56.5,63.2,69.5,75.5,81.4,87.1,incl]

    ##CALCULATE THE SIZE OF THE GRID BASED ON THE RESAMPLING FACTOR

      igrid=ceil(2*maxr/pixsize/dist/resamp)

    ##MAKE SURE THE NUMBER OF PIXELS IS ODD (FOR CENTRING)
      
      if(igrid%2==0):
        igrid=igrid+1
      
      self.mod_store['igrid'][2]=igrid

    ##PRINT THE START MESSAGE AND QUERY FOR START
      pad=float(self.mod['pixpad'][3].get())/100.
      pad=ceil(igrid*pad)
      newigrid=igrid+2*pad

      self.nrun=len(h_array)*len(rho_array)*len(beta_array)  ##NUMBER OF MODELS TO BE RUN

      message="" 
      message=message+"\n\n"
      message=message+"Will run "+str(self.nrun)+" models at "+str(newigrid) +" x "+ str(newigrid) +" pixels.\n"
      
      answer=tkMessageBox.askquestion("Parameters",message)

      ##IF TO CONTINUE

      if(answer == 'yes'):

        ##GET DIRECTORY NAME
        self.dirname=tkFileDialog.askdirectory()
        ##??? check ERRROR

        if (self.dirname != ""):
          if (not os.path.exists(self.dirname)): 
            os.makedirs(self.dirname)

        ##CREATE THE COUNTER WINDOW

          ##WRITE THE MODEL SUMMARY INFORMATION TO THE LOG FILE

          f=open(self.dirname+"/"+self.mod_store['logfile'][2],'w')
          f.write("-----RUNNING CODE-------\n\n")
          
          
          f.write("-----MODEL DIMENSIONS-------\n\n")
          
          sline="  Distance = "+str(dist)+" AU\n"
          f.write(sline)
          sline="  Pixel Scale = "+str(pixsize)+" arcsec\n"
          f.write(sline)
          sline="  Correponding grid is "+str(newigrid)+" x "+str(newigrid)+" pixels\n\n"
          f.write(sline)
          
          f.write("-----GRID PARAMETERS-------\n\n")
          
          
          sline="  Running a total of "+str(self.nrun)+" model(s)\n\n"
          f.write(sline)
          
          sline="  "+str(len(h_array))+" different values for h\n"
          f.write(sline)

          sline="          h values: "
          for hh in h_array:
            sline=sline+" "+str(hh)
          sline=sline+"\n\n"
          f.write(sline)
          
          if(self.ss_options_store['fixh']==2):
            sline="  "+"h values scaled from a radius of "+str(r_for_h_in_AU)+" AU\n"
          if(self.ss_options_store['fixh']==1):
            sline="  "+"h values scaled from a radius R*\n"
          f.write(sline)
          
          
          sline="  "+str(len(beta_array))+" different values for beta\n"
          f.write(sline)
          sline="          beta values: "
          for bb in beta_array:
            sline=sline+" "+str(bb)
          sline=sline+"\n\n"
          f.write(sline)

          
          if(self.ss_options_store['fixrho']==2):
            sline="   "+str(len(rho_array))+" different values for rho\n"
            f.write(sline)
            sline="          rho values: "
            for rr in rho_array:
              sline=sline+" "+str(rr)
            sline=sline+"\n\n"
            f.write(sline)

            sline="  "+"rho values scaled from "+str(r_for_rho_in_AU)+" AU\n"
            f.write(sline)
          
          if(self.ss_options_store['fixrho']==1):
            sline="   "+str(len(rho_array))+" different values for rho\n"
            f.write(sline)
            sline="          rho values: "
            for rr in rho_array:
              sline=sline+" "+str(rr)
            sline=sline+"\n\n"
            f.write(sline)

            sline="  "+"rho values scaled from R*\n"
            f.write(sline)
          
          f.write("\n")
          if(self.ss_options_store['fixrho']==3):
            sline="  Values for rho fixed by the condition\n"
            f.write(sline)
            sline="  tau= "+str(tau_for_arbitrary_theta)+" at an angle of "+str(theta_for_arbitrary_tau)+" degrees\n"
            f.write(sline)
          
          if(self.ss_options_store['fixrho']==4):
            sline="  Values for rho fixed by a total disk mass of "+str(disk_mass)+" solar masses\n"
            f.write(sline)

          f.write("\n")

          sline="---Model Parameters---\n\n"
          f.write(sline)
          sline="  Number of Photons = "+str( ntot  )+"\n"
          f.write(sline)
          sline="  User Specified Angle of Inclination = "+str( incl )+" degrees\n\n"
          f.write(sline)

          sline="---Dust Parameters---\n\n"
          f.write(sline)
          sline="  Dust Model = "+str( dust_model )+"\n"
          f.write(sline)
          sline="  Wavelength = "+str( wavelength )+" microns \n\n"
          f.write(sline)

          sline="---Fixed Disk Parameters---\n\n"
          f.write(sline)
          sline="  Disk Model = "+str( disk_model  )+"\n"
          f.write(sline)
          sline="  Stellar Radius = "+str( rstar  )+" solar radii\n"
          f.write(sline)
          sline="  Minimum Disk Radius = "+str( rmin  )+" AU\n"
          f.write(sline)
          sline="  Maximum Disk Radius = "+str( maxr  )+" AU\n\n"
          f.write(sline)

          sline="---Advanced  Parameters---\n\n"
          f.write(sline)
          sline="  kappa_ext = "+str( kappa_ext  )+"\n"
          f.write(sline)
          sline="  tau_critical = "+str( tau_critical  )+"\n"
          f.write(sline)
          sline="  tmax = "+str( tmax  )+"\n"
          f.write(sline)
          sline="  Range for User Specified Angle of Inclination = "+str( irange  )+"\n"
          f.write(sline)
          sline="  Maximum Number of Scatters = "+str( maxsca  )+"\n"
          f.write(sline)
          sline="  Step Size for Photon Path = "+str( dpath  )+"\n\n"
          f.write(sline)

          
          sline="---STARTING GRID OF MODELS---\n\n"

          ##UPDATE THE COUNTER WINDOW

          self.counterwin=counter_window(len(h_array)*len(rho_array)*len(beta_array))
          self.counterwin.pack(side=TOP,anchor=W,fill=X)
          self.counterwin.update_idletasks()
          self.counterwin.update()

          start_time=time.time()
          ii=0

          ##NOW CYCLE THROUGH THE LOOP IN H, RHO AND BETA

          ##KEEP TRACK OF THE DENSITY VALUES FOR SCALED DENSITY
          scaled_rho=zeros((len(h_array),len(beta_array)))
          nh=-1

          for h in h_array:
            nh=nh+1
            for rho in rho_array:
              nb=-1
          
              for beta in beta_array:
                nb=nb+1
                ii=ii+1

                if(self.ss_options['fixh'][1].get()==2):
                  h_for_arbitrary_r_in_AU=h
                if (self.ss_options['fixrho'][1].get()==2):
                  rho_for_arbitrary_r=rho

              ##UPDATE THE COUNTER WINDOW, AND AFTER THE FIRST MODEL ESTIMATE COMPLETION TIME

                #self.counterwin.num.set(ii)
                self.counterwin.numdone['text']=str(ii)
                if(ii > 1):
                  current_time=time.time()
                  time_left=(current_time-start_time)/(ii-1.)*(len(h_array)*len(rho_array)*len(beta_array)-ii+1.0)
                  #self.counterwin.time.set(str(int(time_left/60.)))
                  self.counterwin.timeleft['text']=str(int(time_left/60))
                self.update_idletasks()
                self.counterwin.update()

                alpha=beta+1
        

                ##WRITE THE SINGLE MODEL INFORMATION TO THE LOG FILE

                sline="Running Model "+str(ii)+" of "+str(self.nrun)+"\n"
                f.write(sline)
                
                sline="  alpha = "+str(  alpha )+"\n"
                f.write(sline)
                sline="  beta = "+str(  beta  )+"\n"
                f.write(sline)
                sline="  h = "+str(  h  )+"\n"
                f.write(sline)
                
                
                if(self.ss_options_store['fixrho']==1):
                  sline="  rho = "+str(rho)+"\n"
                  f.write(sline)
                if(self.ss_options_store['fixrho']==2):
                  sline="  rho = "+str(rho+for_arbitrary_r)+"\n"
                  f.write(sline)
                

              ##CALL THE SIMULATION CODE AND ASSIGN THE VARIABLES
                                     
                params = MC_simulation(dust_model,wavelength,
                                       disk_model,
                                       alpha,beta,
                                       rmin,maxr, #Rmin_in_AU,Rmax_in_AU,
                                       ntot,
                                       igrid,
                                       incl,
                                       irange,
                                       dpath,
                                       
                                       rho,
                                       disk_mass,
                                       rho_for_arbitrary_r,
                                       r_for_rho_in_AU,
                                       tau_for_arbitrary_theta,
                                       theta_for_arbitrary_tau,
                                        
                                       h,
                                       h_for_arbitrary_r_in_AU, 
                                       r_for_h_in_AU,
                                       
                                       kappa_ext,
                                       rstar,
                                      
                                       -1, # half_image_size_in_AU
                                       maxsca,
                                       tau_critical,
                                       tmax)                  
                params.run_simulation()

                self.twod=params.SS_disk.twod_density_distribution
                self.theta_array,self.r_array,self.z_array=params.SS_disk.get_lines_for_tau_from_star(1.0,d_theta=0.1)
                
                new_rho=params.SS_disk.rho_0
                if(self.ss_options_store['fixrho']==3 or (self.ss_options_store['fixrho']==4)):
                  rhoflag=1
                else:
                  rhoflag=0


              ##DETERMINE OUTPUT FILE NAME: USED SCALED RHO IF APPROPRIATE. 
              
                if(rhoflag==1):
                  outputfile="run_"+str("%5.3e" % h)+"_"+str("%5.3e" % new_rho)+"_"+str("%5.3e" % beta)
                  print nh,nb
                  scaled_rho[nh][nb]=new_rho
                else:
                  outputfile="run_"+str("%5.3e" % h)+"_"+str("%5.3e" % rho)+"_"+str("%5.3e" % beta)

                sline=outputfile
                f.write(sline)

              ##CALCULATE PADDING AROUND THE IMAGE AND ASSIGN THE DATA

                pad=float(self.mod['pixpad'][3].get())/100.
                pad=ceil(igrid*pad)

                data={ "I":numpy.zeros((11,igrid+2*pad,igrid+2*pad)),
                       "Q":numpy.zeros((11,igrid+2*pad,igrid+2*pad)),
                       "U":numpy.zeros((11,igrid+2*pad,igrid+2*pad)),
                       "V":numpy.zeros((11,igrid+2*pad,igrid+2*pad)),
                       "PI":numpy.zeros((11,igrid+2*pad,igrid+2*pad)),
                       }
              
                data["I"][:,pad-1:pad+igrid-1,pad-1:pad+igrid-1]=params.image_I/float(resamp)**2
                data["Q"][:,pad-1:pad+igrid-1,pad-1:pad+igrid-1]=params.image_Q/float(resamp)**2
                data["U"][:,pad-1:pad+igrid-1,pad-1:pad+igrid-1]=params.image_U/float(resamp)**2
                data["V"][:,pad-1:pad+igrid-1,pad-1:pad+igrid-1]=params.image_V/float(resamp)**2
              
 
              ##WRITE TO THE OUTPUT FILE: USED SCALED RHO IF APPROPRIATE. 
                if(self.savedat.get() == 1):
                  
                  if(rhoflag==1):
                    self.writeout_textfile(self.dirname+'/'+outputfile+".dat",data,self.angles,new_rho,h,beta)
                  else:
                    self.writeout_textfile(self.dirname+'/'+outputfile+".dat",data,self.angles,rho,h,beta)
                  
              ##DEFINE THE PI VECTOR AND SAVE THE DATA CUBE AND 2D PROFILE TO A PKL FILE
              
                data['PI']=sqrt(data["Q"]**2+data["U"]**2)
              
                myObject = [data,self.twod,self.r_array,self.z_array]
                f1 = open(self.dirname+'/'+outputfile+'.pkl','w')
                pickle.dump(myObject,f1)
                f1.close()
             

          ##AND FINISH UP THE LOG FILE
          f.write("\n -------FINISHED SIMULATIONS SUCCESSFULLY---------\n")
          f.close


          
        ##SAVE THE GENERAL RUN INFORMATION TO A SINGLE PICKLE FILE
          igrid=igrid+2*pad

          self.mod_store['igrid'][2]=igrid

          if(rhoflag==0):
            scaled_rho=[]


          myObject= [h_array,rho_array,beta_array,self.mod_store,self.dust_store,self.ss_store,self.angles,self.ss_options_store,self.ss_var_store,scaled_rho]
          f1 = open(self.dirname+'/'+'run_info.pkl','w')
          pickle.dump(myObject,f1)
          f1.close()

          self.counterwin.destroy()

          tkMessageBox.showinfo("","Finished Running Models")

   ##----------------------------------------------------------
      
  def writeout_textfile(self,outputfile,data,angles,rho,h,beta):
 
   """

   Routine to write to a text file. 


   """


    ##THIS LETS US USE UTF-8 CHARACTERS, INCLUDING RHO, ETC.

   f=codecs.open(outputfile,'w',encoding="utf8")

   ##PPRINT OUT PARAMETERS FROM MODEL

   f.write("##  PARAMETERS ##")
   f.write("## ")

   f.write("## outputfile = "+outputfile+"\n")              

   ##MODEL PARAMETERS
   for k in self.mod_store:
     sline="## "+self.mod_store[k][0]+" = "+str(self.mod_store[k][2])+" "+self.mod_store[k][1]+"\n"
     f.write(sline)

   ##DUST PARAMETERS
   for k in self.dust_store:
     sline="## "+self.dust_store[k][0]+" = "+str(self.dust_store[k][2])+" "+self.dust_store[k][1]+"\n"
     f.write(sline)

   ##DISK PARAMETERS
   for k in self.ss_store:
     sline="## "+self.ss_store[k][0]+" = "+str(self.ss_store[k][2])+" "+self.ss_store[k][1]+"\n"
     f.write(sline)

   ##H OPTIONS
   if(self.ss_options_store['fixh']==0):
     sline="  "+"h values scaled from a radius R*"+"\n"
     f.write(sline)
   if(self.ss_options_store['fixh']==1):
     sline="  "+"h values scaled from a radius R = "+self.ss_options_store['h0_val']+"\n"
     f.write(sline)

   ##RHO OPTIONS
   if(self.ss_options_store['fixrho']==0):
     sline="  "+"rho values scaled from a radius R*"+"\n"
     f.write(sline)
   if(self.ss_options_store['fixrho']==1):
     sline="  "+"rho values scaled from a radius R = "+self.ss_options_store['rho_val']+"\n"
     f.write(sline)
   if(self.ss_options_store['fixrho']==3):
     sline="  Values for rho fixed by the condition"+"\n"
     f.write(sline)
     f.write(sline)
     sline="  tau= "+self.ss_options_store['tau_val']+" at an angle of "+self.ss_options_store['ang_val']+" degrees"+"\n"
     f.write(sline)
     f.write(sline)
   if(self.ss_options_store['fixrho']==3):
     sline="  Values for rho fixed by the condition"+"\n"
     f.write(sline)
     sline="  disk_mass= "+self.ss_options_store['mass_val']+" solar masses"+"\n"
     f.write(sline)

   f.write("##\n")

   f.write("## rho_0 =  "+str(rho)+"\n")                      
   f.write("## alpha =  "+str(beta+1)+"\n")                      
   f.write("## beta =  "+str(beta)+"\n")                      
   f.write("## h_0 =  "+str(h)+"\n")                        
   
   f.write("##\n")

   ##PRINT OUT PARAMETERS FROM DUST
   
   f.write("##     Angle     X    Y     I             Q            U               V            \n")

   for k in range(11):
       for i in range(int(self.mod_store['igrid'][2])):
           for j in range(int(self.mod_store['igrid'][2])):
             strn='%12.4f %5d %5d %12.4e %12.4e %12.4e %12.4e \n' % (angles[k], i, j,
                                                                            data["I"][k][j][i],
                                                                            data["Q"][k][j][i],
                                                                            data["U"][k][j][i],
                                                                            data["V"][k][j][i])

             f.write(strn)

   f.close



a=MC_parameter_window() 
a.pack() 
a.mainloop() 


