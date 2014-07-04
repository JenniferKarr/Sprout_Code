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

   This is the GUI for viewing models created by run_models.py

   Version 1.0 (release) Last Updated May 25 2013



"""

from Tkinter import * 
import tkMessageBox 
import tkFileDialog 
from pylab import * 
import scipy
import scipy.stats
import numpy
import numpy.ma as ma
from numpy import outer
import pickle
import pyfits
import Pmw
import os
import string
import sys


import convolve_psf as c
import congrid
from textfile_window import *

##------------------------------------------------------------------------

class cmap_window(Frame):

 def __init__(self,master=None):

     """

     This routine creates a popup window to select colour maps, and will display
     a view of the available colour maps. 

     """


     Frame.__init__(self,master)

     ##CURRENT SPECTRAL MAP AND LIST OF OPTIONS
     self.cmap='spectral'  
     
     ##LIST OF COLOUR MAPS

     self.cmaps=['spectral','autumn','bone','cool','copper','flag','gray','hot','hsv','jet','pink','prism','spring','summer','winter']
     self.cmaps_num=range(15)

     self.col=IntVar() ##VARIABE TO HOLD COLOUR CHOICE
     
     ##CREATE WINDOW AND GET LIST

     top=Toplevel(self) 

     fm=Frame(top,relief=RAISED, bd=1)
     mb2=Button(fm,text='Done',command=top.destroy).pack(side=LEFT,anchor=W)
     mb2=Button(fm,text='Show Colours',command=self.cmap_show).pack(side=LEFT,anchor=W)
     fm.pack()

     fm=Frame(top,relief=RAISED, bd=1)
     i=0
     for i in range(len(self.cmaps)):
       col=self.cmaps[i]
       mb2=Radiobutton(fm,text=col,variable=self.col,value=i,command=self.cmap_get).pack(side=TOP,anchor=W)
     fm.pack()
     
     ##SET TO DEFAULT VALUE

     self.col.set(0)

##------------------------------------------------------------------------

 def cmap_get(self):

     ##SET THE COLOUR MAP


     self.cmap=self.cmaps[self.col.get()]

##------------------------------------------------------------------------

 def cmap_show(self):

     """
     A routine to show the availible colour maps in graphical form. 

     """

     ##SHOW THE AVAILABLE COLOUR MAPS GRAPHICALLY

   ##CREATE A NEW FIGURE AND CLEAR IT

     fig=figure(3)
     clf()
        
     ##CREATE AN ARRAY THAT VARIES FROM ZERO TO ONE IN THE X AXIS,
     ##AND IS CONSTANT IN THE Y AXIS

     a=outer(ones(10),arange(0,1,0.01))

     i=0
     ##FOR EACH COLOUR 
     for col in self.cmaps:
         subplot(len(self.cmaps),1,i)

         ##THIS PLOTS A BAND STRETCHING OVER THE WHOLE COLOUR MAP

         plot=imshow(a,aspect='auto',cmap=get_cmap(col),origin="lower")

         ##NO TICKMARKS
         plot.axes.get_xaxis().set_ticks([])
         plot.axes.get_yaxis().set_ticks([])

         ##AND LABELS IT WITH THE COLOUR NAME
         annotate(col, (0,1),color='white')
         i=i+1
     ##AND RESET THE FIGURE NUMBER
     fig=figure(1)

     if(string.lower(sys.platform) != "darwin"):
       show()


##------------------------------------------------------------------------

class parameter_view(Frame): 
    def __init__(self,master=None): 


        """
        This is the class for creating a GUI to run a grid of models. 


        """


        Frame.__init__(self,master) 
     
        ##INITIALIZE VARIABLES
        
        ##GENERAL VARIABLES

        self.convol=IntVar()               ##FLAG FOR CONVOLUTION WITH PSF
        self.masksize=StringVar()          ##MASK SIZE IN PIXELS
        self.loaded=0                      ##HAS THE DATA BEEN LOADED
        self.image_log=IntVar()            ##FLAG FOR LOG IMAGE
        self.contour_loaded=0              ##HAS TEH CONTOUR BEEN LOADED
        self.which_pol=StringVar()         ##VARIABLE FOR WHICH DISPLAY MODE TO USE


        ##CONTOUR VARIABLES

        self.contour=IntVar()              ##FLAG FOR PLOTTING CONTOURS
        self.contour_num=IntVar()          ##NUMBER OF CONTOURS
        self.contour_low=DoubleVar()       ##LOW CONTOUR VALUE
        self.contour_high=DoubleVar()      ##HIGH CONTOUR VALUE
        self.contour_log=IntVar()          ##FLAG FOR LOG CONTOURS
        self.contour_cust=IntVar()         ##FLAG FOR CUSTOM CONTOURS
        self.contour_levels=StringVar()    ##VALUES OF CONTOURS
        self.contour_colour=StringVar()    ##COLOUR OF CONTOURS
        self.rotvar=DoubleVar()            ##ROTATION OF THE CONTOURS

        ##H, BETA, RHO AND ANGLE VARIABES FOR SCROLLING

        self.h_val=StringVar()             ##CURRENT H VALUE
        self.beta_val=StringVar()          ##CURRENT BETA VALUE
        self.rho_val=StringVar()           ##CURRENT RHO VALUE
        self.angle_val=DoubleVar()         ##CURRENT ANGLE VALUE

        ##VARIABLES FOR BATCH PLOT

        self.p_angle_val1=StringVar()      ##FIRST ANGLE FOR BATCH PLOT
        self.p_angle_val2=StringVar()      ##SECOND ANGLE FOR BATCH PLOT
        self.p_h_val1=StringVar()          ##FIRST H FOR BATCH PLOT
        self.p_h_val2=StringVar()          ##SECOND H FOR BATCH PLOT
        self.p_rho_val1=StringVar()        ##FIRST RHO FOR BATCH PLOT
        self.p_rho_val2=StringVar()        ##SECOND RHO FOR BATCH PLOT
        self.p_beta_val1=StringVar()       ##FIRST BETA FOR BATCH PLOT
        self.p_beta_val2=StringVar()       ##SECOND BETA FOR BATCH PLOT
        self.choice1=IntVar()              ##FIRST VARIABLE FOR BATCH PLOT
        self.choice2=IntVar()              ##SECOND VARIABLE FOR BATCH PLOT


        ##POLARIZATION VECTOR VARIABLES

        self.showvect=IntVar()            ##FLAG FOR SHOWING POLARIZATION VECTORS
        self.vectsamp=StringVar()         ##VECTOR RESAMPLING FACTOR
        self.vectscale=StringVar()        ##VECTOR SCALING FACTOR
        self.vectcolour=StringVar()       ##VECTOR COLOUR


        ##PROFILE VARIABLES

        self.prof=IntVar()                ##FLAG TO SELECT SLICE OR WEDGE
        self.whichax=IntVar()             ##VARIABLE TO SELECT WHICH AXIS TO TAKE THE PROFILE
        self.profangle=StringVar()        ##ANGLE IN DEGREES FOR THE WEDGE PROFILE
        self.profwidth=StringVar()        ##WIDTH IN PIXELS FOR THE SLICE PROFILE
        self.profunit=StringVar()         ##VARIABILE TO SET THE UNITS (DEGREES OR PIXELS) 
        self.whichlab=StringVar()         ##VARIABILE TO SET THE LABEL (WIDTH/ANGLE)
        self.proflogx=IntVar()            ##FLAG FOR LOG SCALE IN X AXIS
        self.proflogy=IntVar()            ##FLAG FOR LOG SCALE IN Y AXIS
        self.datascale=StringVar()        ##FACTOR BY WHICH TO SCALE THE DATA

        
        self.xo=0  ##VARIABLES FOR ROTATION
        self.yo=0

        ##LABEL AND TYPES FOR PLOTS

        self.plottype={"PI Image":0,"I Image":1,"Q Image":2,"U Image":3,"V Image":4,"Density Distribution":5,
                       "Radial Profile (I)":6,"Radial Profile (PI)":7}

        self.plotlab={"PI Image":"PI","I Image":"I","Q Image":"Q","U Image":"U","V Image":"V","Density Distribution":"rho",
                       "Radial Profile (I)":"profi","Radial Profile (PI)":"profpi"}

        


        ##PSF VARIABLES

        self.whichpsf=IntVar()               ##WHICH TYPE OF PSF
        self.psfunit=StringVar()             ##UNITS LABEL FOR PSF INFO
        self.psflab=StringVar()              ##LABABEL FOR PSF INFO
        self.psfinfo=StringVar()             ##USER DEFINED INFO FOR PSF
        self.goodpsf=0                       ##FLAG FOR GOOD PSF
        self.oldpsfinfo="a"                  ##VARIABLE TO CHECK PSF VALUES
        self.convflag=1                      ##FLAG TO CHOOSE CONVOLUTION METHOD. SET HERE.

        ##----------------------------------------------------------------------

        ##SET UP THE GUI

        ##TOP BAR - QUITTING AND LOADING OPTIONS

        fm=Frame(self, relief=RAISED, bd=1) 
        mb=Button(fm,text='Quit',command=self.quit).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Help',command=self.show_help).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Load Simulation',command=self.load_simulation_dir).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Load Observations',command=self.load_contour_data).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Save Current File',command=self.save_current).pack(side=LEFT,anchor=W)
        fm.pack(side=TOP,anchor=W,fill=X) 
 

        ##OPTIONS FOR IMAGE DISPLAY TYPE

        fm=Frame(self, relief=RAISED, bd=1) 
        mb=Label(fm,text="Display Format ",width=12).pack(side=LEFT,anchor=W)
        opt=OptionMenu(fm,self.which_pol,"PI Image","I Image","Q Image","U Image","V Image","Density Distribution",
                       "Radial Profile (I)","Radial Profile (PI)",command=self.update_plot)
        opt.pack(side=LEFT,anchor=W)
        mb=Checkbutton(fm,text="Show Data Contour?  ",variable=self.contour, command=self.start_loop).pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Contour Rotation (deg)').pack(side=LEFT,anchor=W)
        mb=Entry(fm,textvariable=self.rotvar,width=5).pack(side=LEFT,anchor=W)
        fm.pack(side=TOP,anchor=W,fill=X) 


       ##SECOND BAR - GLOBAL OPTIONS LINE 1 (PSF)

        fm=Frame(self, relief=RAISED, bd=1) 
        mb=Checkbutton(fm,text="Convolve with PSF?",variable=self.convol, command=self.start_loop).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="Gaussian PSF",variable=self.whichpsf,command=self.set_convol_parms,value=2).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="File PSF",variable=self.whichpsf,command=self.set_convol_parms,value=1).pack(side=LEFT,anchor=W)
        self.convtype=Label(fm,  textvariable=self.psflab, width=8)
        self.convtype.pack(side=LEFT,anchor=W)
        mb=Entry(fm,textvariable=self.psfinfo,width=20).pack(side=LEFT,anchor=W)
        self.convtype=Label(fm,  textvariable=self.psfunit, width=8)
        self.convtype.pack(side=LEFT,anchor=W)
        fm.pack(side=TOP,anchor=W,fill=X) 

        self.psfunit.set('pixels')
        self.whichpsf.set(2)
        self.psfinfo.set('8')
        self.psflab.set('FWHM=')

        ##SECOND BAR - GLOBAL OPTIONS LINE 1 (VECTORS, MASK)
        fm=Frame(self, relief=RAISED, bd=1) 
        mb=Checkbutton(fm, text="Show Polarization Vectors?  ",variable=self.showvect, command=self.start_loop).pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Mask Size (pixels)').pack(side=LEFT,anchor=W)
        mb=Entry(fm,textvariable=self.masksize,width=3).pack(side=LEFT,anchor=W)
        fm.pack(side=TOP,anchor=W,fill=X) 

        ##THIRD BAR - GLOBAL OPTIONS LINE 2 (CONTOURS)

        ##SET ROTATION ANGLE AND MASK SIZE DEFAULTS
       
        fm=Frame(self) 
        mb=Button(fm,text='Plot/Refresh Plot',command=self.start_loop).pack(side=LEFT,anchor=W)
        fm.pack(side=BOTTOM) 

        ##SET DEFAULT VALUES
        self.which_pol.set('PI Image')
        self.rotvar.set(0)
        self.masksize.set('10')


        ##----------------------------------------------------------------------

        ##USE PMW NOTEBOOK TABS TO ORGANIZE THE OPTIONS.  FOUR TABS, ONE FOR 
        ##THE SCROLLING OPTIONS, ONE FOR IMAGE DISPLAY OPTIONS, ONE FOR
        ##SNAPSHOT MODE, AND ONE FOR PROFILE DISPLAY OPTIONS

        notebook=Pmw.NoteBook(self)  ##INITIALIZE THE NOTEBOOK

        ##FIRST NOTEBOOK WINDOW; SCROLLING OPTIONS

        page1 = notebook.add('Model Controls')
        notebook.tab('Model Controls').focus_set()
 
        ##TITLE

        fm=Frame(page1)
        Label(fm, text="Model Controls", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()


        ##SCALEBAR FOR ANGLE. INITIALIZE TO ZERO, TO BE UPDATED WHEN MODEL IS LOADED

        fm=Frame(page1, relief=RAISED, bd=1) 
        self.ang_label = Label(fm, text="Angle", width=8).pack(side=LEFT,anchor=W)
        self.ang_value = Label(fm, text="", width=10)
        self.ang_value.pack(side=LEFT,anchor=W)
        fm.pack()
        self.angle_scale = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_plot, length=300)
        self.angle_scale.pack()
        fm.pack()

        #SCALEBAR FOR H_0. INITIALIZE TO ZERO, TO BE UPDATED WHEN MODEL IS LOADED

        fm=Frame(page1, relief=RAISED, bd=1) 
        self.h_label = Label(fm, text="h_0", width=8).pack(side=LEFT,anchor=W)
        self.h_value = Label(fm, text="", width=10)
        self.h_value.pack(side=LEFT,anchor=W)
        fm.pack()
        self.h_scale = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_plot, length=300)
        self.h_scale.pack()
        fm.pack()
        
        ##SCALEBAR FOR RHO. INITIALIZE TO ZERO, TO BE UPDATED WHEN MODEL IS LOADED
        
        fm=Frame(page1, relief=RAISED, bd=1) 
        self.rho_label = Label(fm, text=u"\u03c1_0", width=8).pack(side=LEFT,anchor=W)
        self.rho_value = Label(fm, text="", width=10)
        self.rho_value.pack(side=LEFT,anchor=W)
        fm.pack()
        self.rho_scale = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_plot, length=300)
        self.rho_scale.pack()
        fm.pack()
        
        ##SCALEBAR FOR BETA. INITIALIZE TO ZERO, TO BE UPDATED WHEN MODEL IS LOADED
        
        fm=Frame(page1, relief=RAISED, bd=1) 
        self.beta_label = Label(fm, text=u"\u03b2_0", width=8).pack(side=LEFT,anchor=W)
        self.beta_value = Label(fm, text="", width=10)
        self.beta_value.pack(side=LEFT,anchor=W)
        fm.pack()
        self.beta_scale = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_plot, length=300)
        self.beta_scale.pack()
        fm.pack()

        ##OPEN THE DISPLAY WINDOW AND SET UP THE MOUSE BINDINGS FOR ROTATING THE CONTOUR INTERACTIVELY

        ##CREATE THE FIGURE
        self.fig = plt.figure()
        ax=self.fig.add_subplot(111)
        

        ##ON CLICKING THE MOUSE BUTTON 
        def onclick(event):

            ##READ THE X AND Y POSITIONS OF THE START POINT
            self.xo=event.xdata
            self.yo=event.ydata

        ##ON RELEASING THE MOUSE BUTTON
        def offclick(event):
            ##READ THE X AND Y POSITIONS OF THE END POINT
            self.xn=event.xdata
            self.yn=event.ydata

            ##SET THE ROTATION FLAG TO TELL START_LOOP TO UPDATE
            self.rotit=1

            ##CALCULATE THE CHANGE IN ANGLE, AND ADD IT TO ROTVAR

            angle1=180*numpy.arctan(self.yo/self.xo)/numpy.pi
            angle2=180*numpy.arctan(self.yn/self.xn)/numpy.pi
            self.rotangle=self.rotangle+angle2-angle1
            self.rotvar.set(self.rotangle)

            ##UPDATE THE PLOT
            self.start_loop()
        
        ##THIS CONNECTS THE TWO ABOVE ROUTINES TO THE MOUSE ACTIONS    

        cid = self.fig.canvas.mpl_connect('button_press_event', onclick)
        cid = self.fig.canvas.mpl_connect('button_release_event', offclick)
    

        ##----------------------------------------------------------------------

        ##SECOND NOTEBOOK TAB - IMAGE OPTIONS

        page2 = notebook.add('Image Options')
        notebook.tab('Image Options').focus_set()

        ##TITLE

        fm=Frame(page2)
        Label(fm, text="Image Options", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()

        ##SCALEBAR FOR SCALE.  INITIALIZE TO FULL RANGE

        bound=Frame(page2,relief=RAISED, bd=1)

        fm=Frame(bound)
        self.low_label = Label(fm, text="Lower Limit", width=10).pack(side=LEFT,anchor=W)
        self.low_lim = Scale(fm, from_=0, to=100, 
                               orient=HORIZONTAL, command=self.update_plot, length=300)
        self.low_lim.pack()
        self.low_lim.set(0)
        fm.pack()

        fm=Frame(bound)
        self.high_label = Label(fm, text="Upper Limit", width=10).pack(side=LEFT,anchor=W)
        self.high_lim = Scale(fm, from_=0, to=100, 
                               orient=HORIZONTAL, command=self.update_plot, length=300)
        self.high_lim.pack()
        self.high_lim.set(100)
        fm.pack()

        fm=Frame(bound)
        self.padd = Label(fm, text=" ", width=10).pack(side=LEFT,anchor=W)
        fm.pack()

        ##IMAGE SCALE AND COLOUR OPTIONS

        fm=Frame(bound)
        mb=Checkbutton(fm,text="Logarithmic Image Scale ",variable=self.image_log,command=self.start_loop).pack(side=LEFT,anchor=W)
        mb=Button(fm,text="Change Colour Map",command=self.get_new_cmap).pack(side=LEFT,anchor=W)
        fm.pack()
        bound.pack(fill=X)

        self.image_log.set(0)

        ##CONTOUR OPTIONS 

        bound=Frame(page2,relief=RAISED, bd=1)
        fm=Frame(bound)
        Label(fm, text="Contour Options", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()
        ##LOG OR LINEAR CONTOURS, CONTOUR COLOURS
        fm=Frame(bound)
        mb=Checkbutton(fm,text="Logarithmic Contours     ",variable=self.contour_log).pack(side=LEFT,anchor=W)
        self.contour_label = Label(fm, text="  Contour Colour", width=12).pack(side=LEFT,anchor=W)
        self.contour_label = Entry(fm, textvariable=self.contour_colour, width=10).pack(side=LEFT,anchor=W)
        fm.pack()
        bound.pack(fill=X)
        self.contour_colour.set('white')

        ##GENERATE CUSTOM CONTOURS, VIEW DATA STATISTICS
        fm=Frame(bound)
        mb=Checkbutton(fm,text="Custom Contours",variable=self.contour_cust,command=self.start_loop).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Generate',command=self.generate_contours).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Data Statistics',command=self.calculate_data_stats).pack(side=LEFT,anchor=W)
        fm.pack()
        self.contour_cust.set(0)

        ##CHOOSE PARAMETERS FOR CONTOURS
        fm=Frame(bound)
        self.contour_label = Label(fm, text="Number of Contours", width=20).pack(side=LEFT,anchor=W)
        self.contour_label = Entry(fm, textvariable=self.contour_num, width=10).pack(side=LEFT,anchor=W)
        self.contour_label = Label(fm, text="Contour Limits", width=15).pack(side=LEFT,anchor=W)
        self.contour_label = Entry(fm, textvariable=self.contour_low, width=10).pack(side=LEFT,anchor=W)
        self.contour_label = Entry(fm, textvariable=self.contour_high, width=10).pack(side=LEFT,anchor=W)
        fm.pack()


        ##SHOW OR EDIT THE CONTOURS MANUALLY
        fm=Frame(bound)
        self.contour_label = Label(fm, text="Levels", width=15).pack(side=LEFT,anchor=W)
        self.contour_label = Entry(fm, textvariable=self.contour_levels, width=40).pack(side=LEFT,anchor=W)

        ##A BIT FO PADDING BEFORE FOR VISUAL CLARITY
        fm.pack()
        fm=Frame(bound)
        self.padd = Label(fm, text=" ", width=10).pack(side=LEFT,anchor=W)
        fm.pack()
        bound.pack(fill=X)

        ##POLARIZATION VECTOR OPTIONS
        fm=Frame(page2)
        Label(fm, text="Polarization Options", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()

        fm=Frame(page2, relief=RAISED, bd=1) 
        mb=Label(fm,text='Vector Resampling').pack(side=LEFT,anchor=W)
        mb=Entry(fm, textvariable=self.vectsamp, width=3).pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Length Scale').pack(side=LEFT,anchor=W)
        mb=Entry(fm, textvariable=self.vectscale, width=3).pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Vector Colour').pack(side=LEFT,anchor=W)
        mb=Entry(fm, textvariable=self.vectcolour, width=8).pack(side=LEFT,anchor=W)
        fm.pack()
        self.vectsamp.set('1')
        self.vectscale.set('1')
        self.vectcolour.set('white')

        ##----------------------------------------------------------------------

        ##THIRD PAGE FOR THE NOTEBOOK; SNAPSHOT VIEW

        page3 = notebook.add('Snapshot View')
        notebook.tab('Snapshot View').focus_set()

        ##TITLE
        fm=Frame(page3)
        Label(fm, text="Snapshot View", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()

        ##BUTTONS FOR MAKING THE IMAGES
        fm=Frame(page3)
        mb=Button(fm,text='Make Figs',command=self.make_all_figs).pack(side=LEFT,anchor=W)
        mb=Button(fm,text='Save All Figs',command=self.cycle_all).pack(side=LEFT,anchor=W)
        fm.pack()

        ##ANGLE VARIABLE SCALES. SET TO ZERO, AS THE DEFAULTS WILL BE RESET WHEN THE 
        ##SIMULATION IS LOADED

        fm=Frame(page3, relief=RAISED, bd=1) 
        mb=Radiobutton(fm,variable=self.choice1,value=1).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,variable=self.choice2,value=1).pack(side=LEFT,anchor=W)
        self.p_ang_label = Label(fm, text="Angle Range", width=15).pack(side=LEFT,anchor=W)

        self.p_ang_value1 = Label(fm,text="",width=10)
        self.p_ang_value1.pack(side=LEFT,anchor=W)

        self.p_angle_scale1 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_angle_scale1.pack(side=LEFT,anchor=W)
        self.p_ang_value2 = Label(fm,text="",width=10)
        self.p_ang_value2.pack(side=LEFT,anchor=W)

        self.p_angle_scale2 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_angle_scale2.pack(side=LEFT,anchor=W)
        fm.pack()

        ##H0 VARIABLE SCALES SET TO ZERO, AS THE DEFAULTS WILL BE RESET WHEN THE 
        ##SIMULATION IS LOADED

        fm=Frame(page3, relief=RAISED, bd=1) 
        mb=Radiobutton(fm,variable=self.choice1,value=2).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,variable=self.choice2,value=2).pack(side=LEFT,anchor=W)
        self.p_h_label = Label(fm, text="h Range", width=15).pack(side=LEFT,anchor=W)
        self.p_h_value1 = Label(fm,text="",width=10)
        self.p_h_value1.pack(side=LEFT,anchor=W)

        self.p_h_scale1 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_h_scale1.pack(side=LEFT,anchor=W)
        self.p_h_value2 = Label(fm,text="",width=10)
        self.p_h_value2.pack(side=LEFT,anchor=W)

        self.p_h_scale2 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_h_scale2.pack(side=LEFT,anchor=W)
        fm.pack()

        ##RHO VARIABLE SCALES SET TO ZERO, AS THE DEFAULTS WILL BE RESET WHEN THE 
        ##SIMULATION IS LOADED

        fm=Frame(page3, relief=RAISED, bd=1) 
        mb=Radiobutton(fm,variable=self.choice1,value=3).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,variable=self.choice2,value=3).pack(side=LEFT,anchor=W)
        self.p_rho_label = Label(fm, text=u"\u03c1_0 Range", width=15).pack(side=LEFT,anchor=W)
        self.p_rho_value1 = Label(fm,text="",width=10)
        self.p_rho_value1.pack(side=LEFT,anchor=W)
        self.p_rho_scale1 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_rho_scale1.pack(side=LEFT,anchor=W)
        self.p_rho_value2 = Label(fm,text="",width=10)
        self.p_rho_value2.pack(side=LEFT,anchor=W)
        self.p_rho_scale2 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_rho_scale2.pack(side=LEFT,anchor=W)
        fm.pack()

        ##BETA VARIABLE SCALES SET TO ZERO, AS THE DEFAULTS WILL BE RESET WHEN THE 
        ##SIMULATION IS LOADED

        fm=Frame(page3, relief=RAISED, bd=1) 
        mb=Radiobutton(fm,variable=self.choice1,value=4).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,variable=self.choice2,value=4).pack(side=LEFT,anchor=W)
        self.p_beta_label = Label(fm, text=u"\u03b2_0 Range", width=15).pack(side=LEFT,anchor=W)
        self.p_beta_value1 = Label(fm,text="",width=10)
        self.p_beta_value1.pack(side=LEFT,anchor=W)
        self.p_beta_scale1 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_beta_scale1.pack(side=LEFT,anchor=W)
        self.p_beta_value2 = Label(fm,text="",width=10)
        self.p_beta_value2.pack(side=LEFT,anchor=W)
        self.p_beta_scale2 = Scale(fm, from_=0, to=0, 
                        orient=HORIZONTAL, command=self.update_figs, length=150)
        self.p_beta_scale2.pack(side=LEFT,anchor=W)
        fm.pack()

        ##DEFAULT VALUES FOR THE RADIO BUTTONS

        self.choice1.set(1)
        self.choice2.set(2)
        
        ##----------------------------------------------------------------------

        ##FOURTH NOTEBOOK PAGE, PROFILE OPTIONS 

        page4 = notebook.add('Profile Options')
        notebook.tab('Profile Options').focus_set()

        ##TITLE
        fm=Frame(page4)
        Label(fm, text="Profile Options", width=50,background="honeydew3").pack(side=LEFT,anchor=W)
        fm.pack()

        ##SAVE FUNCTION
        fm=Frame(page4)
        mb=Button(fm,text='Save Current Profile',command=self.save_radial_prof).pack(side=LEFT,anchor=W)
        fm.pack()


        ##CHOOSE TYPE OF PROFILE

        fm=Frame(page4, relief=RAISED, bd=1) 
        mb = Label(fm, text="Profile Type", width=10).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="Wedge",variable=self.prof,command=self.set_prof_parms,value=1).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="Slice",variable=self.prof,command=self.set_prof_parms,value=2).pack(side=LEFT,anchor=W)        
        self.profangle_lab=Label(fm,  textvariable=self.whichlab, width=15)
        self.profangle_lab.pack(side=LEFT,anchor=W)
        self.profangle_ent=Entry(fm, textvariable=self.profangle, width=10).pack(side=LEFT,anchor=W)
        self.profangle_unit=Label(fm,  textvariable=self.profunit, width=10)
        self.profangle_unit.pack(side=LEFT,anchor=W)
        fm.pack()

        ##CHOOSE THE AXIS
        fm=Frame(page4, relief=RAISED, bd=1) 
        mb = Label(fm, text="Axis", width=10).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="+X",variable=self.whichax,command=self.set_prof_parms,value=1).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="+Y",variable=self.whichax,command=self.set_prof_parms,value=2).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="-X",variable=self.whichax,command=self.set_prof_parms,value=3).pack(side=LEFT,anchor=W)
        mb=Radiobutton(fm,text="-Y",variable=self.whichax,command=self.set_prof_parms,value=4).pack(side=LEFT,anchor=W)
        fm.pack()

        ##LOGARITHMIC OR LINEAR, SSCALING
        fm=Frame(page4, relief=RAISED, bd=1) 
        mb=Checkbutton(fm,text="Logarithmic Plot Scale (X)",variable=self.proflogx, command=self.set_prof_parms).pack(side=LEFT,anchor=W)
        mb=Checkbutton(fm,text="Logarithmic Plot Scale (Y)",variable=self.proflogy, command=self.set_prof_parms).pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Scale data by',width=15).pack(side=LEFT,anchor=W)
        Entry(fm, textvariable=self.datascale, width=10).pack(side=LEFT,anchor=W)
        fm.pack()


        self.modelplot=StringVar()
        self.dataplot=StringVar()
        self.modelcol=StringVar()
        self.datacol=StringVar()

        ##PLOTTING OPTIONS

        fm=Frame(page4, relief=RAISED, bd=1) 
        mb = Label(fm, text="Plot Type (Model)", width=15).pack(side=LEFT,anchor=W)
        opt=OptionMenu(fm,self.modelplot,'Diamond','Circle','Square','Triangle','Plus','None',command=self.update_plot)
        opt.pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Plot Color (Model)',width=15).pack(side=LEFT,anchor=W)
        Entry(fm, textvariable=self.modelcol, width=10).pack(side=LEFT,anchor=W)
        fm.pack()

        fm=Frame(page4, relief=RAISED, bd=1) 
        mb = Label(fm, text="Plot Type (Data)", width=15).pack(side=LEFT,anchor=W)
        opt=OptionMenu(fm,self.dataplot,'Diamond','Circle','Square','Triangle','Plus','None',command=self.update_plot)
        opt.pack(side=LEFT,anchor=W)
        mb=Label(fm,text='Plot Color (Data)',width=15).pack(side=LEFT,anchor=W)
        Entry(fm, textvariable=self.datacol, width=10).pack(side=LEFT,anchor=W)
        fm.pack()



        ##SET DEFAULT VALUES FOR PROFILES

        self.datascale.set('1')
        self.whichlab.set('Angle of Wedge')
        self.prof.set(1)
        self.whichax.set(1)
        self.profangle.set('5')
        self.profunit.set('degrees')
        self.proflogx.set(1)
        self.proflogy.set(0)
        self.modelplot.set('None')
        self.dataplot.set('Diamond')
        self.modelcol.set('black')
        self.datacol.set('green')


        ##NOW PACK THE NOTEBOOK AND SET THE SIZE TO CONTAIN ALL THE CONTENTS

        notebook.pack(fill = 'both', expand = 1, padx = 10, pady = 10)
        notebook.setnaturalsize()

##------------------------------------------------------------------------
    def set_convol_parms(self):

      """
      Set the appropriate labels and defults 
        

      """


      ##SET THE PARAMETERS APPROPRIATELY

      convtype=self.whichpsf.get()
      if(convtype==1):
        self.psflab.set('PSF File=')
        self.psfunit.set('')
        self.psfinfo.set('psf_refstar_cropped.fits')
      if(convtype==2):
          self.psflab.set('FWHM=')
          self.psfunit.set('pixels')
          self.psfinfo.set('8')

##------------------------------------------------------------------------

    def show_disk_profile(self):

        """
        Displays a disk density profiles. 

        """

        ##DISPLAY A DISK DENSITY PROFILE 

        #clf()

        ##SET REASONABLE VALUES FOR MIN AND MAX OF DISPLAY

        mmin=self.twod[0,self.twod.shape[1]-1]
        mmax=self.twod[0,self.twod.shape[1]/10]
        try:
          self.cmap=get_cmap(self.cc.cmap)
        except:
            self.cmap=cm.spectral

        self.cmap.set_bad('black',1) ##MISSING OR BAD PIXELS ARE BLACK

        ##DO THE PLOT

        self.plotit=imshow(self.twod,origin='lower',interpolation='nearest',extent=[0,self.maxr,0,self.maxr],vmin=mmin,vmax=mmax,cmap=self.cmap)
        plot(self.r_array,self.z_array,'w-' )

        if(string.lower(sys.platform) != "darwin"):
          show()

     ##--------------------------------------------------------------------------------------------

    def generate_contours(self):

        """

        Create custom contour levels and format as a string, based on user  input


        """

        ##ROUTINE TO CREATE CUSTOM CONTOURS

        contourstring=""    ##STRING CONTAINING THE CONTOUR VALUES, SEPARATED BY COMMAS
        
        ##mmin          MINIMUM VALUE FOR CONTOURS
        ##mmax          MAXIMUM VALUE FOR CONTOURS
        ##ncont         NUMBER OF CONTOURS
        ##levels        ARRAY OF CONTOUR LEVELS

        if self.contour.get() == 1:
            ##GET MIN/MAX/NUMBER
            mmin=self.contour_low.get()
            mmax=self.contour_high.get()
            ncont=self.contour_num.get()
            ##CHECK FOR VALID NUMBERS AND RESET IF INPUT IS INVALID
            if (mmin == mmax or mmin > mmax or ncont<=0):
                tkMessageBox.showwarning("Warning", "Incorrect Contour Parameters. \nValues must be positive, minimum value must be less than maximum")
                levels=0
                contourstring=""
        
            else:
                ##LINEAR CONTOURS
                if(self.contour_log.get() == 0):
                    levels=numpy.linspace(mmin,mmax,ncont)
                    print numpy.linspace(mmin,mmax,ncont)
                    #print "A",levels
                ##LOGARITHMIC CONTOURS
                else:

                    if(mmin > 0):
                        levels=10**numpy.linspace(start=log10(mmin),stop=log10(mmax),num=ncont,endpoint="True")
                        numpy.linspace(start=log10(mmin),stop=log10(mmax),num=ncont,endpoint="True")
                        print "B",levels
                    #GET LOGARITHMIC INTERVALS FOR VALUES < 0
                    if(mmin <=0):
                        lrange=abs(mmax-mmin)
                        mmin1=lrange/20.
                        mmax1=lrange
                        levels=10**numpy.linspace(start=log10(mmin1),stop=log10(mmax1),num=ncont,endpoint="True")-mmin1+mmin
                        print numpy.linspace(start=log10(mmin1),stop=log10(mmax1),num=ncont,endpoint="True")-mmin1+mmin
                        print "C",levels
                print "AAA", self.contour_log.get(),levels
                ##CREATE THE STRING

                for val in levels:
                    print "VAL",val
                    contourstring=contourstring+str("%5.3e" % val)+','

                ##STRIP OF THE FINAL COMMA
                contourstring=contourstring[:-1]
                print contourstring

           ##write THE STRING IN THE WINDOW
            self.contour_levels.set(contourstring)



##------------------------------------------------------------------------

    def make_all_figs(self):

        """
        
        Generates the snapshot mode images. 

        """

        ##h          ARRAY OF H VALUES FOR THE PLOT
        ##rho        ARRAY OF RHO VALUES FOR THE PLOT
        ##beta       ARRAY OF BETA VALUES FOR THE PLOT
        ##angle_ind  ARRAY OF ANGLES  FOR THE PLOT
        ##ncols      NUMBER OF COLUMNS IN THE PLOT
        ##xlabel     STRING FOR THE X AXIS LABEL
        ##nrows      NUMBER OF ROWS IN THE PLOT
        ##ylabel     STRING FOR THE Y AXIS LABEL
        ##sgn        SET TO -1 OR +1 DEPENDING ON WHETHER THE FIRST INDEX IS >= THAN THE SECOND FOR A VARIABLE
        ##ii         COUNTER FOR THE # OF PLOTS DONE SO FAR
        ##i,j,k,l    LOOP VARIABLES
        
        ##THIS ROUTINE GENERATE THE SNAPSHOT MODE IMAGE GRID.  
  
        ##CLOSE THE WINDOW IF ONE IS ALREADY OPEN

        close(4)
        

        ##HAS THE DATA BEEN LOADED? IF NOT, EXIT
        if self.loaded == 0:
            tkMessageBox.showwarning("", "Load Simulation Data")
            return
 
       ##MAKE SURE THAT THE THERE ARE TWO UNIQUE IMAGES CHOSEN

        if(self.choice1.get() == self.choice2.get()):
            tkMessageBox.showwarning("Warning", "Choose Two Different Variables")
            return


        ##SET THE VALUE OF EACH VARIABLE TO THE FIRST VALUE FROM THE SLIDER FOR EACH VARIABLE
        ##THIS WILL BE THE VALUE FOR THE UNSELECTED PARAMETERS

        angle_ind=[int(self.p_angle_scale1.get())]
        h=[self.h_array[int(self.p_h_scale1.get())]]
        rho=[self.rho_array[self.p_rho_scale1.get()]]
        beta=[self.beta_array[self.p_beta_scale1.get()]]
  

        ##NOW GET THE FIRST AND SECOND CHOICES FOR THE VARIABLES, CALCULATE THE NUMBER OF IMAGES TO PLOT
        ##AND SET THE AXIS LABELS


        ##SGN WILL BE POSITIVE IF THE FIRST INDEX IS SMALLER THAN THE LAST AND NEGATIVE OTHERWISE (SETS THE DIRECTION OF THE 
        ##PLOTS CORRECTLY. IF THE FIRST INDEX=SECOND INDEX, DEFAULTS TO 1
        ##THEN CALCULATE THE ARRAY OF VALUES FOR EACH VARIABLE, AND USE SGN TO GET THE ORDER RIGHT. 
        ##CALCULATE THE NUMBER OF COLUMNS AND THE LABEL FOR THE AXIS


        ##THIS CODE GETS THE DIRECTION OF THE VARIABLES TO PLOT, AND SETS THE ARRAYS OF VALUES
        ##IT'S THE SAME BASIC CODE, REPEATED FOR DIFFERENT FIRST VARIABLES AND SECOND VARIABLES.

        ##FOR THE FIRST VARIABLE

        if(self.choice1.get()==1):
          try:
            sgn=abs(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))/(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))
          except:
            sgn=1
          angle_ind=int(self.p_angle_scale1.get())+sgn*arange(abs(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))+1)
          ncols=len(angle_ind)
          if(sgn > 0):
            xlabel="Angle="+str(self.angles[min(angle_ind)])+" to "+str(self.angles[max(angle_ind)])
          if(sgn < 0):
            xlabel="Angle="+str(self.angles[max(angle_ind)])+" to "+str(self.angles[min(angle_ind)])

        if(self.choice1.get()==2):
          try:
            sgn=abs(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))/(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))
          except:
            sgn=1
          h=self.h_array[int(self.p_h_scale1.get())+sgn*arange(abs(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))+1)]
          ncols=len(h)
          if(sgn > 0):
            xlabel="h="+str("%5.2f" % min(h))+" to "+str("%5.2f" % max(h))
          if(sgn < 0):
            xlabel="h="+str("%5.2f" % max(h))+" to "+str("%5.2f" % min(h))

        if(self.choice1.get()==3):
          try:
            sgn=abs(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))/(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))
          except:
            sgn=1
          rho=self.rho_array[int(self.p_rho_scale1.get())+sgn*arange(abs(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))+1)]
          ncols=len(rho)
          if(sgn > 0):
            xlabel=u"\u03c1_0="+str("%5.2e" % min(rho))+" to "+str("%5.2e" % max(rho))
          if(sgn < 0):
            xlabel=u"\u03c1_0="+str("%5.2e" % max(rho))+" to "+str("%5.2e" % min(rho))

        if(self.choice1.get()==4):
          try:
            sgn=abs(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))/(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))
          except:
            sgn=1
          beta=self.beta_array[int(self.p_beta_scale1.get())+sgn*arange(abs(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))+1)]
          ncols=len(beta)
          if(sgn > 0):
            xlabel=u"\u03b2_0="+str("%5.2f" % min(beta))+" to "+str("%5.2f" % max(beta))
          if(sgn < 0):
            xlabel=u"\u03b2_0="+str("%5.2f" % max(beta))+" to "+str("%5.2f" % min(beta))

        ##FOR THE SECOND VARIABLE DO THE SAME THING, BUT FOR NROWS AND YLABEL

        if(self.choice2.get()==1):
          try:
            sgn=abs(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))/(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))
          except:
            sgn=1
          angle_ind=int(self.p_angle_scale1.get())+sgn*arange(abs(int(self.p_angle_scale2.get())-int(self.p_angle_scale1.get()))+1)
          angle_ind=angle_ind[::-1]
          nrows=len(angle_ind)
          if(sgn > 0):
            ylabel="Angle="+str(self.angles[min(angle_ind)])+" to "+str(self.angles[max(angle_ind)])
          if(sgn < 0):
            ylabel="Angle="+str(self.angles[max(angle_ind)])+" to "+str(self.angles[min(angle_ind)])

        if(self.choice2.get()==2):
          try:
            sgn=abs(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))/(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))
          except:
            sgn=1
          h=self.h_array[int(self.p_h_scale1.get())+sgn*arange(abs(int(self.p_h_scale2.get())-int(self.p_h_scale1.get()))+1)]
          h=h[::-1]
          nrows=len(h)
          if(sgn > 0):
            ylabel="h="+str("%5.2f" % min(h))+" to "+str("%5.2f" % max(h))
          if(sgn < 0):
            ylabel="h="+str("%5.2f" % max(h))+" to "+str("%5.2f" % min(h))

        if(self.choice2.get()==3):
          try:
            sgn=abs(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))/(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))
          except:
            sgn=1
          rho=self.rho_array[int(self.p_rho_scale1.get())+sgn*arange(abs(int(self.p_rho_scale2.get())-int(self.p_rho_scale1.get()))+1)]
          rho=rho[::-1]
          nrows=len(rho)
          if(sgn > 0):
            ylabel=u"\u03c1_0="+str("%5.2e" % min(rho))+" to "+str("%5.2e" % max(rho))
          if(sgn < 0):
            ylabel=u"\u03c1_0="+str("%5.2e" % max(rho))+" to "+str("%5.2e" % min(rho))

        if(self.choice2.get()==4):
          try:
            sgn=abs(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))/(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))
          except:
            sgn=1
          beta=self.beta_array[int(self.p_beta_scale1.get())+sgn*arange(abs(int(self.p_beta_scale2.get())-int(self.p_beta_scale1.get()))+1)]
          beta=beta[::-1]
          nrows=len(beta)
          if(sgn > 0):
            ylabel=u"\u03b2_0="+str("%5.2f" % min(beta))+" to "+str("%5.2f" % max(beta))
          if(sgn < 0):
            ylabel=u"\u03b2_0="+str("%5.2f" % max(beta))+" to "+str("%5.2f" % min(beta))


        ##SET THE SIZE OF THE WINDOW

        self.setsize(nrows,ncols)
  
        ##TWO VERSIONS OF THE LOOP, TO SORT OUT THE X AND Y AXES IN THE RIGHT ORDER, AND THE RIGHT VALUES. 

        clf()
        if(self.choice1.get() < self.choice2.get()):
                                
            ii=0
                  ##GENERATE THE PLOTS
                              
            for l in range(len(beta)):
                for k in range(len(rho)):
                    for j in range(len(h)):
                        for i in range(len(angle_ind)):
                            
                            ii=ii+1
                            ax=subplot(nrows,ncols,ii)
                            if(self.rhoflag==0):
                                self.plot_file(h[j],rho[k],beta[l],angle_ind[i],1,ax)  
                            else:
                                self.plot_file(h[j],self.scaled_rho[j][l],beta[l],angle_ind[i],1,ax)  

            ##LABEL THE AXES
            text(-1.*(ncols-1)/2.+0.1,0.1,xlabel,
                 transform = ax.transAxes,color='white')

            text(-1.*(ncols-1)+0.1,0.5*(nrows+1),ylabel,
                 transform = ax.transAxes,color='white',rotation=90)


        if(self.choice1.get() > self.choice2.get()):
          
            
            ii=0
                  ##GENERATE THE PLOTS
  
            
            for i in range(len(angle_ind)):
                for j in range(len(h)):
                    for k in range(len(rho)):
                        for l in range(len(beta)):

                            ii=ii+1
                            ax=subplot(nrows,ncols,ii)
                            if(self.rhoflag==0):
                                self.plot_file(h[j],rho[k],beta[l],angle_ind[i],1,ax)  
                            else:
                                self.plot_file(h[j],self.scaled_rho[j][l],beta[l],angle_ind[i],1,ax)  

                         
            ##LABEL THE AXES
            text(-1.*(ncols-1)/2.+0.1,0.1,xlabel,
                 transform = ax.transAxes,color='white')

            text(-1.*(ncols-1)+0.1,0.5*(nrows+1),ylabel,
                 transform = ax.transAxes,color='white',rotation=90)


        ##RESET THE FIGURE
                    
        figure(1)
##------------------------------------------------------------------------

    def setsize(self,nrows,ncols):

        """
        A simple estimation of the size of a window for the batch plots, to keep it to a printable size.

        """

        ##CALCULATE THE X-Y SIZE OF THE FIGURE FOR BATCH PLOTS, AND
        ##SET THE PLOTTING PARAMETERS. THE LARGEST AXIS WILL HAVE A
        ##SIZE OF 8 INCHES AT 75 DPI.
        
        ##sz1      SIZE OF HORIZONTAL AXIS
        ##sz2      SIZE OF VERTICAL AXIS

        if (nrows <= ncols):
            sz1=8
            sz2=8.*float(nrows)/float(ncols)
        if(nrows > ncols):
            sz1=8.*float(ncols)/float(nrows)
            sz2=8
        figure(4,figsize=(sz1,sz2))
        subplots_adjust(wspace=0.0,hspace=0.0,left=0,right=1,top=1,bottom=0)
        clf()


##------------------------------------------------------------------------

    def update_plot(self,value):

        ##SIMPLE WRAPPER TO CALL START_LOOP

        if self.loaded == 1:
            self.start_loop()

##------------------------------------------------------------------------
    def cycle_all(self):

        """
        
        Cycle through all the models and save an image for each, for the given viewing mode. 
        The file name is the same as the pkl file, but with the image name as prefix, and the angle
        appended to the name. In addition, a final plot of the observational data in the same form is produced. 


        """

        ##THIS ROUTINE CYCLES THROUGH ALL THE MODELS, AND GENERATES AND SAVES AN IMAGE FOR EACH. THE FILE NAME
        ##IS THE SAME AS THE PKL FILE, BUT WITH IMAGE AS PREFIX, AND THE ANGLE APPENDED TO THE NAME. THE FINAL IMAGE
        ##PLOTS THE DATA AS AN IMAGE IN THE SAME FASHION AS THE REST OF THE PLOTS. 


        ##h            CURRENT VALUE OF H IN LOOP
        ##beta         CURRENT VALUE OF H IN LOOP
        ##rho          CURRENT VALUE OF H IN LOOP
        ##angle_ind    CURRENT VALUE OF H IN LOOP

        ax=subplot(1,1,1)

        hn=-1
        for h in self.h_array:
            hn=hn+1
            bn=-1
            for beta in self.beta_array:
                bn=bn+1
                for rho in self.rho_array:
                    for angle_ind in range(11):
                        clf()
                        ax=subplot(1,1,1)
                        if(self.rhoflag==0):
                            self.plot_file(h,rho,beta,angle_ind,0,ax)
                            plotfile="image_"+str("%5.3e" % h)+"_"+str("%5.3e" % rho)+"_"+str("%5.3e" % beta)+str("%2i" % self.angles[angle_ind])+self.plotlab[self.which_pol.get()]+".jpg"

                        else:
                            self.plot_file(h,self.scaled_rho[hn][bn],beta,angle_ind,0,ax)

                            plotfile="image_"+str("%5.3e" % h)+"_"+str("%5.3e" % self.scaled_rho[hn][bn])+"_"+str("%5.3e" % beta)+str("%2i" % self.angles[angle_ind])+self.plotlab[self.which_pol.get()]+".jpg"

                        savefig(self.dirname+"/"+plotfile)
                        
        if (self.contour_loaded==1):
            imshow(self.cont_mask,interpolation="nearest",origin='lower',cmap=self.cmap,extent=[numpy.min(self.xaxis),numpy.max(self.xaxis),numpy.min(self.yaxis),numpy.max(self.yaxis)])

            colorbar()

            savefig("self.dirname+"/"+data.jpg")

        if(string.lower(sys.platform) != "darwin"):
          show()



##------------------------------------------------------------------------

    def update_figs(self,val):

        """
        Gets the values from the GUI for batch plots and sets the
        window labels, before running the plotting routine. This is called
        everytime the sliders are moved or buttons selected in snapshot mode.


        """

        ##h_var1         FORMATTED LABEL FOR FIRST H SLIDER
        ##rho_var1       FORMATTED LABEL FOR FIRST H SLIDER
        ##beta_var1      FORMATTED LABEL FOR FIRST H SLIDER
        ##ang_var1       FORMATTED LABEL FOR FIRST H SLIDER
        ##h_var2         FORMATTED LABEL FOR SECOND H SLIDER
        ##rho_var2       FORMATTED LABEL FOR SECOND H SLIDER
        ##beta_var2      FORMATTED LABEL FOR SECOND H SLIDER
        ##ang_var2       FORMATTED LABEL FOR SECOND H SLIDER


        ##MAKE SURE ITS LOADED.

        if self.loaded != 1:
          return

        if(self.choice1.get()== self.choice2.get()==1):
          tkMessageBox.showwarning("", "Choose 2 different variables for Snapshot Mode")
          return

        ##READ THE VALUES FROM THE SCALE BARS IN THE SNAPSHOT MODE, FORMAT THE VALUES
        ##PROPERLY, AND SET THE LABELS FOR THE SLIDERS 

        h_var1=str("%8.6f" % self.h_array[self.p_h_scale1.get()])
        beta_var1=str("%4.2f" % self.beta_array[self.p_beta_scale1.get()])
        angle_var1=str(self.angles[self.p_angle_scale1.get()])

        h_var2=str("%8.6f" % self.h_array[self.p_h_scale2.get()])
        beta_var2=str("%4.2f" % self.beta_array[self.p_beta_scale2.get()])
        angle_var2=str(self.angles[self.p_angle_scale2.get()])

        if(self.rhoflag==0):
            rho_var1=str("%5.3e"  % self.rho_array[self.p_rho_scale1.get()])
            rho_var2=str("%5.3e"  % self.rho_array[self.p_rho_scale2.get()])
        else:
            rho_var1==""
            rho_var2==""
            
        self.p_ang_value1["text"]=angle_var1
        self.p_h_value1["text"]=h_var1
        self.p_rho_value1["text"]=rho_var1

        ##FOR THE NON CHOSEN VALUES IN SNAPSHOT MODE, ONLY DISPLAY THE 
        ##VALUE IN THE FIRST SLIDER

        self.p_ang_value2['text']=""
        self.p_h_value2['text']=""	
        self.p_beta_value2['text']=""
        self.p_rho_value2['text']=""

        if(self.choice1.get()==1 or self.choice2.get()==1):
          self.p_ang_value2["text"]=angle_var2
        if(self.choice1.get()==2 or self.choice2.get()==2):
          self.p_h_value2["text"]=h_var2
        if(self.choice1.get()==3 or self.choice2.get()==3):
          self.p_beta_value2["text"]=beta_var2
        if(self.choice1.get()==3 or self.choice2.get()==3):
          self.p_rho_value2["text"]=rho_var2



##------------------------------------------------------------------------

    def start_loop(self):

        """

        This routine is called everytime the scale bars are moved, a button is clicked
        or the plot is refreshed, for normal plots. It reads the values for the scale bars,
        set the labels, and calls the plotting routine.

        Error checking for values is done here. 

        """


        ##GET THE DESIRED IMAGE VIEW

        self.pol=self.plottype[self.which_pol.get()]

        ###ERROR MESSAGES

        ##CHECK TO SEE IF THE IMAGE IS LOADED

        if(self.loaded == 0):
            #NO ERROR MESSAGE HERE BECAUSE OF REPETITION
            return

        ##CHECK FOR NUMERIC VALUES/POSITIVE/INTEGERS AS APPROPRIATE

        ##MASKSIZE MUST BE > 0
        try:
          radval=float(self.masksize.get())
        except:
          tkMessageBox.showwarning("", "Mask Size Must Be a Number >= 0")
          return

        if(radval < 0):
          tkMessageBox.showwarning("", "Mask Size Must Be a Number >= 0")
          return

        ##VECTOR RESAMPLING MUST BE A POSITIVE INTEGER

        try:
          z=int(self.vectsamp.get())
        except: 
          tkMessageBox.showwarning("", "Vector Resampling Must Be a Positive Integer")
          return
        if(z <= 0):
          tkMessageBox.showwarning("", "Vector Resampling Must Be a Positive Integer")
          return

        try:
          z=float(self.vectscale.get())
        except: 
          tkMessageBox.showwarning("", "Vector Rescaling Must Be a Positive Float")
          return
        if(z <= 0):
          tkMessageBox.showwarning("", "Vector Resampling Must Be a Positive Float")
          return

        
        ##PROFILE WIDTH MUST BE A POSITIVE NUMBER
        try:
          z=float(self.profangle.get())
        except:
          tkMessageBox.showwarning("", "Profile Width Must Be a Positive Number")
          return

        if(z <=0):
          tkMessageBox.showwarning("", "Profile Width Must Be a Positive Number")
          return

        ##SCALEING MUST BE A POSITIVE NUMBER

        try:
          z=float(self.datascale.get())
        except:
          tkMessageBox.showwarning("", "Scale Value for Profile Must be a Positive Number")
          return

        if(z <=0):
          tkMessageBox.showwarning("", "Scale Value for Profile Must be a Positive Number")
          return
  
        
        ##READ VALUES AND SET LABLES ON THE SCROLLBARS

        h=self.h_array[self.h_scale.get()]
        beta=self.beta_array[self.beta_scale.get()]
        angle_ind=self.angle_scale.get()
        self.angle_ind=angle_ind

        if(self.rhoflag==0):
            rho=self.rho_array[self.rho_scale.get()]
        else:
            rho=self.scaled_rho[self.h_scale.get()][self.beta_scale.get()]

        #print "AA",self.rhoflag,self.scaled_rho[self.h_scale.get()][self.beta_scale.get()],rho

        self.ang_value["text"]=str(self.angles[angle_ind])
        self.h_value["text"]=str("%8.6f" % h)
        self.beta_value["text"]=str("%4.2f" % beta)
        
        if(self.rhoflag==0):
            self.rho_value["text"]=str("%5.3e"  % rho)
        else:
            self.rho_value["text"]=""


        ##CLEAR THE WINDOW, CREATE A SUBPLOT, AND CALL THE PLOTTING ROUTINE

        clf()
        ax=subplot(1,1,1)

        self.plot_file(h,rho,beta,angle_ind,0,ax)

##------------------------------------------------------------------------

    def plot_file(self,h,rho,beta,angle_ind,batch,ax):

        """

        This is the workhorse routine for plotting the images. There are two modes, 
        one for snapshot mode, and one for single images. 

        """


        if(string.lower(sys.platform) != "darwin"):  ##MAC SPECIFIC
          ion()

        ##GET THE FILE NAME AND LOAD THE FILE

        self.outputfile="run_"+str("%5.3e" % h)+"_"+str("%5.3e" % rho)+"_"+str("%5.3e" % beta)

        f = open(self.dirname+"/"+self.outputfile+'.pkl','r')
        myObject=pickle.load(f)
        f.close()

        ##TITLES FOR TWO DIFFERENT CASES [INDEPENDENT OR SCALED RHO]

        titlehead=self.which_pol.get()


        

        title=titlehead+"  h0="+str("%5.2e" % h)+"  "+u"\u03b2_0="+str("%5.2e" % beta)+"  "+u"\u03c1_0="+str("%5.2e" % rho)+"   angle="+str("%2i" % self.angles[angle_ind])


        ##SET THE DATA TYPES

        alldata=myObject[0]
        self.alldata=alldata
        self.twod=myObject[1]
        self.r_array=myObject[2]
        self.z_array=myObject[3]

        ##DENSITY DISTRIBUTION CASE - SET THE LABELS AND CALL THE PLOTTING ROUTINE

        if(self.pol==5):

          #clf()
          self.show_disk_profile()

          if(batch==0):
            suptitle(title,fontsize=12)
            colorbar()
            xlabel('AU')
            ylabel('AU')

          if(batch==1):
            ax.set_xticklabels([])
            ax.set_yticklabels([])

          return

        ##GET THE MASK SIZE, IF ANY

        try:
            self.radval=float(self.masksize.get())
        except:
            self.radval=0

        ##GET THE CONTOUR DATA AND REBIN FOR THE DATA SIZE


        ##CASE FOR DISPLAYING CONTOUR, AND CONTOUR DATA IS ALREADY LOADED

        if self.contour.get() == 1 and self.contour_loaded==1:
            
            ##ROTATE IF APPROPRIATE
            self.rotangle=self.rotvar.get()

            self.rotcontdata=scipy.ndimage.interpolation.rotate(self.contresamp,-self.rotangle,reshape=False)
            self.rotit=0

            ##RESAMPLE AND CROP TO IMAGE SIZE AND DISPLAY
                
            self.cropcontdata=self.rotcontdata[(self.npix_new-self.npix)/2:(self.npix_new+self.npix)/2,(self.npix_new-self.npix)/2:(self.npix_new+self.npix)/2]
            
            self.cont_mask=ma.array(data=self.cropcontdata,mask=[self.rad < self.radval])

        ##NOW GET THE CONVOLUTION KERNEL, IF NEEDED.  

        ##TWO SEPARATE CASES - ONE FOR READING THE PSF FROM A FITS FILE
        ##ONE FOR CALULATING A GAUSSIAN 
        
        ##IN BOTH CASES, ONLY RECALCULATE IF THE PARAMETERS HAVE CHANGED,
        ##EXIT IF THE READ/CALCULATE ROUTINES FAIL

        if(self.convol.get()==1):

            if(self.whichpsf.get()==1):  
                if(self.oldpsfinfo != self.psfinfo.get()):  
                    self.read_psf()
                elif(self.oldpsfinfo == self.psfinfo.get() and self.goodpsf==0): 
                    tkMessageBox.showwarning("Warning", "Error Reading PSF File")
                    return

            if(self.whichpsf.get()==2):  
                if(self.oldpsfinfo != self.psfinfo.get()):
                    self.calc_psf()
                elif(self.oldpsfinfo != self.psfinfo.get()):
                    tkMessageBox.showwarning("Warning", "Error Calculating PSF")
                    return

        ##IF THE PSF DIDN'T WORK
        if(self.goodpsf==0 and self.convol==1):
            return


        ##CASE FOR I PROFILE - GET THE INFORMATION, AND CALL THE PLOTTING ROUTINE

        if(self.pol==6):
          #clf()
          self.PI=alldata["I"][angle_ind]
          self.set_prof_parms()
          if self.convol.get() == 1:

            self.PI=c.convolve(self.PI,self.psf,self.convflag,data_pixel_binning=self.resamp)
          

          self.calculate_radial_prof()
          if(batch==0):
            suptitle(title,fontsize=12)
            xlabel('AU')
            ylabel('Flux')

          if(batch==1):
            ax.set_xticklabels([])
            ax.set_yticklabels([])


          return

        ##CASE FOR PI PROFILE - GET THE INFORMATION AND CALL THE PLOTTING ROUTINE

        if(self.pol==7):
          #clf()
          self.PI=alldata["PI"][angle_ind]
          self.set_prof_parms()
          if self.convol.get() == 1:
            self.PI=c.convolve(self.PI,self.psf,self.convflag,data_pixel_binning=self.resamp)
          self.calculate_radial_prof()
          #suptitle(title,fontsize=16)

          if(batch==0):
            suptitle(title,fontsize=12)
            xlabel('AU')
            ylabel('Flux')

          if(batch==1):
            ax.set_xticklabels([])
            ax.set_yticklabels([])

          return


        ##NOW THE CASE FOR POLARIZATION IMAGES
        ##GET THE CORRECT POLARIZATON AND ANGLE

        PI=alldata[self.polar[self.pol]][angle_ind]

        ##KLUDGE FOR BRIGHT CENTRAL POINT AND DISPLAY
        if(self.pol==1):
            mm=PI.max()
            ind=where(PI < mm)
            mm1=PI[ind].max()
            if(mm > mm1*10.):
                ind=where(PI == mm)
                PI[ind]=0

                    
        ##CONVOLVE WITH THE PSF IF NEEDED 
        if self.convol.get() == 1:
            PI=c.convolve(PI,self.psf,self.convflag,data_pixel_binning=self.resamp)
        self.PI=PI
  
        ##SET THE MASK

        if(self.radval != 0):
          PI_mask=ma.array(data=PI,mask=[self.rad < self.radval])
          a=ma.MaskedArray.ravel(PI_mask[where(PI_mask != 0)])

        else:
          PI_mask=PI
          a=numpy.ravel(PI_mask[where(PI_mask != 0)])

        ##GET THE DATA LIMITS AND CALCULATE THE MIN/MAX VALUES
        per1=self.low_lim.get()
        per2=self.high_lim.get()
        
        self.minn=scipy.stats.mstats.scoreatpercentile(a,float(per1))
        self.maxx=scipy.stats.mstats.scoreatpercentile(a,float(per2))
        

        ##THRESHOLD THE IMAGE
        show_image=numpy.clip(PI_mask,self.minn,self.maxx)

        ##CONVERT TO LOG SCALE IF APPROPRIATE 

        if self.image_log.get() == 1:
            show_image=numpy.log10(show_image-numpy.min(show_image)+(numpy.max(show_image)-numpy.min(show_image))*0.1)

        ##SET THE COLOUR MAP AND MASK COLOUR

        try:
          self.cmap=get_cmap(self.cc.cmap)
        except:
            self.cmap=cm.spectral

        self.cmap.set_bad('black',1)

        ##PLOT CURRENT IMAGE

        parmplot=imshow(transpose(show_image),interpolation="nearest",origin='lower',cmap=self.cmap,extent=[numpy.min(self.xaxis),numpy.max(self.xaxis),numpy.min(self.yaxis),numpy.max(self.yaxis)])


        ##IS IT A BATCH JOB, OR A SINGLE IMAGE - SET TITLES AND TICKS AS NEEDED. 
        if(batch==0):
            colorbar()
            suptitle(title,fontsize=12)
            xlabel('AU')
            ylabel('AU')

        if(batch==1):
            ax.set_xticklabels([])
            ax.set_yticklabels([])

        if(self.showvect.get()==1):
                ##PLOT POLARIZATION VECTORS IF NEEDED

          resamp=int(self.vectsamp.get())


          ###USE THE CONVOLVED IMAGES IF THIS OPTION CHOSEN

          if(self.convol.get()==1):

              ii=c.convolve(alldata['I'][angle_ind],self.psf  ,self.convflag,data_pixel_binning=self.resamp)
              qq=c.convolve(alldata['Q'][angle_ind] ,self.psf ,self.convflag,data_pixel_binning=self.resamp)
              uu=c.convolve(-alldata['U'][angle_ind],self.psf ,self.convflag,data_pixel_binning=self.resamp)

              pp=numpy.sqrt(qq**2+uu**2)/ii

              ##THIS SET THE VALUES TO ZERO FOR AREAS OF LOW SIGNAL
              ##TO NOISE. WE USE THE MINIMUM VALUE OF THE UN-CONVOLVED
              ##IMAGE TO SET A REASONABLE CUTOFF

              ##MAKE A MASKED I IMAGE
              PI=alldata["I"][angle_ind]
              if(self.radval != 0):
                  PI_mask=ma.array(data=PI,mask=[self.rad < self.radval])
                  a=ma.MaskedArray.ravel(PI_mask[where(PI_mask != 0)])

              else:
                  PI_mask=PI
                  a=numpy.ravel(PI_mask[where(PI_mask != 0)])

              mina=scipy.stats.mstats.scoreatpercentile(a,float(3))

              #ind=where(PI_mask != 0)
              #llim=abs([ind]).min()
              ind=where(PI_mask < mina)
              pp[ind]=0

          ##OR THE SIMPLER, UNCONVOLVED OPTION
          if(self.convol.get()==0):
              ii=alldata['I'][angle_ind]  
              qq=alldata['Q'][angle_ind]  
              uu=-alldata['U'][angle_ind] 

              pp=numpy.sqrt(qq**2+uu**2)/ii

          pa=0.5*numpy.arctan2(uu,qq)

          dx=pp*numpy.sin(pa)
          dy=pp*numpy.cos(pa)

          xx,yy=meshgrid(self.xaxis,self.yaxis)

          xnew=xx[::resamp,::resamp]
          ynew=yy[::resamp,::resamp]
          dxnew=dx[::resamp,::resamp]
          dynew=dy[::resamp,::resamp]

          vscale=0.05/float(self.vectscale.get())
          vwidth=float(self.vectscale.get())*0.5

          plot=quiver(ynew,xnew,dynew,dxnew,pivot='middle',units='x',scale=vscale,width=vwidth,color=self.vectcolour.get(),headwidth=0)
          qk = quiverkey(plot, 0.05,0.1, 1, "5 %",
                            labelpos='E',labelcolor=self.vectcolour.get())

          ##NOW PLOT CONTOURS
          ##CASE FOR CUSTOM CONTOURS - READ THEM 

        if self.contour.get() == 1 and self.contour_loaded==1:
                
            ##CASE FOR CUSTOM CONTOURS - READ THEM 
            if self.contour_cust.get() == 1:

                levels=[]
                contourstring=self.contour_levels.get()

                if(contourstring != ""):
                    levelsarray=contourstring.rsplit(',')

                    try:
                        for val in levelsarray:
                            levels.append(float(val))
                            self.contour_plot=contour(self.cropcontdata,levels,origin='lower',colors=self.contour_colour.get(),linewidths=2,extent=[numpy.min(self.xaxis),numpy.max(self.xaxis),numpy.min(self.yaxis),numpy.max(self.yaxis)])

                    except:
                        tkMessageBox.showwarning("Warning", "There appears to be a problem with the contours")


       
            ##DEFAULT CONTOURS
            else:
                col="white"
                if(self.contour_colour.get() != ""):
                    col=self.contour_colour.get()
                self.contour_plot=contour(self.cropcontdata,origin='lower',colors=col,linewidths=2,extent=[numpy.min(self.xaxis),numpy.max(self.xaxis),numpy.min(self.yaxis),numpy.max(self.yaxis)])

                ##NOW PUT THE LEVELS IN THE CONTOUR ENTRY UNDER IMAGE OPTIONS

                levels=self.contour_plot.levels
                contourstring=""
                ##CREATE THE STRING
                for val in levels:
                    contourstring=contourstring+str(val)+','

                ##STRIP OF THE FINAL COMMA
                contourstring=contourstring[:-1]

           ##write THE STRING IN THE WINDOW
            self.contour_levels.set(contourstring)

        if(string.lower(sys.platform) != "darwin"):
            show()
            
    ##--------------------------------------------------------------------------------------------

    def read_psf(self):

        """

        Read the PSF from a file. 

        """

        try:
            fname=self.psfinfo.get()
            f=pyfits.open(fname)
            psf=f[0].data
            f.close()

            ##CONGRID THE PSF TO MATCH THE DATA
            nx_psf,ny_psf       = psf.shape
            psf=congrid.congrid(psf,(int(nx_psf/self.resamp),int(ny_psf/self.resamp)),centre=True)

            self.psf=psf/sum(psf)
            self.goodpsf=1
            self.oldpsfinfo=self.psfinfo.get()

        except:
            tkMessageBox.showwarning("Warning", "Error Reading PSF")
            self.oldpsfinfo=self.psfinfo.get()
            self.goodpsf=0

    ##--------------------------------------------------------------------------------------------

    def calc_psf(self):

        """

        Calculate a Gaussian PSF.

        """

        try:


            #GET FWHM AND CONVERT SO SIGMA 
            fwhm=float(self.psfinfo.get())
            sigma=fwhm/sqrt(8*log(2)) 

            #SET SIZE OF BOX, MAKING SURE IT'S ODD

            boxsize=int(fwhm*4)
            if(boxsize % 2)==0:
                boxsize=boxsize+1
            midpoint=boxsize/2
            
            #NOW CALCULATE PSF AND NORMALIZE

            self.psf=zeros((boxsize,boxsize))

            for i in range((boxsize)):
                for j in range((boxsize)):
                    nn=(i-midpoint)
                    mm=(j-midpoint)
                    self.psf[i][j]=exp(-(nn**2+mm**2)/(2*sigma**2))

            self.psf/=sum(self.psf)
            self.oldpsfinfo=self.psfinfo.get()

            ##WE NOW HAVE A GOOD PSF
            self.goodpsf=1

        except:
            tkMessageBox.showwarning("Warning", "Error Calculating PSF")
            self.goodpsf=0
            self.oldpsfinfo=self.psfinfo.get()

    ##--------------------------------------------------------------------------------------------

    def save_current(self):

        """

        Save the current model and inclination to a FITS file and current view to a PS image. 

        """

      ##ASSIGN THE DATA

        self.data=self.alldata

        ##CONVOLVE THE IMAGES IF NECESSARY
        if(self.convol.get()==1):
 
            self.data['PI'][self.angle_ind]=c.convolve(self.data['PI'][self.angle_ind],self.psf,self.convflag,data_pixel_binning=self.resamp)
            self.data['I'][self.angle_ind]=c.convolve(self.data['I'][self.angle_ind],self.psf,self.convflag,data_pixel_binning=self.resamp)
            self.data['Q'][self.angle_ind]=c.convolve(self.data['Q'][self.angle_ind],self.psf,self.convflag,data_pixel_binning=self.resamp)
            self.data['U'][self.angle_ind]=c.convolve(self.data['U'][self.angle_ind],self.psf,self.convflag,data_pixel_binning=self.resamp)
            self.data['V'][self.angle_ind]=c.convolve(self.data['V'][self.angle_ind],self.psf,self.convflag,data_pixel_binning=self.resamp)

      ##MAKE THE OUTPUT FILE NAME

        sangle= "%i" % self.angle
        outfile = self.outputfile+"_"+sangle+".fits"

      ##GET CURRENT VARIABLES

        h=self.h_array[self.h_scale.get()]
        beta=self.beta_array[self.beta_scale.get()]
        
        if(self.rhoflag==0):
            rho=self.rho_array[self.rho_scale.get()]
        else:
            rho=self.scaled_rho[self.h_scale.get()][self.beta_scale.get()]

      ##CREATE APPROPRIATE HEADERS AND WRITE

        ext=pyfits.PrimaryHDU(data=None)

      ##PUT THE MODEL PARAMETERS AND COMMENTS IN THE MAIN HEADER

        ##MODEL PARAMETERS
        for k in self.mod:
          sline="## "+self.mod[k][0]+" = "+str(self.mod[k][2])+" "+self.mod[k][1]
          ext.header.add_comment(sline)
        
        ##DUST PARAMETERS
        for k in self.dust:
          sline="## "+self.dust[k][0]+" = "+str(self.dust[k][2])+" "+self.dust[k][1]
          #ext.header.add_comment(sline)
            
        ##DISK PARAMETERS
        for k in self.ss:
          sline="## "+self.ss[k][0]+" = "+str(self.ss[k][2])+" "+self.ss[k][1]
          ext.header.add_comment(sline)
        
        ##H OPTIONS
        if(self.ss_options['fixh']==0):
          sline="  "+"h values scaled from a radius R*"
          ext.header.add_comment(sline)
        if(self.ss_options['fixh']==1):
          sline="  "+"h values scaled from a radius R = "+self.ss_options['h0_val']
          ext.header.add_comment(sline)
        
        ##RHO OPTIONS
        if(self.ss_options['fixrho']==0):
          sline="  "+"rho values scaled from a radius R*"
          ext.header.add_comment(sline)
        if(self.ss_options['fixrho']==1):
          sline="  "+"rho values scaled from a radius R = "+self.ss_options['rho_val']
          ext.header.add_comment(sline)
        if(self.ss_options['fixrho']==3):
          sline="  Values for rho fixed by the condition"
          ext.header.add_comment(sline)
          sline="  tau= "+self.ss_options['tau_val']+" at an angle of "+self.ss_options['ang_val']+" degrees"
          ext.header.add_comment(sline)
        if(self.ss_options['fixrho']==3):
          sline="  Values for rho fixed by the condition"
          ext.header.add_comment(sline)
          sline="  disk_mass= "+self.ss_options['mass_val']+" solar masses"
          ext.header.add_comment(sline)
        
          sline="## rho_0 =  "+str(rho)   
          ext.header.add_comment(sline)                   
          sline="## alpha =  "+str(beta+1)  
          ext.header.add_comment(sline)                    
          sline="## beta =  "+str(beta)   
          ext.header.add_comment(sline)                   
          sline="## h_0 =  "+str(h)    
          ext.header.add_comment(sline)                    
        
        
        ext_sub={}
        extlist=[ext]  
        
        self.add_header(ext)
   
        for k,v in self.data.iteritems():
            print "K",k
            thisarray=transpose(self.data[k][self.angle_ind])
            ext_sub[k]=pyfits.ImageHDU(data=thisarray[:,:],name=k)
            self.add_header(ext_sub[k])
            extlist.append(ext_sub[k])

        thdulist=pyfits.HDUList(extlist)
        thdulist.writeto(outfile)

      ##NOW SAVE THE IMAGE

        savefig(self.dirname+"/"+self.outputfile+".ps")
        

    ##--------------------------------------------------------------------------------------------

    def add_header(self,ext):

        """
        Create a FITS header with the appropriate units. 

        """

        ##CREATES A FITS HEARDER FOR THE DATA

        ext.header.update("CRVAL1",0)
        ext.header.update("CRVAL2",0)
        ext.header.update("CRPIX1",(self.npix+1.)/2.)
        ext.header.update("CRPIX2",(self.npix+1.)/2.)
        ext.header.update("CRTYPE1",'X')
        ext.header.update("CRTYPE2",'Y')
        ext.header.update("CDELT1",1)
        ext.header.update("CDELT2",1)
        ext.header.update("BUNIT","Number of Photons")
        ext.header.update("CUNIT1",'AU')
        ext.header.update("CUNIT2",'AU')

    ##--------------------------------------------------------------------------------------------

    def load_contour_data(self):

        """
        
        Load Observational data from a file. 

        """

        ##ROUTINE TO LOAD THE CONTOUR DATA. 

        ##RESET THE VARIABLES. 
        self.contour_loaded=0
        self.rotangle=0

        ##MAKE SURE SIMULATION DATA IS LOADED (NEEDED FOR RESAMPLING). 

        if(self.loaded==0):
            tkMessageBox.showwarning("Warning", "Load Simulation Data First")
            return   

        ##OPEN FILE AND READ THE DATA
        self.contfile=tkFileDialog.askopenfilename()

        f=pyfits.open(self.contfile)

        self.contdata=f[0].data

        ##RESAMPLE THE DATA TO THE SIMULATION RESOLUTION
        self.npix_new=int(floor(self.contdata.shape[0]/self.resamp))
        self.npix_new=int(floor(self.contdata.shape[1]/self.resamp))

        self.contresamp=congrid.congrid(self.contdata,(self.npix_new,self.npix_new),
                                                 centre=True)

        f.close()
        ##SUCCESSFULLY LOADED

        self.contour_loaded=1

    ##--------------------------------------------------------------------------------------------

    def calculate_data_stats(self):

        """
        
        Calculate the data statistics and plot a simple histogram in a new window. 

        """

        if(self.contour_loaded==0):
            return

        fig=figure(5)

        subplot(1,1,1)

        hh=hist(self.cropcontdata,bins=30,range=[numpy.min(self.cropcontdata),numpy.max(self.cropcontdata)])
        figure(1)

    ##--------------------------------------------------------------------------------------------

    def get_new_cmap(self):

        """

        Get the new colour map
        
        """

        self.cc=cmap_window()
        a=self.start_loop()

    ##--------------------------------------------------------------------------------------------

    def load_simulation_dir(self):


        """

        Routine to load the main file for a model. Sets the values for the variable ranges, 
        configure the various menu items, labels and scrollbars. 

        """

        ##READ IN THE FILE FOR THE SIMULATION DATA, AND INITIALIZE THE WINDOWS

        ##RESET THE VARIABLES. 
        self.loaded=0
        self.rotit=0

        ##READ IN MODEL INFORMATION FROM RUN_INFO.PKL

        self.dirname=tkFileDialog.askdirectory()
        if(self.dirname==""):
            return

        try:
            f = open(self.dirname+'/run_info.pkl','r')
            myObject=pickle.load(f)
            f.close()

        except:
            tkMessageBox.showwarning("Warning", "This does not appear to be the correct type of directory!")
            return

        ##ASSIGN VALUES

        self.h_array=myObject[0]
        self.rho_array=myObject[1]
        self.beta_array=myObject[2]
        
        self.mod=myObject[3]
        self.dust=myObject[4]
        self.ss=myObject[5]

        self.angles=myObject[6]

        self.npix=int(self.mod['igrid'][2])
        self.dist=float(self.mod['dist'][2])
        self.pixscale=float(self.mod['pixscale'][2])
        self.resamp=float(self.mod['resamp'][2])
        self.angle=float(self.mod['incl'][2])


        ##NAMES OF POLARIZATION STATES
        self.polar=['PI','I','Q','U','V']

        ##CALCULATE THE NUMBER OF X AND Y PIXELS

        self.xaxis=(numpy.arange(self.npix)-self.npix/2)*self.resamp*self.dist*self.pixscale
        self.yaxis=(numpy.arange(self.npix)-self.npix/2)*self.resamp*self.dist*self.pixscale

        ##INITIALIZE MASK RADIUS ARRAY

        self.rad=numpy.zeros((self.npix,self.npix))
        for i in range(self.npix):
            for j in range(self.npix):
                self.rad[i][j]=sqrt((i-float(self.npix)/2.+1)**2+(j-float(self.npix)/2.+1)**2)


        #CONFIURE THE SCALE BARS TO THE CORRECT DATA.  

        self.angle_scale.config(from_=0,to=len(self.angles)-1)
        self.h_scale.config(from_=0,to=len(self.h_array)-1)
        self.rho_scale.config(from_=0,to=len(self.rho_array)-1)
        self.beta_scale.config(from_=0,to=len(self.beta_array)-1)

        self.scaled_rho=myObject[9]

        if(len(self.scaled_rho)>0):
            self.rhoflag=1
        else:
            self.rhoflag=0


        ##SET THE VALUES IN THE SCALE BARS

        self.angle_val.set(self.angles[0])
        h_var=str("%5.3f" %  self.h_array[0])
        self.h_val.set(h_var)
        beta_var=str("%5.3f" % self.beta_array[0])
        self.beta_val.set(beta_var)

        if(self.rhoflag==0):
            rho_var=str("%5.3e" % self.rho_array[0])
        else:
            rho_var=""

        self.rho_val.set(rho_var)

        #INITIALIZE THE SCALE BARS FOR SNAPSHOT MODE. 

        self.p_angle_scale1.config(from_=0,to=len(self.angles)-1)
        self.p_h_scale1.config(from_=0,to=len(self.h_array)-1)
        self.p_rho_scale1.config(from_=0,to=len(self.rho_array)-1)
        self.p_beta_scale1.config(from_=0,to=len(self.beta_array)-1)

        self.p_angle_scale2.config(from_=0,to=len(self.angles)-1)
        self.p_h_scale2.config(from_=0,to=len(self.h_array)-1)
        self.p_rho_scale2.config(from_=0,to=len(self.rho_array)-1)
        self.p_beta_scale2.config(from_=0,to=len(self.beta_array)-1)

        self.p_angle_val1.set(self.angles[0])
        self.p_rho_val1.set(rho_var)
        self.p_h_val1.set(h_var)
        self.p_beta_val1.set(beta_var)

        self.p_angle_val2.set(self.angles[0])
        self.p_h_val2.set(h_var)
        self.p_rho_val2.set(rho_var)
        self.p_beta_val2.set(beta_var)


        #GENERAL VARIABLES

        self.dpath=float(self.mod['dpath'][2])

        self.disk_model=self.ss['disk_model'][2]
        self.r=float(self.ss['r'][2])
        self.minr=float(self.ss['minr'][2])
        self.maxr=float(self.ss['maxr'][2])
        self.taucrit=float(self.ss['tau_critical'][2])
        self.kappa_ext=float(self.dust['kappa_ext'][2])
        self.dust_model=self.dust['dust_model'][2]
        self.wavelength=float(self.dust['wave'][2])
        

        self.ss_options=myObject[7]

        if(self.ss_options['fixh']==1):
            self.h_for_arbitrary_r_in_AU   = -1.
            self.r_for_h_in_AU             =  1.
        else:
            self.r_for_h_in_AU             =  self.ss_options['h0_val']
      


    ##VALUES OF RHO AT R


        if(self.ss_options['fixrho']==1):
          self.rho_for_arbitrary_r       = -1.
          self.r_for_rho_in_AU           = -1.
          self.tau_for_arbitrary_theta   = -1.
          self.theta_for_arbitrary_tau   = -1.
          self.disk_mass= -1.

        elif (self.ss_options['fixrho']==2):
          self.r_for_rho_in_AU           =  self.ss_options['rho_val']
          self.tau_for_arbitrary_theta   = -1.
          self.theta_for_arbitrary_tau   = -1.
          self.disk_mass= -1.

        elif (self.ss_options['fixrho']==3):
          self.rho_for_arbitrary_r       = -1
          self.theta_for_arbitrary_tau   = self.ss_options['ang_val']
          self.tau_for_arbitrary_theta   = self.ss_options['tau_val']
          self.r_for_rho_in_AU           = -1.
          self.disk_mass= -1.

        else:
          self.rho_for_arbitrary_r       = -1
          self.theta_for_arbitrary_tau   = -1
          self.tau_for_arbitrary_theta   = -1
          self.r_for_rho_in_AU           = -1.
          self.disk_mass=self.ss_options['mass_val']
        
        self.loaded=1
  ##---------------------------------------------------------------------

    def show_help(self): 


        """
        
        Setup the help window and display the text. 

        """

        self.helpwin=textfile_window("Help File","help_viewer.txt")
        self.helpwin.pack(side=TOP,anchor=W,fill=X)
        self.helpwin.update()

  ##---------------------------------------------------------------------
    def set_prof_parms(self):

      """

      Set the labels for the plot profile window, depending on the
      choice of type of profile.

      """

      proftype=self.prof.get()
      if(proftype==1):
        self.whichlab.set('Angle of Wedge')
        self.profunit.set('degrees')
        self.profangle.set('5')
      if(proftype==2):
        self.whichlab.set('Width of Slice')
        self.profunit.set('pixels')
        self.profangle.set('5')
      if(proftype==3):
        self.whichlab.set('N/A')
        self.profunit.set('')
        self.profangle.set('')

  ##---------------------------------------------------------------------

    def save_radial_prof(self):

      """

      Save radial profile to a text file. 
      

      """

      f = open(self.dirname+"/"+self.outputfile+'_prof.txt','w')

      if(self.contour.get()==1):
        for i in range(len(self.xprof)):
          line=str("%10.5e" % self.xprof[i])+" "+str("%10.5e" % self.modelprof[i])+" "+str("%10.5e" % self.dataprof[i])
          f.write(line)
      if(self.contour.get()==0):
        for i in range(len(self.xprof)):
          line=str("%10.5e" % self.xprof[i])+" "+str("%10.5e" % self.modelprof[i])+"\n"
          f.write(line)

      f.close()

  ##---------------------------------------------------------------------

    def calculate_radial_prof(self):

      """

      Calculate the radial profile. 

      """

     ## GET THE TYPE OF PROFILE, AND SET THE CENTRE POINT. 

      proftype=self.prof.get()
      profwidth=float(self.profangle.get())
      ax=self.whichax.get()
      cx=ceil(self.npix/2.)
      cy=ceil(self.npix/2.)
      logitx=self.proflogx.get()
      logity=self.proflogy.get()

      ##GET THE LENGTH OF THE PROFILE IN PIXELS

      if(ax==1 or ax==3):
        proflen=int(self.npix/2.)+1

      if(ax==2 or ax==4):
        proflen=int(self.npix/2.)+1

      ##CREATE THE ARRAY
      prof=zeros((proflen))
      npoint=zeros((proflen))
      dataprof=zeros((proflen))
      datapoint=zeros((proflen))

      factor=float(self.datascale.get())

      ##WEDGE PROFILE - AVERAGE OVER A RANGE OF ANGLES. CASES FOR WHICH AXIS. 
      if(proftype==1):
      
        startang=radians(-profwidth/2.+90.*float((ax-1)))
        endang=radians(profwidth/2.+90.*float((ax-1)))
        while(startang < 0):
          startang=startang+radians(360)

        for i in range(self.npix):
          for j in range(self.npix):
            r=int(sqrt((i-cx)**2+(j-cy)**2))
            ang=arctan2(i-cx,j-cy)
            if(r < proflen and ((ax>1 and ang > startang and ang < endang) or (ax==1 and (ang > startang or ang < endang)))):
              prof[r]=prof[r]+self.PI[i,j]
              npoint[r]=npoint[r]+1.
              if(self.contour.get()==1):
                dataprof[r]=dataprof[r]+self.cropcontdata[i,j]
                datapoint[r]=datapoint[r]+1.
 
        for i in range(len(prof)):
          if(npoint[i] > 0):
            prof[i]=prof[i]/npoint[i]
          if(self.contour.get()==1):
            if(datapoint[i] > 0):
              dataprof[i]=dataprof[i]/datapoint[i]
          
        
      ##SLICE PROFILE - AVERAGE OVER A WIDTH, CASES FOR WHICH AXIS. 
      if(proftype==2):
        ##POSITIVE x DIRECTION
        if(ax==1):
          for i in range(proflen):
            for j in range(int(profwidth)):
              ii=i+cx-1
              jj=j+cy-int(profwidth/2.)
              prof[i]=prof[i]+self.PI[ii,jj]
              npoint[i]=npoint[i]+1.
              if(self.contour.get()==1):
                dataprof[i]=dataprof[i]+self.cropcontdata[ii,jj]
                datapoint[i]=datapoint[i]+1.

        if(ax==2):
          for i in range(proflen):
            for j in range(int(profwidth)):
              ii=j+cx-int(profwidth/2.)
              jj=i+cy-1
              prof[i]=prof[i]+self.PI[ii,jj]
              npoint[i]=npoint[i]+1.
              if(self.contour.get()==1):
                dataprof[i]=dataprof[i]+self.cropcontdata[ii,jj]
                datapoint[i]=datapoint[i]+1.

        if(ax==3):
          for i in range(proflen):
            for j in range(int(profwidth)):
              ii=cx-i
              jj=j+cy-int(profwidth/2.)
              prof[i]=prof[i]+self.PI[ii,jj]
              npoint[i]=npoint[i]+1.
              if(self.contour.get()==1):
                dataprof[i]=dataprof[i]+self.cropcontdata[ii,jj]
                datapoint[i]=datapoint[i]+1.


        if(ax==4):
          for i in range(proflen):
            for j in range(int(profwidth)):
              ii=j+cx-int(profwidth/2.)
              jj=cx-i
              prof[i]=prof[i]+self.PI[ii,jj]
              npoint[i]=npoint[i]+1.
              if(self.contour.get()==1):
                dataprof[i]=dataprof[i]+self.cropcontdata[ii,jj]
                datapoint[i]=datapoint[i]+1.




      ##GET COLOURS AND POINT TYPES
 
      self.ptype={'Diamond':"d",'Circle':"o",'Square':"s",'Triangle':"2",'Plus':"+",'X':"x",'None':""}

      modelsym=self.ptype[self.modelplot.get()]
      datasym=self.ptype[self.dataplot.get()]


        ##NOW THE PLOTS
      #clf()
      xx=arange(len(prof))*self.resamp*self.dist*self.pixscale
      xx1=arange(len(prof))


      ##SCALE THE VALUES IF REQUIRED
      dataprof=dataprof*factor

      ##LINEAR CASE, MASK OUT MASKED REGION
      if(logitx==0 and logity==0):
        ind=where(xx1 >= self.radval & (xx1 > 0))
        self.plotit=plot(xx[ind],prof[ind],marker=modelsym,color=self.modelcol.get())
        if(self.contour.get()==1):
         self.plotit=plot(xx[ind],dataprof[ind],marker=datasym,color=self.datacol.get())

      ##XLOG, YLINEAR,  MASK OUT MASKED REGION AND 0 POINT
      if(logitx==1 and logity==0):
        ind=where((xx1 >= self.radval) & (xx1 > 0))
        self.plotit=semilogx(xx[ind],prof[ind],marker=modelsym,color=self.modelcol.get())
        if(self.contour.get()==1):
          self.plotit=semilogx(xx[ind],dataprof[ind],marker=datasym,color=self.datacol.get())

      ##YLOG, XLINEAR, MASK OUT MASKED REGION AND DATA POINTS <=0
      if(logitx==0 and logity==1):
        ind=where((prof > 0) & (xx1 > self.radval) & (xx1 > 0))
        self.plotit=semilogy(xx[ind],prof[ind],marker=modelsym,color=self.modelcol.get())
        if(self.contour.get()==1):
          ind=where((dataprof > 0) & (xx1 > self.radval) & (xx1 > 0))
          self.plotit=semilogy(xx[ind],dataprof[ind],marker=datasym,color=self.datacol.get())

      if(logitx==1 and logity==1):
        ind=where((prof > 0) & (xx1 > self.radval) & (xx1 >0))
        self.plotit=loglog(xx[ind],prof[ind],marker=modelsym,color=self.modelcol.get())
        if(self.contour.get()==1):
          ind=where((dataprof > 0) & (xx1 > self.radval) & (xx1 >0))
          self.plotit=loglog(xx[ind],dataprof[ind],marker=datasym,color=self.datacol.get())

      
      ##SET GLOBAL VARIABLES
     
      self.xprof=xx
      self.modelprof=prof

      if(self.contour.get==1):
        self.dataprof=dataprof

        if(string.lower(sys.platform) != "darwin"):
            show()

 

  ##---------------------------------------------------------------------
       
      
  ##---------------------------------------------------------------------

##RUN THE CODE
a=parameter_view() 
a.pack() 
a.mainloop() 
