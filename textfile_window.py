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


"""

from Tkinter import * 

class textfile_window(Frame): 
  def __init__(self,title,textfile,master=None): 
    Frame.__init__(self,master) 

    """

    Class to show the contents of a text file, with a quit button and scrollbar. 
    Two parameters, the title, and the name of the file. 


    To call: self.helpwin=textfile_window(title,textfile)
             self.helpwin.pack(side=TOP,anchor=W,fill=X)
             self.helpwin.update()

    """

    ##CREATE WINDOW, ADD CLOSE BUTTON

    top=Toplevel(self) 
    top.title(title) 

    fm=Frame(top, relief=RAISED, bd=1) 
    mb2=Button(fm,text='Close',command=self.destroy).pack(side=LEFT,anchor=W)
    fm.pack(side=TOP,anchor=W,fill=X) 

    ##CREATE THE TEXT WINDOW AND ADD THE SCROLLBAR
    fm=Frame(top, relief=RAISED, bd=1)
    self.text=Text(fm,height=30,width=100)
    self.scroll=Scrollbar(fm)
    self.text.focus_set()
    self.scroll.pack(side=RIGHT,fill=Y)
    self.text.pack(side=LEFT,fill=Y)
    self.scroll.config(command=self.text.yview)
    self.text.config(yscrollcommand=self.scroll.set)
    fm.pack(side=TOP,anchor=W,fill=X)

    ##OPEN THE TEXT FILE AND DISPLAY

    try:
      myFile=open(textfile,'r')
      myText= myFile.read() 
      myFile.close() 
    except:
      myText="Cannot Find File "+textfile+"\n"

    self.text.insert(0.0,myText) 

##-----------------------------------------------------------------------
