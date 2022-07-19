#!/usr/bin/python
import tkinter as Tk
from tkinter import filedialog

import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch, simple_norm
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
#NavigationToolbar2Tk) 

from os.path import basename
from astropy.io import fits
import numpy as np
import subprocess
from time import sleep

# Donuts part (toko)
import don11
import common
filename_params = 'donut.par'
DS9_NAME = 'KGO'
import pyds9 # pip install pyds9
from pyds9 import DS9 as ds9 
import parse as prs # pip install parse

class Ds9_stuff:
  def __init__(self):
     print('It\'s DS9 constructor. Nothing more')
     self.ds9 = -1

  def find_ds9(self):
     try:
       ds9_list = pyds9.ds9_targets()
       if(len(ds9_list) > 1): self.ds9 = ds9('DS9:'+DS9_NAME)
       else: self.ds9 = ds9()
#       print('ds9=',self.ds9)
     except Exception as e:
       print('Ds9_stuff:find_ds9 exception ',e)
       self.ds9 = -1
       return -1
     return 0

  def view(self, fullpath):
    if self.find_ds9() < 0:
      subprocess.call('./ds9_start.sh')
      sleep(1)
      if self.find_ds9() < 0: return
      else: self.ds9.set('file ' + fullpath)
    else:
      try:
        self.ds9.set('file ' + fullpath)
      except Exception as e:
        print('ds9.set file '+fullpath + 'error:', e)


  def find_all_circle_regions(self, regions_string):
    reg = regions_string
    circ_list = []
    pos = reg.find('circle')
    while pos >= 0:
      circ = prs.search('circle({:g},{:g},{:g}', reg, pos=pos)
      circ_list.append(circ)
      pos = reg.find('circle',pos+1)
    return circ_list

  def get_circle_region(self, which=0):   # number in the list, m.b 0,1, ... -1,-2
     if self.find_ds9() < 0:
       print('No ds9 found')
       return 0,0,0
     try:
        regs = self.ds9.get('region').replace('\"','')
     except Exception as e:
        print('Ds9_stuff:find_ds9 exception ',e)
        return -1
     circ_list = self.find_all_circle_regions(regs)
     if which in range(-len(circ_list), len(circ_list)):
       circ = circ_list[which]
       return circ[0], circ[1], circ[2] # ordinary image x_centre, y_centre, radius
     return 0,0,0


class Display():
  def __init__(self, parent):
        self.fig = Figure(figsize=(8,4), dpi=100, facecolor='black')
        px = 0.01
        self.fig.subplots_adjust(bottom=px, top=1-px, left=px, right=1.0-px, wspace=px)
        self.ax_left  = self.fig.add_subplot(121)
        self.ax_left.axis('off')
        self.ax_right = self.fig.add_subplot(122)
        self.ax_right.axis('off')
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)  # creating the Tkinter canvas containing the Matplotlib...
        self.canvas.get_tk_widget().pack()


  def plot(self, img, side='right'): 

#        norm = simple_norm(img, 'sqrt', percent=99.99)
        if side.upper() != 'LEFT':
#          self.ax_right.imshow(img, norm=norm, origin='lower', cmap='Greys_r'), interpolation='nearest')
          self.ax_right.imshow(img, origin='lower', cmap='Greys_r', interpolation='nearest')
        else:
          self.ax_left.imshow(img, origin='lower', cmap='Greys_r', interpolation='nearest')

        self.canvas.draw() 

# ------------------ Main Window -------------------

class donutsGUI():

  def select_fits(self):
    fits_full_name = Tk.filedialog.askopenfilename(initialdir = ".",\
                    title = "Select file",\
                    filetypes = (("fits files","*.fit*"),\
                    ("all files","*.*")))
    print(fits_full_name)
    if len(fits_full_name) == 0:
       return

    try:
      hdulist = fits.open(fits_full_name, uint=False)  # working with float data
      nhdu = len(hdulist)
      self.img_tot  = hdulist[nhdu-1].data
      self.extract_button.configure(state='normal')
      self.get_center_button.configure(state='normal')
      self.extract()
#     head = hdulist[0].header
    except Exception as e:
      print('Error:', e)
      self.fits_tkvar.set('error')
      return

    name = basename(fits_full_name).split('.')[0]
    self.fits_tkvar.set(name)
    self.ds9_stuff.view(fits_full_name)

  def extract(self):
    if self.img_tot is None: return
    self.impix = don11.extract(self.img_tot, self.xc_tkvar.get(), self.yc_tkvar.get(), common.fovpix)
    self.display.plot(self.impix, side='left')
    self.fit_button.configure(state='normal')

  def init_par(self):
    don11.init(filename_params)
    self.xc_tkvar.set(common.donpar.xc)
    self.yc_tkvar.set(common.donpar.yc)

    inv_dict = {v: k for k, v in self.focus_pos_choices.items()}

    focus_pos_name = inv_dict.get(common.donpar.efoc, None)  # intra +1  extra -1
    self.focus_pos.set(focus_pos_name)

  def fit(self):
    # correct params:
    common.donpar.seeing = self.seeing0_tkvar.get()
    common.donpar.efoc   = self.focus_pos_choices.get(self.focus_pos.get(),1)

    print('Fit it!')
    if self.impix is None:
      print('impix is empty. Why????')
      return
    if self.display_tkvar.get():
      display = self.display.plot
    else: display = None
    zres, impix, chi2 = don11.fit(self.impix, display=display)
    num = 0
    for idx in self.zernike_selected_indices:
      self.zernike_selected_list[num].set(round(zres[idx],3))
      num += 1
    self.rms_tkvar.set(round(chi2*100,2))
    self.focus_mkm_tkvar.set(round(zres[3]*(-2)*8**2,  3))
    print('chi2 =', chi2, 'zres =', zres)
    self.display.plot(impix, side='right')

  def extract_from_ds9(self):
      x,y,rad = self.ds9_stuff.get_circle_region()
      if rad <= 0: return
      self.xc_tkvar.set(x)
      self.yc_tkvar.set(y)
      self.extract()

  def __init__(self, parent):
    self.parent = parent
# Input:
    self.seeing0_tkvar = Tk.DoubleVar(value=0.5)
    self.xc_tkvar      = Tk.IntVar(value=300)
    self.yc_tkvar      = Tk.IntVar(value=300)
    self.fits_tkvar    = Tk.StringVar(value='select fits:')
    self.display_tkvar = Tk.BooleanVar(value='False')

# Output:
    self.seeing_tkvar    = Tk.DoubleVar(value=0.0)
    self.rms_tkvar       = Tk.DoubleVar(value=0.0)
    self.focus_mkm_tkvar = Tk.DoubleVar(value=0.0)

    pdx = 5; pdy = 5
    frame_upper = Tk.LabelFrame(parent, text='Input', borderwidth=3, relief=Tk.RIDGE, padx=pdx, pady=pdy)
    frame_upper.pack(side='top', fill='both')
    self.img_tot = None
    self.impix = None

    lwidth = 7
    seeing0_label = Tk.Label(frame_upper, text = '~seeing,\"')
    seeing0_entry = Tk.Entry(frame_upper, textvariable=self.seeing0_tkvar, width=5)
    xc_label = Tk.Label(frame_upper, width=lwidth, text = 'xc')
    xc_entry = Tk.Entry(frame_upper, textvariable=self.xc_tkvar, width=5)
    yc_label = Tk.Label(frame_upper, width=lwidth, text = 'yc')
    yc_entry = Tk.Entry(frame_upper, textvariable=self.yc_tkvar, width=5)

    self.focus_pos = Tk.StringVar()
    self.focus_pos_choices = {'INTRA':1, 'EXTRA':-1}
    focus_pos_label        = Tk.Label(frame_upper,  width=11, text='focus pos')
    focus_menu             = Tk.OptionMenu(frame_upper, self.focus_pos, *self.focus_pos_choices) # foreground='blue'
    display_intermediate   = Tk.Checkbutton(frame_upper, var=self.display_tkvar, text='display')

    fits_label            = Tk.Label (frame_upper, width=14,     textvariable=self.fits_tkvar, foreground='blue', relief=Tk.SUNKEN)
    select_fits_button    = Tk.Button(frame_upper, text='Select fits', command=self.select_fits)
    self.extract_button   = Tk.Button(frame_upper, width=lwidth, text='From xc,yc', state='disabled', command=self.extract)
    self.get_center_button= Tk.Button(frame_upper, width=lwidth, text='From ds9',   state='disabled', command=self.extract_from_ds9)
    self.fit_button       = Tk.Button(frame_upper, width=lwidth, text='Fit',        state='disabled', command=self.fit)


    self.ds9_stuff = Ds9_stuff()

    row = 0; column = 0;
    for w in [seeing0_label, xc_label, yc_label, focus_pos_label, fits_label, display_intermediate]:
       w.grid(row=row, column=column, padx=1, sticky=Tk.EW)
       column += 1
    row += 1; column = 0;
    for w in [seeing0_entry, xc_entry, yc_entry, focus_menu, select_fits_button, self.get_center_button,\
                                  self.extract_button, self.fit_button]:
       w.grid(row=row, column=column, padx=1, sticky=Tk.EW)
       column += 1

    frame_middle = Tk.Frame(parent, borderwidth=3, relief = Tk.RIDGE, padx=pdx, pady=pdy)
    frame_middle.pack(side='top')
    self.display = Display(frame_middle)

    frame_bottom = Tk.LabelFrame(parent, borderwidth=3, text='Output', relief=Tk.RIDGE, padx=pdx, pady=pdy)
    frame_bottom.pack(side='top', fill='both')

    row = 0; column = 0;
    
    self.zernike_selected_list = []
    zernike_selected_names       =  ['seeing\"', 'defocus', '/ astigm', '| astigm', '| coma', '-- coma', '| tref','/ tref', 'spher']
    self.zernike_selected_indices = [0,          3,             4,               5,        6,         7,     8,     9,  10 ]
    if len(zernike_selected_names) != len(self.zernike_selected_indices):
      sys.exit('Error: zernike_selected_names do not fit the zernike_selected_indices')

    lwidth = 7
    
    Tk.Label(frame_bottom, text='RMS,%', width=lwidth, padx=pdx, pady=pdy).grid(row=row, column=column, padx=3)
    column += 1
    for name in zernike_selected_names:
       label = Tk.Label(frame_bottom, text=name, width=lwidth, padx=pdx, pady=pdy)
       label.grid(row=row, column=column, padx=3)
       column += 1

    row += 1; column = 0
    Tk.Label(frame_bottom, width=lwidth, textvariable=self.rms_tkvar, foreground='blue', relief=Tk.SUNKEN, padx=pdx, pady=pdy).grid(row=row, column=column, padx=2)
    column += 1
    for idx in self.zernike_selected_indices:
      tkvar = Tk.DoubleVar(value=0.0)
      self.zernike_selected_list.append(tkvar)
      label = Tk.Label(frame_bottom, width=lwidth, textvariable=tkvar, foreground='blue', relief=Tk.SUNKEN, padx=pdx, pady=pdy)
      label.grid(row=row, column=column, padx=2)
      column += 1

    sai25focus_label = Tk.Label(frame_bottom,  text='SAI25 defocus, mkm')
    sai25focus_mkm   = Tk.Label(frame_bottom, width=lwidth, textvariable=self.focus_mkm_tkvar, foreground='blue', relief=Tk.SUNKEN)
    row += 1; column = 0;
    sai25focus_label.grid(row=row, columnspan=2, column=column, padx=1, sticky=Tk.EW)
    column += 2
    sai25focus_mkm.grid(row=row, column=column, padx=1, sticky=Tk.EW)
      


# Do something
    self.init_par()


root = Tk.Tk() 


root.title('Donuts') 
gui = donutsGUI(root)
root.mainloop() 
