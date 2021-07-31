#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 19:15:22 2021

@author: mkravche

The main file for generating. To me its entirely obvious what is going on, 
so I am unsure whether to fully annotate this file. 

Order of execution:
	1. Import standard python libraries
	2. Import self written classes
	3. Close existing figures since we are going fullscreen.
	4. Exstantiate the Fields class with the exciter
	
	class Fields(TM_FEM,exciter):
		pass

	# run __init__ from TM_FEM
	F = Fields()
	
	5. Scale the space and time grid to be realistic. Yea I could just use the kwargs I made, but im lazy.
	6. Extract and draw the coordinates defined by the BMP raster.
	7. Setup the figure with the titles. I didn't bother to implmenet blitting since there was no real benefit. 
	Im redrawing all the images anyway. The real performance benefit is initialize once, then reuse and writeover 
	the same image memory. calling imshow each time accrues images and hogs memory resulting in an ever waning computation speed until crash.
	
	8. axis_image. This is a data structure (struct) that holds all the axis and figure pointer handles relevent to
	each of the six displays. It links the function pointer towards the relevant data fields to prevent class clutter.
	This is more to make it more readable.
	self.set_data = self.AxisImage.set_data
        self.set_clim = self.AxisImage.set_clim
        self.autoscale = self.AxisImage.autoscale
        
        9. Alot is happening here. I make a list of `axis_image` which is the datastructure containg all my pointers.
        	it takes the axes handle, makes an AxisImage handle for imshow which it stores.
        	This axes handle is then calling the parent figure handle to display color bar for the specific axis, 
        	while suppling the new AxisImage imshow handle and itself as an argument. 
        	
        	Take note how reference is a function handle. It is called when initialzing imshow for the first time.
        	we then store the function as a handle so we can call it externally at any other time.
        	
        	This is why I created the boiler plate handle functions in FEM_RF. 
        	Numpy doesn't allow references to a memory location as far as I know with standard python.
        	This saves me alot of trivial code implementation.
        	
        	The boundary image is a composite of all the boundaries together. Each field has its own boundary conditions
        	This is them all added together for visualization. In this simple example, there is a little redundancy, but 
        	I kept it this way for more complex cases.
        	
        	
        	
        	figs.AxesImages = []
         	figs.AxesImages.append(axis_image(figs.axes[0][0],F.ex,cmap='seismic'))
		figs.AxesImages.append(axis_image(figs.axes[0][1],F.hzx,cmap='twilight'))
		figs.AxesImages.append(axis_image(figs.axes[0][2],F.dhz_dx,cmap='twilight'))
		figs.AxesImages.append(axis_image(figs.axes[1][0],F.eyT,cmap='seismic'))
		figs.AxesImages.append(axis_image(figs.axes[1][1],F.hzy,cmap='twilight'))
		figs.AxesImages.append(axis_image(figs.axes[1][2],F.boundary_image))
		
	10. Add dithering backround noise that has virtually no influence on the simulation. Breaks up quantization blockyness.
	11. _update(n) function that updates the frames for discrete time n.
		This function deals with setting the color limits and titles to display nicely as well.
	12. Animation functions. This is my special sauce blend of parameters that work for me.
		It requries an indivual install for access to libx264. The apt-get ubuntu repos do not work.
		https://ffmpeg.org/download.html
		
"""
#import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from T_junction import T_junction as Tj
from FEM_RF import sinusoidal_exciter as exciter
from FEM_RF import TM_FEM

from fullfigure import full_figure

# Close existing figures
plt.close('all')

title = "2D FEM Transverse magnetic T-junction simulation"

# Tie excitation function to TM_FEM
# all the functions including __init__ are inherited
# TM_FEM dominates function inheritance
class Fields(TM_FEM,exciter):
    pass

# run __init__ from TM_FEM
F = Fields()

# scale down
F.dx /= 100*2 # 10 for 10cm to 1 cm , 10 for 1 cm to 1 mm
F.dt /=100*2
F.dx *=5
F.dt *=5

# Write T junction on the field
T = Tj(remove_whitespace=False)
F.draw(T)
# Setup animation figure
figs = full_figure(super_title=title,
                   fullscreen=True,
                   subplot_kw=dict( frameon=False,
                                    aspect='equal',
                                    animated=False,
                                    xticks=[],
                                    yticks=[]),)

# Label all titles 
# indexed my [row][col]
figs.axes[0][0].set_title('$e_x$')
figs.axes[0][1].set_title('$h_{zx}$')
figs.axes[0][2].set_title('$\\frac{dh_z}{dx}$')
figs.axes[1][0].set_title('$e_y^T$')
figs.axes[1][1].set_title('$h_{zy}$')
figs.axes[1][2].set_title('Composite \nnpmls boundaries')


# save handle or access ridiculously with figs.axes[0][0].images[0].set_data(...)
# For the most recent member figs.axes[0][0].images[-1].set_data(...)
# keep in mind that repeatedly calling this adds to the images
# lookup members with figs.axes[0][0].__dict__

class axis_image:
    def __init__(self,axis,reference,**imshow_kw):
        """
        Parameters
        ----------
        axis : matplotlib AxisImage
            Handle for imshow
        reference : Function handle which returns working numpy array
            Function handle that returns the sliced, working numpy array
        **imshow_kw : kw arguments for imshow. 
        """
        self.min = 0
        self.max = 0
        self.AxisImage = axis.imshow(reference(),**imshow_kw) # 
        axis.figure.colorbar(self.AxisImage,ax=axis)
        #axis.colorbar()
        self.reference = reference 
        self.set_data = self.AxisImage.set_data
        self.set_clim = self.AxisImage.set_clim
        self.autoscale = self.AxisImage.autoscale
        
    def scale_minmax(self,auto=True):
        a = self.reference().min()
        b = self.reference().max()
        if a < self.min:
            self.min = a
        if b > self.max:
            self.max = b
            
        if auto:
            self.set_clim([self.min,self.max])
        

        
        
# figs.AxesImages.append(figs.axes[0][0].imshow(F.E[:,:,0],cmap='RdBu'))        
figs.AxesImages = []
# figs.AxesImages.append(axis_image(figs.axes[0][0],F.ex,cmap='seismic'))
# figs.AxesImages.append(axis_image(figs.axes[0][1],F.hzx,cmap='seismic'))
# figs.AxesImages.append(axis_image(figs.axes[0][2],F.dhz_dx,cmap='nipy_spectral'))
# figs.AxesImages.append(axis_image(figs.axes[1][0],F.eyT,cmap='seismic'))
# figs.AxesImages.append(axis_image(figs.axes[1][1],F.hzy,cmap='seismic'))
# figs.AxesImages.append(axis_image(figs.axes[1][2],F.boundary_image))
# figs.AxesImages.append(axis_image(figs.axes[0][0],F.ex,cmap='twilight'))
# figs.AxesImages.append(axis_image(figs.axes[0][1],F.hzx,cmap='twilight'))
# figs.AxesImages.append(axis_image(figs.axes[0][2],F.dhz_dx,cmap='nipy_spectral'))
# figs.AxesImages.append(axis_image(figs.axes[1][0],F.eyT,cmap='twilight'))
# figs.AxesImages.append(axis_image(figs.axes[1][1],F.hzy,cmap='twilight'))
# figs.AxesImages.append(axis_image(figs.axes[1][2],F.boundary_image))
figs.AxesImages.append(axis_image(figs.axes[0][0],F.ex,cmap='seismic'))
figs.AxesImages.append(axis_image(figs.axes[0][1],F.hzx,cmap='twilight'))
figs.AxesImages.append(axis_image(figs.axes[0][2],F.dhz_dx,cmap='twilight'))
figs.AxesImages.append(axis_image(figs.axes[1][0],F.eyT,cmap='seismic'))
figs.AxesImages.append(axis_image(figs.axes[1][1],F.hzy,cmap='twilight'))
figs.AxesImages.append(axis_image(figs.axes[1][2],F.boundary_image))

# Dither
eps = 1e-6
dither = np.random.random(F.E[:,:,1].shape) * eps
F.E[:,:,1] *= dither

# As before in previous projects, we do not re-call imshow again and access the 
# image to redraw it.
#  access with figs.AxesImages[0].set_data(mm)
#                               a.get_array().data
#              figs.AxesImages[0].set_clim([mm.min(),mm.max()])
#                               a.get_clim()
# The axis_image class links too the functions directly too cut down on syntax clutter

def _update(n):
    bb = title+"\n"+ "frame: " + "{:10}, t={:10.3f}"+'ns'+', frequency= '+ str(F.freq/1e9)+ 'Ghz'+ ', feedline_width = '+ str(round(100*9*F.dx,1)) + 'cm'
    #print(bb)
    bb = bb.format(n,round(1e9*n*F.dt,3))
    figs.fig.suptitle(bb)
    F.update_fields(n)
    
    g = 0
    for axe in figs.AxesImages:
        
        axe.set_data(axe.reference())
        #axe.set_clim(0,1)
        if g == 0 or g == 3:
            #print(g, "minmax scale")
            #axe.scale_minmax()
            axe.set_clim([-50,50])
        elif g == 1 or g == 4:
            axe.set_clim([-1.1,1.1])
        elif g == 2:
            axe.set_clim([-0.01,0.01])
        else:
            #print(g, "autoscale")
            axe.autoscale()
        g +=1
        
        
   

    
nmax = 2100*2

myanimation = animation.FuncAnimation(figs.fig, _update,frames=nmax)

# The feed width is 9*dx pixels wide in the example made
# def get_addr(addr):
#     return [x for x in globals().values() if id(x) == addr]

plt.rcParams['animation.ffmpeg_path'] = r"/opt/ffmpeg/ffmpeg-4.2.1-linux-64/ffmpeg"
def save_animation(filename,fps=None,dpi=None,bitrate=28000):
    Writer = animation.writers['ffmpeg']
    #0 bitrate is lossless for libx264
    writer = Writer(fps=fps,metadata=dict(artist='Mikhail Kravchenko'),
                    bitrate=bitrate,
                    extra_args=['-vcodec', 'libx264'])
    myanimation.save(filename,writer=writer,dpi=dpi)

import datetime
d = datetime.datetime.now()
f = 'T-junction-'+ '{date:%Y-%m-%d_%H_%M_%S}' +'.mkv'
f = f.format(date=d)
save_animation(f,fps=30,dpi=200)

