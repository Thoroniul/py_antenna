#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 17:52:38 2021

@author: mkravche
"""

import numpy as np
from PIL import Image

class T_junction:
"""
T_junction(image_path="T_junction.bmp",remove_whitespace=True)
INPUT
image.bmp -> This is the input image which I kept the same as the raster size of the simulation.
		It doesn't have to be the same size. The coordinates are extracted from the image.
		In this version there are only ideal walls and source blocks.
		In this version: 
		The black pixels are the wall and are detected by finding pixels with B=0
		The 'pink' pixels are the source coordinates and are for pixels with G=127
		BMP has BGR pixel orientation.
		
remove_whitespace -> This subtracts the minimum value from the coordinates as an offset. 
			This is an appendage for when I command the coordinates with a local handler object.
			
BUG Note there is a slight bug in which it only takes the minimum from the wall coordinates and applys it to 
both the source and walls. If by chance a source pixel precedes a wall, then we may end up offsetting the source 
off the raster. This is corrected in the C++ version.
"""
    def __init__(self,image_path="T_junction.bmp",remove_whitespace=True):
        self.image = Image.open(image_path)
        self.image.load()
        
        #Takes black pixels
        self.mask = np.array(self.image)[:,:,0] == 0
        
        #Return indexes of black pixel location
        self.coords = list(np.where(self.mask))
        
        # Takes pink pixels as source
        self.source = list(np.where(np.array(self.image)[:,:,1] == 127)) #Fast way for pink
        
        #Subtract whitespace by subtracting offset
        # This assumes the source in pink is preceded by black pixels.
        if remove_whitespace:
            self.offset = np.min(self.coords,axis=1)
            self.coords[0] -= self.offset[0]
            self.coords[1] -= self.offset[1]
      
            self.source[0] -= self.offset[0]
            self.source[1] -= self.offest[1]
            
    def T(self):
        # Functions that returns a tuple wrapping the transposed coordinates
        return (self.coords[1],self.coords[0])
            
