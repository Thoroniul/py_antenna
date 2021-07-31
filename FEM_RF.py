#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 17:52:22 2021

@author: mkravche

	This file contains the classes that define the EM fields. 
	This is me consolidating (or attempting to consoldiate) and documenting various code and notes I have lazily accreted in the last few years.
	I designed it so that the exciter is linked separatly. 

	Contains:

	class sinusoidal_exciter(dt,frequency):       ------>       s.exciter(n) -> (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)
	A signal sinusoidal signal generator that returns a floating point number for a point n.
	The (1.0 - exp( (-n/k/2)^2 )*Z  is a ramp up the sinusoid where Z is the impedance of the wave in free space. aimp is air impedance -_-.
	When dealing with multi-material boundary transisions account for the difference in wave impedance to properly observe snell's law.


	It is consolidated with the TM_FEM class in the main executing function using the boilerplate:
	When doing so, it takes the __init__function of TM_FEM which has precidence. This works out nicely so we do not define the constants twice.
			class Fields(TM_FEM,exciter):
			    pass

	It is applied the field with the defined helper function as part of the field manager crudly named 	    
	INPUT:
	dt -> difference time step (s)
	frequency -> frequency (hz)
	n -> An integer corresponding to the frame number. Realistically it will take floats as well, nothing stopping it. 

	OUTPUT:
	-> SIGNAL (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)


		
	class TM_FEM(npmls=8,m=100,n=100,  frequency=2.4e9,dx=None,dt=None): 
		
	INPUT:
	npmls is the layers of the perfect absorbing material
	m,n = stage raster size
	frequency = frequency of RF
	dx = cell size, defaults to lambda/2/6.0 (m) (this is an arbitrary gridsize so that dx and dt are constrained to the speed of light) 
	dt = time cell size (s)

	OUTPUT:
	Any field. It really is meant to be a data storage class with functions contained to operate on it.


	VECTORS:
	s.E[:,:,0] --> Ex
	s.E[:,:,1] --> Ey.T

	s.cE[:,:,0] --> caex
	s.cE[:,:,1] --> cbex
	s.cE[:,:,2] --> caey.T
	s.cE[:,:,3] --> cbey.T

	s.hz[:,:,0] --> hz
	s.hz[:,:,1] --> hzx
	s.hz[:,:,2] --> hzy
	s.hz[:,:,3] --> dahzx
	s.hz[:,:,4] --> dbhzx
	s.hz[:,:,5] --> dahzy
	s.hz[:,:,6] --> dbhzy

	ACCESS FUNCTIONS:
	# function handles for pseudo pass by reference. No reason one can't access the vectors directly in python. 
	To link to the animator it was much easier to provide a function as a parameter than to give a numpy array.
	def ex(self):
		return self.E[:,:,0]

	def eyT(self):
		return self.E[:,:,1].T

	def hzx(self):
		return self.hz[:,:,1]

	def dhz_dx(self):
		i,j = self.theater
		return self.hz[:i,1:j,0] - self.hz[:i,:j-1,0]

	def hzy(self):
		return self.hz[:,:,2]

	def boundary_image(self):
		return self.cE[:,:100,3].T + self.cE[:,:100,1]
		
	OPERATION FUNCTIONS:
	These functions operate on the field vectors.

	def update_fields(self,n=0,excitation=True):
	update_fields takes a positive time n and induces the changes in the field.

	D=epsion*E a change in the electromagnetic field induces a displacement field in the dielectric medium.
	The magnetic moment of the dielectric induces a change in the magnetic field H (A/m). 
	This magnetic field induces a current with Ampere's law and influences and induction B (N/A-m).
	H = B/muo

	E (V/m) (N/C)
	D (C/m^2)
	B (N/Am)
	H (A/m)

	To avoid dealing with complex numbers, I kept Ex and Ey separate in line with my reference textbook.

	The exciter is linked here by the self.draw function. Which "draws" the extracted coordinates from the BMP.

	if excitation: 
	    s = self.source
	    self.E[s[1],s[0],1] = self.exciter(n)

	self.source must be exstantiated. This is done by the self.draw function.
	Since the field is transpoed, the coordinate indicies are swapped.

		
	def s._define_PML_constants():
	Hidden helper function to exstantiate the Berenger PML constants for a "perfectly transparent" boundary around the edges of the FEM.

	http://www.sam.math.ethz.ch/~hiptmair/Seminars/ABC/slides/Krish.pdf

	def make_boundaries(self,A,border_vec,orientation='lr',padding='',inline = True):

	INPUT
	A -> Field Vector operated on
	border_bec -> material PML vector 
	orientation='lr' or 'tb' -> boundary to apply along vertical boundaries (top and bottom) or lateral boundaresi( left and right)
	padding -> adds a gap depending on whether a characther exists in the padding string. t-> top, b-> bottom, l-> left, r-> right.
	inline-> true controls the output type. If inline it operates on the numpy vector directly
	OUTPUT
	A -> Outputs the same A if the inline is False. Otherwised it modifies it directly.

	def draw(T):
	'draws' the coordinates extracted from the bitmap onto the field equations.
	It draws a value of -1 for a perfect reflector for the ca field (ca (ex, ey) ) 
	During the difference equation operating the FEM, this will flip the sign and send the wave in the opposite direction.
	It draws 0 for (cb (ex, ey) ). Regardless of any change in the field. These will have a charge 0 which set our boundary
	condition to the standing wave in the cavity between walls.

	These draw idealistic boundaries in a raster. Further improvements can be added by adding lossy materials, anisotropic materials, 
	and materials of different dielectric and wave propogation constants. For now this is just a filler.

	INPUT 
	T -> class containing the extracted coordinates of the BMP image.

	

"""

import numpy as np

class sinusoidal_exciter:
    def __init__(self,dt,frequency):
    """
    class sinusoidal_exciter(dt,frequency):       ------>       s.exciter(n) -> (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)
	A signal sinusoidal signal generator that returns a floating point number for a point n.
	The (1.0 - exp( (-n/k/2)^2 )*Z  is a ramp up the sinusoid where Z is the impedance of the wave in free space. aimp is air impedance -_-.
	When dealing with multi-material boundary transisions account for the difference in wave impedance to properly observe snell's law.


	It is consolidated with the TM_FEM class in the main executing function using the boilerplate:
	When doing so, it takes the __init__function of TM_FEM which has precidence. This works out nicely so we do not define the constants twice.
			class Fields(TM_FEM,exciter):
			    pass

	It is applied the field with the defined helper function as part of the field manager crudly named 	    
	INPUT:
	dt -> difference time step (s)
	frequency -> frequency (hz)
	n -> An integer corresponding to the frame number. Realistically it will take floats as well, nothing stopping it. 

	OUTPUT:
	-> SIGNAL (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)	"""
        self.pi = np.pi
        self.freq=frequency
        self.epso = 8.854e-12                      # Permittivity of free space
        self.muo = 4.0e-7*self.pi             
        self.aimp = np.sqrt(self.muo/self.epso)    # Wave Impedance in free space | air impedance
        #self.butt = 'butt'
    def exciter(self,n):
    	"""
    	s.exciter(n) -> (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)
    	INPUT:
    	n -> An integer corresponding to the frame number. Realistically it will take floats as well, nothing stopping it. 
    	OUTPUT:
    	-> SIGNAL (1.0 - exp( (-n/k)^2 ) * Z * sin(2pi*f*dt*n)
    	"""
        return (1.0 - np.exp(-((n/10.0/2)**2)))*self.aimp*np.sin(2.0*self.pi*self.freq*self.dt*n)
    


class TM_FEM:
    def __init__(self, npmls=8,m=100,n=100,  frequency=2.4e9,dx=None,dt=None):
        """
        class TM_FEM(npmls=8,m=100,n=100,  frequency=2.4e9,dx=None,dt=None): 
        	
        INPUT:
        npmls is the layers of the perfect absorbing material
        m,n = stage raster size
        frequency = frequency of RF
        dx = cell size, defaults to lambda/2/6.0 (m) (this is an arbitrary gridsize so that dx and dt are constrained to the speed of light) 
        dt = time cell size (s)
        
        OUTPUT:
        Any field. It really is meant to be a data storage class with functions contained to operate on it.


	VECTORS:
        s.E[:,:,0] --> Ex
        s.E[:,:,1] --> Ey.T
        
        s.cE[:,:,0] --> caex
        s.cE[:,:,1] --> cbex
        s.cE[:,:,2] --> caey.T
        s.cE[:,:,3] --> cbey.T

        s.hz[:,:,0] --> hz
        s.hz[:,:,1] --> hzx
        s.hz[:,:,2] --> hzy
        s.hz[:,:,3] --> dahzx
        s.hz[:,:,4] --> dbhzx
        s.hz[:,:,5] --> dahzy
        s.hz[:,:,6] --> dbhzy
        
        ACCESS FUNCTIONS:
	# function handles for pseudo pass by reference. No reason one can't access the vectors directly in python. 
	To link to the animator it was much easier to provide a function as a parameter than to give a numpy array.
	def ex(self):
		return self.E[:,:,0]

	def eyT(self):
		return self.E[:,:,1].T

	def hzx(self):
		return self.hz[:,:,1]

	def dhz_dx(self):
		i,j = self.theater
		return self.hz[:i,1:j,0] - self.hz[:i,:j-1,0]

	def hzy(self):
		return self.hz[:,:,2]

	def boundary_image(self):
		return self.cE[:,:100,3].T + self.cE[:,:100,1]
		
	OPERATION FUNCTIONS:
	These functions operate on the field vectors.
	
	def update_fields(self,n=0,excitation=True):
	update_fields takes a positive time n and induces the changes in the field.
	
	D=epsion*E a change in the electromagnetic field induces a displacement field in the dielectric medium.
	The magnetic moment of the dielectric induces a change in the magnetic field H (A/m). 
	This magnetic field induces a current with Ampere's law and influences and induction B (N/A-m).
	H = B/muo
	
	E (V/m) (N/C)
	D (C/m^2)
	B (N/Am)
	H (A/m)
	
	To avoid dealing with complex numbers, I kept Ex and Ey separate in line with my reference textbook.
	
	The exciter is linked here by the self.draw function. Which "draws" the extracted coordinates from the BMP.
	
	if excitation: 
            s = self.source
            self.E[s[1],s[0],1] = self.exciter(n)
        
        self.source must be exstantiated. This is done by the self.draw function.
        Since the field is transpoed, the coordinate indicies are swapped.
        
                
        def s._define_PML_constants():
        Hidden helper function to exstantiate the Berenger PML constants for a "perfectly transparent" boundary around the edges of the FEM.
        
        http://www.sam.math.ethz.ch/~hiptmair/Seminars/ABC/slides/Krish.pdf
        
	def make_boundaries(self,A,border_vec,orientation='lr',padding='',inline = True):
    	
    	INPUT
    	A -> Field Vector operated on
    	border_bec -> material PML vector 
    	orientation='lr' or 'tb' -> boundary to apply along vertical boundaries (top and bottom) or lateral boundaresi( left and right)
    	padding -> adds a gap depending on whether a characther exists in the padding string. t-> top, b-> bottom, l-> left, r-> right.
    	inline-> true controls the output type. If inline it operates on the numpy vector directly
    	OUTPUT
    	A -> Outputs the same A if the inline is False. Otherwised it modifies it directly.
    	
    	def draw(T):
    	'draws' the coordinates extracted from the bitmap onto the field equations.
    	It draws a value of -1 for a perfect reflector for the ca field (ca (ex, ey) ) 
    	During the difference equation operating the FEM, this will flip the sign and send the wave in the opposite direction.
    	It draws 0 for (cb (ex, ey) ). Regardless of any change in the field. These will have a charge 0 which set our boundary
    	condition to the standing wave in the cavity between walls.
    	
    	These draw idealistic boundaries in a raster. Further improvements can be added by adding lossy materials, anisotropic materials, 
    	and materials of different dielectric and wave propogation constants. For now this is just a filler.
    	
    	INPUT 
    	T -> class containing the extracted coordinates of the BMP image.
        """
        #self.butt = 'butt'
        self.theater = (m,n)
        self.npmls = 8
        
        self.pi = np.pi
        self.muo = 4.0*self.pi*1.0e-7
        self.epso = 8.854e-12        # Permittivity of ree space
        self.co = 1.0/np.sqrt(self.muo*self.epso) # Speed of light in free space
        self.aimp = np.sqrt(self.muo/self.epso)    # Wave Impedance in free space | air impedance
        self.freq = frequency
        self.lamb = self.co/frequency
        
        self.dx = self.lamb/2/6.0
        self.dt = self.dx/self.co/2.0
        if dx:
            self.dx = dx        # Finite Difference Time Domain (FDTD) cell size
        if dt:
            self.dt = dt/self.co/2.0  # Time Step for speed of light
        
        #Startup
        self._define_PML_constants();
        self._init_fields()
        self._init_boundaries()
        
    def update_fields(self,n=0,excitation=True):
    """
    	def update_fields(self,n=0,excitation=True):
	update_fields takes a positive time n and induces the changes in the field.
	
	D=epsion*E a change in the electromagnetic field induces a displacement field in the dielectric medium.
	The magnetic moment of the dielectric induces a change in the magnetic field H (A/m). 
	This magnetic field induces a current with Ampere's law and influences and induction B (N/A-m).
	H = B/muo
	
	E (V/m) (N/C)
	D (C/m^2)
	B (N/Am)
	H (A/m)
	
	The exciter is linked here.
	if excitation: 
            s = self.source
            self.E[s[1],s[0],1] = self.exciter(n)
        
        self.source must be exstantiated. This is done by the self.draw function.
        Since the field is transpoed, the coordinate indicies are swapped.
    """
        # Update Ex
        #100x99
        i,j = self.theater
        delta_hz = self.hz[:i,1:j,0] - self.hz[:i,:j-1,0]
        
        self.E[:i,1:j,0] =  self.cE[:i,1:j,0]*self.E[:i,1:j,0] + \
                            self.cE[:i,1:j,1]*delta_hz
        # EY update
        # 99x100
        delta_hz_c = self.hz[:i-1,:j,0] - self.hz[1:i,:j,0]
        Ey =    self.cE[:,:,2].T[1:i,:j]*self.E[:,:,1].T[1:i,:j] + \
                self.cE[:,:,3].T[1:i,:j]*delta_hz_c
        self.E[:j,1:i,1] = Ey.T
        
                            
        if excitation: 
            s = self.source
            self.E[s[1],s[0],1] = self.exciter(n)
            
        # HZ update 
        # 100 x 100
        delta_ey_T = self.E[:,:,1].T[:i,:j] - self.E[:,:,1].T[1:i+1,:j]
        self.hz[:i,:j,1] =  self.hz[:i,:j,3] * self.hz[:i,:j,1] + \
                            self.hz[:i,:j,4] * delta_ey_T
        
        delta_ex = self.E[:i,1:j+1,0] - self.E[:i,:j,0]
        self.hz[:i,:j,2] =  self.hz[:i,:j,5] * self.hz[:i,:j,2] + \
                            self.hz[:i,:j,6] * delta_ex
        
        # hz = hzx + hzy
        self.hz[:,:,0] = self.hz[:,:,1] + self.hz[:,:,2]
            
            
    def _define_PML_constants(self):
        """
        def s._define_PML_constants():
        Hidden helper function to exstantiate the Berenger PML constants for a "perfectly transparent" boundary around the edges of the FEM.
        
        https://en.wikipedia.org/wiki/Perfectly_matched_layer
        http://www.sam.math.ethz.ch/~hiptmair/Seminars/ABC/slides/Krish.pdf
        """
        # Berenger PML material constants
        sigmax = -3.0*self.epso*self.co*np.log(1e-15)/(2.0*self.dx*self.npmls)
        rhomax = sigmax*(self.aimp**2)
        
        m = np.arange(1,self.npmls+1)
        self.sig = sigmax*((m-0.5)/(self.npmls+0.5))**2
        self.rho = rhomax*(m/(self.npmls+0.5))**2
        
        class PML:
            def __init__(self):
                pass
        self.PML = PML()
            
        self.PML.re = self.sig*self.dt/self.epso
        self.PML.rm = self.rho*self.dt/self.muo
        self.PML.ca = np.exp(-self.PML.re)
        self.PML.cb = -1*(np.exp(-self.PML.re)-1.0)/self.sig/self.dx
        self.PML.da = np.exp(-self.PML.rm)
        self.PML.db = -(np.exp(-self.PML.rm)-1.0)/self.rho/self.dx
        
    def _init_fields(self):
        """
        initializes the 13 field vectors.
        """
        m,n = self.theater
        
        """
        E[:,:,0] --> Ex
        E[:,:,1] --> Ey.T
        """
        self.E = np.zeros( (m,n+1,2) )
        
        """
        cE[:,:,0] --> caex
        cE[:,:,1] --> cbex
        cE[:,:,2] --> caey.T
        cE[:,:,3] --> cbey.T
        """
        cE = np.ones( (m,n+1,4) )  # cae
        cE[:,:,[1,3]] = self.dt/self.epso/self.dx # 
        self.cE = cE
        
        """
        hz[:,:,0] --> hz
        hz[:,:,1] --> hzx
        hz[:,:,2] --> hzy
        hz[:,:,3] --> dahzx
        hz[:,:,4] --> dbhzx
        hz[:,:,5] --> dahzy
        hz[:,:,6] --> dbhzy
        """
        hz = np.zeros( (m,n,7) )
        hz[:,:,[3,5]] = 1 # dahz
        hz[:,:,[4,6]] = self.dt/self.muo/self.dx #dbhz
        self.hz = hz
    
    def _init_boundaries(self):
    	"""
    	Calls makeboundaries to set up the boundaries for all the field vectors.
    	
    	Note that caey and cbey (cE[:,:,[1,3]]) are transposed.
    	"""
        # caex = make_boundaries(caex, ca,padding='lrt')
        # cbex = make_boundaries(cbex, cb,padding='lrt')
        # caey = make_boundaries(caey, ca, orientation='tb',padding='ltb')
        # cbey = make_boundaries(cbey, cb, orientation='tb',padding='ltb')
        # dahzx = make_boundaries(dahzx,da,orientation='tb',padding='l')
        # dbhzx = make_boundaries(dbhzx,db,orientation='tb',padding='l')
        # dahzy = make_boundaries(dahzy,da,padding='t')
        # dbhzy = make_boundaries(dbhzy,db,padding='t')
        
        self.make_boundaries(self.cE[:,:,0], self.PML.ca, padding='lrt' )
        self.make_boundaries(self.cE[:,:,2], self.PML.ca, padding='lrt' )
        self.make_boundaries(self.cE[:,:,1], self.PML.cb, padding='lrt' )
        self.make_boundaries(self.cE[:,:,3], self.PML.cb, padding='lrt' )
        
        self.make_boundaries(self.hz[:,:,3],self.PML.da,orientation='tb',padding='l')
        self.make_boundaries(self.hz[:,:,4],self.PML.db,orientation='tb',padding='l')
        self.make_boundaries(self.hz[:,:,5],self.PML.da,padding='t')
        self.make_boundaries(self.hz[:,:,6],self.PML.db,padding='t')
    

        
    def make_boundaries(self,A,border_vec,orientation='lr',padding='',inline = True):
    	"""
    	def make_boundaries(self,A,border_vec,orientation='lr',padding='',inline = True):
    	
    	INPUT
    	A -> Field Vector operated on
    	border_bec -> material PML vector 
    	orientation='lr' or 'tb' -> boundary to apply along vertical boundaries (top and bottom) or lateral boundaresi( left and right)
    	padding -> adds a gap depending on whether a characther exists in the padding string. t-> top, b-> bottom, l-> left, r-> right.
    	inline-> true controls the output type. If inline it operates on the numpy vector directly
    	OUTPUT
    	A -> Outputs the same A if the inline is False. Otherwised it modifies it directly.
    	"""
        m,n = A.shape
        npmls = len(border_vec)
        
        io = 0
        ii = m
        
        jo = 0
        ji = n
        
        if 't' in padding:
            io+=1
        if 'b' in padding:
            ii-=1
        if 'l' in padding:
            jo+=1
        if 'r' in padding:
            ji-=1
        
        if orientation == 'lr':
            A[io:ii,jo:jo+npmls] = np.repeat([np.flip(border_vec)],ii-io,axis=0)
            A[io:ii,ji-npmls:ji] = np.repeat([border_vec],ii-io,axis=0)
        if orientation == 'tb':
            A[io:io+npmls,jo:ji] = np.repeat([np.flip(border_vec)],ji-jo,axis=0).T
            A[ii-npmls:ii,jo:ji] = np.repeat([border_vec],ji-jo,axis=0).T
        
        if not inline:
            return A
        
    def draw(self,T):
    	"""
    	def draw(T):
    	'draws' the coordinates extracted from the bitmap onto the field equations.
    	It draws a value of -1 for a perfect reflector for the ca field (ca (ex, ey) ) 
    	During the difference equation operating the FEM, this will flip the sign and send the wave in the opposite direction.
    	It draws 0 for (cb (ex, ey) ). Regardless of any change in the field. These will have a charge 0 which set our boundary
    	condition to the standing wave in the cavity between walls.
    	
    	These draw idealistic boundaries in a raster. Further improvements can be added by adding lossy materials, anisotropic materials, 
    	and materials of different dielectric and wave propogation constants. For now this is just a filler.
    	
    	INPUT 
    	T -> class containing the extracted coordinates of the BMP image.
    	"""
        # Establish materials
        coords = tuple(T.coords)
        coords_T = tuple(T.T())
        self.cE[:,:,0][coords] = -1 #caex
        self.cE[:,:,1][coords] = 0  #cbex
        self.cE[:,:,2][coords_T] = -1 #caey.T
        self.cE[:,:,3][coords_T] = 0  #cbey.T
        
        self.source = tuple(T.source)
        
    # function handles for pseudo pass by reference
    def ex(self):
        return self.E[:,:,0]
    
    def eyT(self):
        return self.E[:,:,1].T
    
    def hzx(self):
        return self.hz[:,:,1]
    
    def dhz_dx(self):
        i,j = self.theater
        return self.hz[:i,1:j,0] - self.hz[:i,:j-1,0]
    
    def hzy(self):
        return self.hz[:,:,2]
    
    def boundary_image(self):
        return self.cE[:,:100,3].T + self.cE[:,:100,1]
    
    
    
    
