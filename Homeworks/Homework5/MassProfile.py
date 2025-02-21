#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: catherinesnyder
"""
#imports all neccesarry things
import numpy as np
import astropy.units as u
import astropy.table as tbl
 

from ReadFile import Read
from CenterOfMass import CenterOfMass
from astropy.constants import G 
import matplotlib.pyplot as plt
import matplotlib
#%matplotlib inline


class MassProfile:
    
    def __init__(self, galaxy, snap):#inititiaes class
        ''' 
            
            PARAMETERS
            ----------
            galaxy : `str`
                name of the galaxy being profiled
            snap : `int; 1, 2, or 3`
                snapshot number
        '''
        
        # add a string of the filenumber to the value “000”
        ilbl = "000" + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename='%s_'%(galaxy) + ilbl + ".txt" #creates filename
        self.time, self.total, self.data = Read(self.filename) #reads data from file
        self.gname = galaxy #names galaxy
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #converts G units
        
    #2    
    def MassEnclosed(self, particle_type, radius_array): 
        '''
        Function to compute the mass enclosed within given radius of COM for
        specific galaxy and specific particle type
        Inputs: particle_type = specified in the file, possible particle types are Dark Matter,
                       Disk Stars, or Bulge Stars
                radius = array of radii
        Outputs: enclosed_mass = mass enclosed of particular particle in particular radii
            
        '''
        self.COM = CenterOfMass(self.filename, 2) #creates COM and imports function
        COMpos = self.COM.COM_P(0.1)/u.kpc #Calls on COM function and gets units
        X_COM = COMpos[0] #gets the COM x pos
        Y_COM = COMpos[1] #gets the COM y pos
        Z_COM = COMpos[2] #gets the COM z pos
        #radius_array = len(radius_array)
        self.index = np.where(self.data['type'] == particle_type) #makes index of particle type
        self.m = self.data['m'][self.index] #store mass
        self.x = self.data['x'][self.index] #store x pos
        self.y = self.data['y'][self.index] #store y pos
        self.z = self.data['z'][self.index] #store z pos
        x_new = self.x - X_COM #get proper x coordinate
        y_new = self.y - Y_COM #get proper y coordinate
        z_new = self.z - Z_COM #get proper z coordinate
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) # gets the raius in kpc
        #index2 = np.where(r_new <= r_max) #index of particles in range
        #x2 = self.x[index2] #index x values in range
        #y2 = self.y[index2] #index y values in range
        #z2 = self.z[index2] #index z values in range
        #m2 = self.m[index2] #index m values in range
        mass_array = np.zeros(len(radius_array)) #creates empty mass array
        
        for i in range(len(radius_array)):#goes through length of radius array
            index2 = np.where(r_new <= radius_array[i])#index of particles in range
            if len(index2)>0: #checks to make sure there are particles of that type
                mass_sum = sum(self.m[index2]) #sums the masses of particles in that radius
                mass_array[i] = mass_sum #adds the mass to array of masses
        return mass_array*1e10*u.Msun #returns mass array with units
    
        
    #3
    def MassEnclosedTotal(self, radii):
        '''
        Function to return an arry of masses representing the total mass enclosed at
        each radius
        Inputs: radii : array of radii
        Outputs: total_enclosed_mass = mass enclosed of all particles in particular radii
        '''
        
        self.MassBulge = self.MassEnclosed(3, radii) #gets mass array for bulge particles
        self.MassHalo = self.MassEnclosed(1, radii) #gets mass array for halo particles
        self.MassDisk = self.MassEnclosed(2, radii) #gets mass array for disk particles
        self.total_enclosed_mass = self.MassBulge + self.MassHalo + self.MassDisk #adds all the masses of each particle type
        return self.total_enclosed_mass #returns total mass array
    #4    
    def HernquistMass(self, radius, h_a=60*u.kpc, Mhalo=1.975):
        '''
        Function that will compute the mass enclosed within a given radius 
        using the theoretical profile
        Inputs: radius = radius for calculation
                a = scale factor
                Mhalo = halo mass
        Outputs: hernquist_mass_profile = halo mass in solar mass units
        '''
        a = Mhalo*1e12*u.Msun #constants, correcting units
        b = radius**2/(h_a +radius)**2 #creates b

    
        hernquist_mass_profile = a*b #finds the mass profile   
    
        return hernquist_mass_profile #returns mass profile
    #5    
    def CircularVelocity(self, particle_type, radii):
        '''
        Function that finds the circular speed using mass enclosed by radii 
        assuming spherical symmetry
        Inputs: particle_type = specified in the file, possible particle types are Dark Matter,
                       Disk Stars, or Bulge Stars
                radius = array of radii
        Outputs: circular_velocity = array of circular speeds in km/s
        '''
        M = self.MassEnclosed(particle_type,radii) #gets mass array for particle type
        circular_velocity = np.zeros(len(M))*u.km/u.s #creates empty velocity array
        for i in range(len(M)): #iterates the length of the mass index
            circular = np.sqrt(self.G*M[i]/(radii[i]*u.kpc)) #calculates circular velocity for radius
            circular_velocity[i] = circular #adds the circular veloocity to the array
        return circular_velocity #returns circular velocity ray
    #6    
    def CircularVelocityTotal(self, radii):
        '''
        Function that finds the total circular velocity of the galaxy with all
        particle types
        Inputs: radii = array of radii
        Outputs: total_circular_velocity = array of circular velocity (in units of km/s) representing the total Vcirc
        created by all the galaxy components (bulge+disk+halo) 
        at each radius of the input array
        '''
        M = self.MassEnclosedTotal(radii) #gets total mass array
        total_circular_velocity = np.zeros(len(M))*u.km/u.s #creates empty velocity array with units
        for i in range(len(M)): #iterates through lenght of mass array
            circular = np.sqrt(self.G*M[i]/(radii[i]*u.kpc)) #calculates the circular velocity
            total_circular_velocity[i] = circular #adds total circular velocity to array
        return total_circular_velocity #returns array
        
    #7
    def HernquistVCirc(self, radius, h_a=60*u.kpc, Mhalo=1.975):
        '''
        Function that finds the circular speed using the Hernquist mass profile
        Inputs: radius = radius for calculation
                a = scale factor
                Mhalo = halo mass
        Outputs: hernquist_circular_velocity =  the circular speed 
        in units of km/s, rounded to two decimal places.
        '''
        M = self.HernquistMass(radius, h_a=60*u.kpc, Mhalo =1.975) #calls on the hernquist mass profile
        hernquist_circular_velocity = np.sqrt(G*M/radius) #calculates circular velocity
        return hernquist_circular_velocity #returns circular velocity
        
   
 #8
 # Plot the Mass Profile for each galaxy
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    #MW plot
    MW_MP = MassProfile("MW", 0) #creates mass profile for Milky Way
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = MW_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, massHalo, 'ro-', label='Halo')  # Plot halo data on the Axes.
    massDisk = MW_MP.MassEnclosed(2, r) # get the enclosed disk masses at each element in 'r
    ax.semilogy( r, massDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    massBulge = MW_MP.MassEnclosed(3, r) # get the enclosed bulge masses at each element in 'r
    ax.semilogy( r, massBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    massTotal = MW_MP.MassEnclosedTotal(r) # get the total enclosed masses at each element in 'r
    ax.semilogy( r, massTotal, 'g^-', label='Total')# Plot totals data on the Axes.
    
    massHern = MW_MP.HernquistMass(r*u.kpc, h_a=60*u.kpc, Mhalo=1.975) # get the hernquist masses at each element in 'r
    ax.semilogy( r, massHern/u.Msun, 'ks-', label='Hernquist') #Plot hernquist data on the Axes
    print('Hernquist Mass:\n',massHern) #print the hernquist
    
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) #makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Mass [Msun]') #label y value
    ax.set_title('MW Component Mass') #label title
    ax.legend() #creates a legend
    plt.show() #show plot
    plt.savefig('MW Component Mass.png') #save the figure to a png
    
    #M33 plot
    M33_MP = MassProfile("M33", 0) #creates mass profile for M33
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = M33_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, massHalo, 'ro-', label='Halo')  # Plot halo data on the Axes.
    massDisk = M33_MP.MassEnclosed(2, r)  # get the enclosed disk masses at each element in 'r
    ax.semilogy( r, massDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    massBulge = M33_MP.MassEnclosed(3, r) # get the enclosed bulge masses at each element in 'r
    ax.semilogy( r, massBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    massTotal = M33_MP.MassEnclosedTotal(r) # get the enclosed total masses at each element in 'r
    ax.semilogy( r, massTotal, 'g^-', label='Total')# Plot total data on the Axes.
    massHern = MW_MP.HernquistMass(r*u.kpc, h_a=60*u.kpc, Mhalo=1.975) # get the hernquist masses at each element in 'r
    ax.semilogy( r, massHern/u.Msun, 'ks-', label='Hernquist')  #Plot hernquist data on the Axes
    print('Hernquist Mass:\n',massHern) #print the hernquist
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) #makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Mass [Msun]') #label y value
    ax.set_title('M33 Component Mass') #label a title
    ax.legend() #creates a legend
    plt.show() #show plot
    fig.savefig('M33 Component Mass.png') #save plot as png

    #M31 plot
    M31_MP = MassProfile("M31", 0) ##creates mass profile for M33
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = M31_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, massHalo, 'ro-', label='Halo') # Plot halo data on the Axes.
    massDisk = M31_MP.MassEnclosed(2, r) # get the enclosed disk masses at each element in 'r
    ax.semilogy( r, massDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    massBulge = M31_MP.MassEnclosed(3, r) # get the enclosed bulge masses at each element in 'r
    ax.semilogy( r, massBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    massTotal = M31_MP.MassEnclosedTotal(r) # get the enclosed total masses at each element in 'r
    ax.semilogy( r, massTotal, 'g^-', label='Total')  # Plot total data on the Axes.
    massHern = MW_MP.HernquistMass(r*u.kpc, h_a=60*u.kpc, Mhalo=1.975) # get the hernquist masses at each element in 'r
    ax.semilogy( r, massHern/u.Msun, 'ks-', label='Hernquist') #Plot hernquist data on the Axes
    print('Hernquist Mass:\n',massHern) #print the hernquist
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) # makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Mass [Msun]') #label y value
    ax.set_title('M31 Component Mass') #label title
    ax.legend() #create a legend
    plt.show() #show plot
    fig.savefig('M31 Component Mass.png') #save plot as png

#9
#Plot the Rotation Curve for each galaxy
    #MW plot
    #MW_MP = MassProfile("MW", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = MW_MP.CircularVelocity(1, r) # get the enclosed halo vel at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, velHalo, 'ro-', label='Halo')  # Plot halo data on the Axes.
    velDisk = MW_MP.CircularVelocity(2, r) # get the enclosed disk vel at each element in 'r
    ax.semilogy( r, velDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    velBulge = MW_MP.CircularVelocity(3, r) # get the enclosed bulge vel at each element in 'r
    ax.semilogy( r, velBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    velTotal = MW_MP.CircularVelocityTotal(r) # get the enclosed total vel at each element in 'r
    ax.semilogy( r, velTotal, 'g^-', label='Total')  # Plot total data on the Axes.
    
    velHern = MW_MP.HernquistVCirc(r, h_a=60*u.kpc, Mhalo=1.975) # get the enclosed hernquist vel at each element in 'r
    ax.semilogy( r, massTotal, 'ks-', label='Hernquist')  # Plot hernquist data on the Axes.
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) #makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Vel [Msun]') #label y value
    ax.set_title('MW Rotation Curve') #label title
    ax.legend() #makes a legend
    plt.show() #show plot
    fig.savefig('MW Rotation Curve.png') #saves plot as png
    
    #M33 plot
    
    #M33_MP = MassProfile("M33", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = M33_MP.CircularVelocity(1, r) # get the enclosed halo vel at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, velHalo, 'ro-', label='Halo')  # Plot halo data on the Axes.
    velDisk = M33_MP.CircularVelocity(2, r) # get the enclosed disk vel at each element in 'r
    ax.semilogy( r, velDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    velBulge = M33_MP.CircularVelocity(3, r) # get the enclosed bulge vel at each element in 'r
    ax.semilogy( r, velBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    velTotal = M33_MP.CircularVelocityTotal(r) # get the enclosed total vel at each element in 'r
    ax.semilogy( r, velTotal, 'g^-', label='Total')  # Plot total data on the Axes.
    
    velHern = M33_MP.HernquistVCirc(r, h_a=60*u.kpc, Mhalo=1.975) # get the enclosed hernquist vel at each element in 'r
    ax.semilogy( r, massTotal, 'ks-', label='Hernquist') # Plot hernquist data on the Axes.
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) #makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Vel [Msun]') #label y value
    ax.set_title('M33 Rotation Curve') #label title
    ax.legend() #makes a legend
    plt.show() #show plot
    fig.savefig('M33 Rotation Curve.png') #saves plot as png
        
    
    #M31 plot
    #M31_MP = MassProfile("MW", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = M31_MP.CircularVelocity(1, r) # get the enclosed halo vel at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.semilogy( r, velHalo, 'ro-', label='Halo')  # Plot halo data on the Axes.
    velDisk = M31_MP.CircularVelocity(2, r) # get the enclosed disk vel at each element in 'r
    ax.semilogy( r, velDisk, 'bx-', label='Disc')  # Plot disk data on the Axes.
    velBulge = M31_MP.CircularVelocity(3, r) # get the enclosed bulge vel at each element in 'r
    ax.semilogy( r, velBulge, 'mv-', label='Bulge')  # Plot bulge data on the Axes.
    velTotal = M31_MP.CircularVelocityTotal(r) # get the enclosed total vel at each element in 'r
    ax.semilogy( r, velTotal, 'g^-', label='Total')  # Plot total data on the Axes.
    
    velHern = M31_MP.HernquistVCirc(r, h_a=60*u.kpc, Mhalo=1.975) # get the enclosed hernquist vel at each element in 'r
    ax.semilogy( r, massTotal, 'ks-', label='Hernquist') # Plot hernquist data on the Axes.
    ax.set_ylim(bottom = 10**9) #limits y
    ax.grid(visible=True) #makes a grid
    ax.set_xlabel('Radius [kpc]') #label x value
    ax.set_ylabel('Enclosed Vel [Msun]') #label y value
    ax.set_title('M31 Rotation Curve') #label title
    ax.legend() #make a legend
    plt.show()  #show plot
    fig.savefig('M31 Rotation Curve.png') #saves plot as png
    

        
        
        
        
