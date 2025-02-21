#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: catherinesnyder
"""

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
    
    def __init__(self, galaxy, snap):
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
        self.filename='%s_'%(galaxy) + ilbl + ".txt"
        self.time, self.total, self.data = Read(self.filename) 
        self.gname = galaxy
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)
        
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
        self.COM = CenterOfMass(self.filename, 2)
        COMpos = self.COM.COM_P(0.1)/u.kpc
        X_COM = COMpos[0]
        Y_COM = COMpos[1]
        Z_COM = COMpos[2]
        #radius_array = len(radius_array)
        self.index = np.where(self.data['type'] == particle_type)
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index] #store x pos
        self.y = self.data['y'][self.index] #store y pos
        self.z = self.data['z'][self.index] #store z pos
        x_new = self.x - X_COM
        y_new = self.y - Y_COM
        z_new = self.z - Z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2) #in kpc
        #index2 = np.where(r_new <= r_max) #index of particles in range
        #x2 = self.x[index2] #index x values in range
        #y2 = self.y[index2] #index y values in range
        #z2 = self.z[index2] #index z values in range
        #m2 = self.m[index2] #index m values in range
        mass_array = np.zeros(len(radius_array))
        
        for i in range(len(radius_array)):
            index2 = np.where(r_new <= radius_array[i])#index of particles in range
            if len(index2)>0:
                mass_sum = sum(self.m[index2])
                mass_array[i] = mass_sum
        return mass_array*1e10*u.Msun
    
        
    #3
    def MassEnclosedTotal(self, radii):
        '''
        Function to return an arry of masses representing the total mass enclosed at
        each radius
        Inputs: radii : array of radii
        Outputs: total_enclosed_mass = mass enclosed of all particles in particular radii
        '''
        
        self.MassBulge = self.MassEnclosed(3, radii)
        self.MassHalo = self.MassEnclosed(1, radii)
        self.MassDisk = self.MassEnclosed(2, radii)
        self.total_enclosed_mass = self.MassBulge + self.MassHalo + self.MassDisk
        return self.total_enclosed_mass
    #4    
    def HernquistMass(radius, h_a=60*u.kpc, Mhalo=1.975):
        '''
        Function that will compute the mass enclosed within a given radius 
        using the theoretical profile
        Inputs: radius = radius for calculation
                a = scale factor
                Mhalo = halo mass
        Outputs: hernquist_mass_profile = halo mass in solar mass units
        '''
        a = Mhalo*1e12*u.Msun #constants, correcting units
        b = radius**2/(h_a +radius)**2

    
        hernquist_mass_profile = a*b #   
    
        return hernquist_mass_profile
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
        M = self.MassEnclosed(particle_type,radii)
        circular_velocity = np.zeros(len(M))
        for i in range(len(M)):
            circular = np.sqrt(G*M[i]/radii[i])
            circular_velocity[i] = circular
        return circular_velocity
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
        M = self.MassEnclosedTotal(radii)
        total_circular_velocity = np.zeros(len(M))
        for i in range(len(M)):
            circular = np.sqrt(G*M[i]/radii[i])
            total_circular_velocity[i] = circular
        return total_circular_velocity
        
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
        M = self.HernquistMass(radius, h_a=60*u.kpc, Mhalo =1.975)
        hernquist_circular_velocity = np.sqrt(G*M/radius)
        return hernquist_circular_velocity
        
   
 #8
 # Plot the Mass Profile for each galaxy
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    #MW plot
    MW_MP = MassProfile("MW", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = MW_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, massHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    massDisk = MW_MP.MassEnclosed(2, r)
    ax.plot( r, massDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    massBulge = MW_MP.MassEnclosed(3, r)
    ax.plot( r, massBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    massTotal = MW_MP.MassEnclosedTotal(r)
    ax.plot( r, massTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('MW Component Mass')
    ax.legend()
    plt.show() 
    plt.savefig('MW Component Mass.png')
    
    #M33 plot
    M33_MP = MassProfile("M33", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = M33_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, massHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    massDisk = M33_MP.MassEnclosed(2, r)
    ax.plot( r, massDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    massBulge = M33_MP.MassEnclosed(3, r)
    ax.plot( r, massBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    massTotal = M33_MP.MassEnclosedTotal(r)
    ax.plot( r, massTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('M33 Component Mass')
    ax.legend()
    plt.show() 
    plt.savefig('M33 Component Mass.png')

    #M31 plot
    M31_MP = MassProfile("M31", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    massHalo = M31_MP.MassEnclosed(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, massHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    massDisk = M31_MP.MassEnclosed(2, r)
    ax.plot( r, massDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    massBulge = M31_MP.MassEnclosed(3, r)
    ax.plot( r, massBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    massTotal = M31_MP.MassEnclosedTotal(r)
    ax.plot( r, massTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('M31 Component Mass')
    ax.legend()
    plt.show() 
    plt.savefig('M31 Component Mass.png')

#9
#Plot the Rotation Curve for each galaxy
    #MW plot
    MW_MP = MassProfile("MW", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = MW_MP.CircularVelocity(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, velHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    velDisk = MW_MP.CircularVelocity(2, r)
    ax.plot( r, velDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    velBulge = MW_MP.CircularVelocity(3, r)
    ax.plot( r, velBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    velTotal = MW_MP.CircularVelocityTotal(r)
    ax.plot( r, velTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('MW Rotation Curve')
    ax.legend()
    plt.show() 
    plt.savefig('MW Rotation Curve.png')
    
    #M33 plot
    
    M33_MP = MassProfile("M33", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = M33_MP.CircularVelocity(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, velHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    velDisk = M33_MP.CircularVelocity(2, r)
    ax.plot( r, velDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    velBulge = M33_MP.CircularVelocity(3, r)
    ax.plot( r, velBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    velTotal = M33_MP.CircularVelocityTotal(r)
    ax.plot( r, velTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('M33 Rotation Curve')
    ax.legend()
    plt.show() 
    plt.savefig('M33 Rotation Curve.png')
        
    
    #M31 plot
    M31_MP = MassProfile("MW", 0)
    
    r = np.arange(0.1, 30.5, 1.5); print(r) # create an array of radii as the input

    velHalo = M31_MP.CircularVelocity(1, r) # get the enclosed halo masses at each element in 'r
    #marry = MW_MP.MassEnclosed(
    fig, ax = plt.subplots()             # Create a figure containing a single Axes.
    ax.plot( r, velHalo, 'ro-', label='Halo')  # Plot some data on the Axes.
    velDisk = M31_MP.CircularVelocity(2, r)
    ax.plot( r, velDisk, 'bx-', label='Disc')  # Plot some data on the Axes.
    velBulge = M31_MP.CircularVelocity(3, r)
    ax.plot( r, velBulge, 'mv-', label='Bulge')  # Plot some data on the Axes.
    velTotal = M31_MP.CircularVelocityTotal(r)
    ax.plot( r, velTotal, 'g^-', label='Total')  # Plot some data on the Axes.
    ax.grid(visible=True)
    ax.set_xlabel('Radius [kpc]')
    ax.set_ylabel('Enclosed Mass [Msun]')
    ax.set_title('M31 Rotation Curve')
    ax.legend()
    plt.show() 
    plt.savefig('M31 Rotation Curve.png')
    

        
        
        
        
