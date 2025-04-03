#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 00:45:11 2025

@author: catherine
"""

"""RESEARCH Assignment 3

Question I will be answering:
How do the particles from the Milky Way and the M33 halo and 
the M31 halo contribute to the final shape of the remnant.


Will look at the triaxiality of the remnant, as well as see how the particles
from each original galaxy contribute to the triaxiality


First we will need to import modules"""
# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile

# for contours
import scipy.optimize as so

#photoutils
import photutils


#Then we will need to define the shape of the remnant in XY, YZ, XZ 
#Can do this by making histograms in each of these coordinates"""

#Using contour plot code from Lab 7, we can make histograms 
#first the code for plotting contours(modified from Lab 7)

def find_confidence_interval(x, pdf, confidence_level):
    return pdf[pdf > x].sum() - confidence_level

def density_contour(xdata, ydata, nbins_x, nbins_y, ax=None, **contour_kwargs):
    """ Create a density contour plot.
    Parameters
    ----------
    xdata : numpy.ndarray
    ydata : numpy.ndarray
    nbins_x : int
        Number of bins along x dimension
    nbins_y : int
        Number of bins along y dimension
    ax : matplotlib.Axes (optional)
        If supplied, plot the contour to this axis. Otherwise, open a new figure
    contour_kwargs : dict
        kwargs to be passed to pyplot.contour()
        
    Example Usage
    -------------
     density_contour(x pos, y pos, contour res, contour res, axis, colors for contours)
     e.g.:
     density_contour(xD, yD, 80, 80, ax=ax, 
         colors=['red','orange', 'yellow', 'orange', 'yellow'])

    """
    H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), density=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    Z = pdf.T
    fmt = {}
    
    ### Adjust Here #### 
    
    # Contour Levels Definitions
    one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
    two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
    three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
    
    # You might need to add a few levels
    my_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.75))


    # Array of Contour levels. Adjust according to the above
    levels = [one_sigma, my_sigma, two_sigma, three_sigma,][::-1]
    
    # contour level labels  Adjust accoding to the above.
    strs = ['0.68','0.75','0.95', '0.99'][::-1]
    ###### 
    
    if ax == None:
        contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
    
    return contour


# Create a COM of object for M31 Dark Matter (particle type=1) Using Code from Homework 4
COMDM31 = CenterOfMass("M31_000.txt",1)
# Create a COM of object for MW Dark Matter (particle type=1) Using Code from Homework 4
COMDMW = CenterOfMass("MW_000.txt",1)

# Compute COM of M31 using disk particles
COMPM31 = COMDM31.COM_P(0.1)
COMVM31 = COMDM31.COM_V(COMPM31[0],COMPM31[1],COMPM31[2])
# Compute COM of MW using disk particles
COMPMW = COMDMW.COM_P(0.1)
COMVMW = COMDMW.COM_V(COMPMW[0],COMPMW[1],COMPMW[2])

# Determine positions of disk particles relative to COM for M31
xDM31 = COMDM31.x - COMPM31[0].value 
yDM31 = COMDM31.y - COMPM31[1].value 
zDM31 = COMDM31.z - COMPM31[2].value 
# Determine positions of disk particles relative to COM for MW
xDMW = COMDMW.x - COMPMW[0].value 
yDMW = COMDMW.y - COMPMW[1].value 
zDMW = COMDMW.z - COMPMW[2].value 

# total magnitude for M31
rtotM31 = np.sqrt(xDM31**2 + yDM31**2 + zDM31**2)
# total magnitude for MW
rtotMW = np.sqrt(xDMW**2 + yDMW**2 + zDMW**2)

# Determine velocities of disk particles relative to COM motion for M31
vxDM31 = COMDM31.vx - COMVM31[0].value 
vyDM31 = COMDM31.vy - COMVM31[1].value 
vzDM31 = COMDM31.vz - COMVM31[2].value 
# Determine velocities of disk particles relative to COM motion for MW
vxDMW = COMDMW.vx - COMVMW[0].value 
vyDMW = COMDMW.vy - COMVMW[1].value 
vzDMW = COMDMW.vz - COMVMW[2].value 

# total velocity for M31
vtotM31 = np.sqrt(vxDM31**2 + vyDM31**2 + vzDM31**2)
# total velocity for MW
vtotMW = np.sqrt(vxDMW**2 + vyDMW**2 + vzDMW**2)

# Arrays for r and v for M31
rM31 = np.array([xDM31,yDM31,zDM31]).T # transposed 
vM31 = np.array([vxDM31,vyDM31,vzDM31]).T
# Arrays for r and v for MW
rMW = np.array([xDMW,yDMW,zDMW]).T # transposed 
vMW = np.array([vxDMW,vyDMW,vzDMW]).T

#code to rotate frame
def RotateFrame(posI,velI):
    """a function that will rotate the position and velocity vectors
    so that the disk angular momentum is aligned with z axis. 
    
    PARAMETERS
    ----------
        posI : `array of floats`
             3D array of positions (x,y,z)
        velI : `array of floats`
             3D array of velocities (vx,vy,vz)
             
    RETURNS
    -------
        pos: `array of floats`
            rotated 3D array of positions (x,y,z) 
            such that disk is in the XY plane
        vel: `array of floats`
            rotated 3D array of velocities (vx,vy,vz) 
            such that disk angular momentum vector
            is in the +z direction 
    """
    
    # compute the angular momentum
    L = np.sum(np.cross(posI,velI), axis=0)
    
    # normalize the angular momentum vector
    L_norm = L/np.sqrt(np.sum(L**2))


    # Set up rotation matrix to map L_norm to
    # z unit vector (disk in xy-plane)
    
    # z unit vector
    z_norm = np.array([0, 0, 1])
    
    # cross product between L and z
    vv = np.cross(L_norm, z_norm)
    s = np.sqrt(np.sum(vv**2))
    
    # dot product between L and z 
    c = np.dot(L_norm, z_norm)
    
    # rotation matrix
    I = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    v_x = np.array([[0, -vv[2], vv[1]], [vv[2], 0, -vv[0]], [-vv[1], vv[0], 0]])
    R = I + v_x + np.dot(v_x, v_x)*(1 - c)/s**2

    # Rotate coordinate system
    pos = np.dot(R, posI.T).T
    vel = np.dot(R, velI.T).T
    
    return pos, vel



#XY plane:    
#for M31
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def FindHistBounds(xIn, yIn, nbins_x=150,nbins_y=150, boundaryThreshold=0.01):
    H, xedges, yedges = np.histogram2d(xIn, yIn, bins=(nbins_x,nbins_y), density=True)
    # NOTE : if you are using the latest version of python, in the above: 
    # instead of normed=True, use density=True
    
    x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

    pdf = (H*(x_bin_sizes*y_bin_sizes))
    
# =============================================================================
#     X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
#     Z = pdf.T
#     fmt = {}
#     
#     ### Adjust Here #### 
#     # Contour Levels Definitions
#     one_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.68))
#     two_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.95))
#     three_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, 0.99))
# 
#     # Contour Levels Definitions
#     my_sigma = so.brentq(find_confidence_interval, 0., 1., args=(pdf, boundaryThreshold))
#     
# =============================================================================
    
    xmax = np.max(pdf,axis=1)
    ymax = np.max(pdf,axis=0)
    xregion, = np.nonzero(xmax > boundaryThreshold)
    xspan = xregion[-1]-xregion[0]
    xBdryIdx = [xregion[0]-xspan,xregion[-1]+xspan]
    xboundary = [xedges[xBdryIdx[0]],xedges[xBdryIdx[-1]]]
    
    yregion, = np.nonzero(ymax > boundaryThreshold)
    yspan = yregion[-1]-yregion[0]
    yBdryIdx = [yregion[0]-yspan,yregion[-1]+yspan]
    yboundary = [yedges[yBdryIdx[0]],yedges[yBdryIdx[-1]]]
    #===================================================================
    fig, ax= plt.subplots(figsize=(12, 10))
    ax.plot(xmax,'r.-')
    ax.plot(xBdryIdx[0],xmax[xBdryIdx[0]],'ro')
    ax.plot(xBdryIdx[-1],xmax[xBdryIdx[-1]],'ro')
    ax.plot(ymax,'b.-')
    ax.plot(yBdryIdx[0],ymax[yBdryIdx[0]],'bo')
    ax.plot(yBdryIdx[-1],ymax[yBdryIdx[-1]],'bo')
    #ax.set_ylim([0, 200])
    ax.grid()
    fig.show
    fig, ax= plt.subplots(figsize=(12, 10))
    #===================================================================
    return (xboundary,yboundary)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother

plt.hist2d(xDM31, yDM31, bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xDM31, yDM31, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('y (kpc) ', fontsize=22)

#set axis limits
#plt.ylim(-40,40)
#plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size
fig.show
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#************************************************************************
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother
xbdry, ybdry = FindHistBounds(xDM31, yDM31, 150)
# print("xbdry: ",round(xbdry),"  ybdry: ",round(ybdry))


plt.hist2d(xDM31, yDM31, bins=150, range=[xbdry, ybdry], norm=LogNorm(), cmap='magma')
#plt.hist2d(xDM31, yDM31, bins=150, norm=LogNorm(), cmap='magma')
cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xDM31, yDM31, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('y (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(ybdry)
# plt.xlim(xbdry)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


fig.show
# Save to a file
plt.savefig('ResearchAssignment3_M31XYHist.png')

#for MW
#************************************************************************
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother
xbdry, ybdry = FindHistBounds(xDM31, zDM31, 150)
# print("xbdry: ",round(xbdry),"  ybdry: ",round(ybdry))


plt.hist2d(xDM31, zDM31, bins=150, range=[xbdry, ybdry], norm=LogNorm(), cmap='magma')
#plt.hist2d(xDM31, yDM31, bins=150, norm=LogNorm(), cmap='magma')
cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xDM31, zDM31, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('z (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(ybdry)
# plt.xlim(xbdry)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


fig.show
# Save to a file
plt.savefig('ResearchAssignment3_M31XYHist.png')

#for MW
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother
plt.hist2d(xDMW, yDMW, bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for MW", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xDMW, yDMW, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('y (kpc) ', fontsize=22)

#set axis limits
xbdry, ybdry = FindHistBounds(xDMW, yDMW, 150)
plt.ylim(ybdry)
plt.xlim(xbdry)
# plt.ylim(-40,40)
# plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size



# Save to a file
plt.savefig('ResearchAssignment3_MW_XY_Hist.png')


    
#YZ plane:
#M31
rnM31, vnM31 = RotateFrame(rM31,vM31)
# Rotated M31 Disk - EDGE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(15, 10))

# plot the particle density for M31 , 2D histogram
# ADD HERE
plt.hist2d(rnM31[:,1], rnM31[:,2], bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin", fontsize=15)


# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('z (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(-10,10)
# plt.xlim(-45,45)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

density_contour(rnM31[:,1], rnM31[:,2], 80, 80, ax=ax, colors=['red','blue','green','yellow'])

# Save to a file
plt.savefig('ResearchAssignment_M31_YZ_Hist.png')

#MW
rnMW, vnMW = RotateFrame(rMW,vMW)
# Rotated M31 Disk - EDGE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(15, 10))

# plot the particle density for M31 , 2D histogram
# ADD HERE
plt.hist2d(rnM31[:,1], rnM31[:,2], bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin", fontsize=15)


# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('z (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(-10,10)
# plt.xlim(-45,45)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

density_contour(rnMW[:,1], rnMW[:,2], 80, 80, ax=ax, colors=['red','blue','green','yellow'])

# Save to a file
plt.savefig('ResearchAssignment_M31_Xy_Hist.png')

    
#XZ plane:
#M31
rnM31, vnM31 = RotateFrame(rM31,vM31)
# Rotated M31 Disk - EDGE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(15, 10))

# plot the particle density for M31 , 2D histogram
# ADD HERE
plt.hist2d(rnM31[:,0], rnM31[:,2], bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin", fontsize=15)


# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('z (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(-10,10)
# plt.xlim(-45,45)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

density_contour(rnM31[:,0], rnM31[:,2], 80, 80, ax=ax, colors=['red','blue','green','yellow'])

# Save to a file
plt.savefig('ResearchAssignment_M31_XZ_Hist.png')
    
#MW
rnMW, vnMW = RotateFrame(rMW,vMW)
# Rotated M31 Disk - EDGE ON

# M31 Disk Density 
fig, ax= plt.subplots(figsize=(15, 10))

# plot the particle density for M31 , 2D histogram
# ADD HERE
plt.hist2d(rnMW[:,0], rnMW[:,2], bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin", fontsize=15)


# Add axis labels
plt.xlabel('x (kpc)', fontsize=22)
plt.ylabel('z (kpc) ', fontsize=22)

#set axis limits
# plt.ylim(-10,10)
# plt.xlim(-45,45)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

density_contour(rnMW[:,0], rnMW[:,2], 80, 80, ax=ax, colors=['red','blue','green','yellow'])

# Save to a file
plt.savefig('ResearchAssignment_MW_XZ_Hist.png')    


#We can then make contour plots to understand the density in the histograms


#We can fit ellipses to these different histograms and decide where the density 
#determines the ellipses of the that specific coordinate orientation


#We can then use the ellipses that we determine, and examine which galaxy contributes
#the most to each apsect of triaxiality of the final remnant



#First we run these functiions with the final remnant(MW +M31) then we can lrun the functions
#on MW and M31 separately and compare



