#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 00:45:11 2025

@author: catherine
"""

"""RESEARCH Assignment 3

Research Topic:
Dark matter halo (density
or halo shape) after a major merger.

Specific question: 
What is the contribution of the MW vs. M31 halo particles to the shape of the merged remnant?
This code will look at where the specific halo
particles from each galaxy end up coalescing in the final remnant.

Will plot the ratio of a specific galaxy's particles to the complete number of particles



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
COMDM31 = CenterOfMass("M31_610.txt",1)
# Create a COM of object for MW Dark Matter (particle type=1) Using Code from Homework 4
COMDMW = CenterOfMass("MW_610.txt",1)

# Compute COM of M31 using disk particles
COMPM31 = COMDM31.COM_P(0.1)
#COMVM31 = COMDM31.COM_V(COMPM31[0],COMPM31[1],COMPM31[2])
# Compute COM of MW using disk particles
COMPMW = COMDMW.COM_P(0.1)
#COMVMW = COMDMW.COM_V(COMPMW[0],COMPMW[1],COMPMW[2])

#Crete a COM of the entire system
COMSYS = (COMPM31 + COMPMW)/2

# Determine positions of disk particles relative to COM for M31
xDM31 = COMDM31.x - COMSYS[0].value 
yDM31 = COMDM31.y - COMSYS[1].value 
zDM31 = COMDM31.z - COMSYS[2].value 
# Determine positions of disk particles relative to COM for MW
xDMW = COMDMW.x - COMSYS[0].value 
yDMW = COMDMW.y - COMSYS[1].value 
zDMW = COMDMW.z - COMSYS[2].value 

#Use numpy's hstack to concatinate positions of particles for entire system
xSYS = np.hstack((xDM31,xDMW))
ySYS = np.hstack((yDM31,yDMW))
zSYS = np.hstack((zDM31,zDMW))


#Create plot for XY plane   
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother
plt.hist2d(xSYS, ySYS, bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xSYS, ySYS, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('y (kpc) ', fontsize=22)

plt.title('Combined M31 - MW System')
#set axis limits
# plt.ylim(-40,40)
# plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size



# Save to a file
plt.savefig('ResearchAssignment3_M31XYHist.png')
#========================================================================




#Create histogram of the percentage of each particle. M31 in red and MW in blue
Hsys, xedges, yedges = np.histogram2d(xSYS, ySYS, bins=150, density=False)
Hsm = np.ma.masked_where(Hsys == 0, Hsys)
X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
Hm31, xedges, yedges = np.histogram2d(xDM31, yDM31, bins=(xedges, yedges),
                                      density=False)
# x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
# y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
Grat = Hm31/Hsys


fig, ax = plt.subplots(figsize=(12, 10))
sysimg = ax.pcolormesh(X, Y, Hsys.T, norm=LogNorm(), cmap='magma')
# ax.pcolormesh(X, Y, Hsys, cmap='magma')
cbar = fig.colorbar(sysimg,ax=ax)
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)
# Add axis labels
ax.set_xlabel('x (kpc) ', fontsize=22)
ax.set_ylabel('y (kpc) ', fontsize=22)
plt.title('Combined M31 - MW System')
fig, ax = plt.subplots(figsize=(12, 10))
# ax.pcolormesh(X, Y, Grat, cmap='bwr')
# gplt = ax.pcolormesh(X, Y, Grat, cmap='plasma')
gplt = ax.pcolormesh(X, Y, Grat.T, cmap='rainbow')
cbar = fig.colorbar(gplt,ax=ax)
cbar.set_label("Fraction of disk particle that are M31", fontsize=15)
# Add axis labels
ax.set_xlabel('x (kpc) ', fontsize=22)
ax.set_ylabel('y (kpc) ', fontsize=22)
plt.title('Combined M31 - MW System Ratio')
plt.show()
#========================================================================
#========================================================================

# =============================================================================
#=======================================================================
#=======================================================================

#XY plane:    
#for M31
cutSysIdx = np.nonzero(zSYS<0)
cutDM31Idx = np.nonzero(zDM31<0)
fig, ax= plt.subplots(figsize=(12, 10))

# ADD HERE
# plot the particle density for M31 using a 2D historgram
# plt.hist2D(pos1,pos2, bins=, norm=LogNorm(), cmap='' )
# cmap options: 
# https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html  
#   e.g. 'magma', 'viridis'
# can modify bin number to make the plot smoother
plt.hist2d(xSYS[cutSysIdx], ySYS[cutSysIdx], bins=150, norm=LogNorm(), cmap='magma')

cbar = plt.colorbar()
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

# ADD HERE
# make the contour plot
# x pos, y pos, contour res, contour res, axis, colors for contours.
# remember to adjust this if there are other contours added
# density_contour(pos1, pos2, res1, res2, ax=ax, colors=[])
density_contour(xSYS, ySYS, 80, 80, ax=ax, colors=['red','blue','green','yellow'])


# Add axis labels
plt.xlabel('x (kpc) ', fontsize=22)
plt.ylabel('y (kpc) ', fontsize=22)
plt.title('Combined M31 - MW System (Cutaway)')
#set axis limits
# plt.ylim(-40,40)
# plt.xlim(-40,40)

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size



# Save to a file
plt.savefig('ResearchAssignment3_CombinedXYHist.png')
#========================================================================
#========================================================================
Hsys, xedges, yedges = np.histogram2d(xSYS[cutSysIdx], ySYS[cutSysIdx], bins=150, density=False)
Hsm = np.ma.masked_where(Hsys == 0, Hsys)
X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
Hm31, xedges, yedges = np.histogram2d(xDM31[cutDM31Idx], yDM31[cutDM31Idx], bins=(xedges, yedges),
                                      density=False)
# x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
# y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
Grat = Hm31/Hsys


fig, ax = plt.subplots(figsize=(12, 10))
sysimg = ax.pcolormesh(X, Y, Hsys.T, norm=LogNorm(), cmap='magma')
# ax.pcolormesh(X, Y, Hsys, cmap='magma')
cbar = fig.colorbar(sysimg,ax=ax)
cbar.set_label("Number of disk particle per bin for M31", fontsize=15)
# Add axis labels
ax.set_xlabel('x (kpc) ', fontsize=22)
ax.set_ylabel('y (kpc) ', fontsize=22)
plt.title('Combined M31 - MW System (Cutaway)')
fig, ax = plt.subplots(figsize=(12, 10))
# ax.pcolormesh(X, Y, Grat, cmap='bwr')
# gplt = ax.pcolormesh(X, Y, Grat, cmap='plasma')
gplt = ax.pcolormesh(X, Y, Grat.T, cmap='rainbow')
cbar = fig.colorbar(gplt,ax=ax)
cbar.set_label("Fraction of disk particle that are M31", fontsize=15)
# Add axis labels
ax.set_xlabel('x (kpc) ', fontsize=22)
ax.set_ylabel('y (kpc) ', fontsize=22)
plt.title('Combined M31 - MW System Ratio (Cutaway)')
plt.show()

