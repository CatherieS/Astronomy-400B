#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 00:45:11 2025

@author: catherine
"""

"""RESEARCH Assignment 3

Question I will be answering:
How do the particles from the Milky Way and the M31 halo contribute to the final remnant.



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
from ConcatenateFiles import concatenate_snap_files

# for contours
import scipy.optimize as so

#photoutils
import photutils


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
        plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=8)

    else:
        contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
        for l, s in zip(contour.levels, strs):
            fmt[l] = s
        ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=8)
    
    return contour
#==============================================================================


halfRange = 1000
bincnt = 51

merged_distance_threshold = 2
ConcatFileNeededFlag = False
snapDir = "C:\\Users\\catherine\\Anaconda\\ASTR400B\\Research\\Research7\\snapshots\\"



 
snapList = np.arange(235,736,100)

figS, axsS = plt.subplots(3, 2, figsize=(8, 10),layout="constrained")
figR, axsR = plt.subplots(3, 2, figsize=(8, 10),layout="constrained")
figG, axsG = plt.subplots(3, 2, figsize=(8, 10),layout="constrained")
figH, axsH = plt.subplots(3, 2, figsize=(8, 10),layout="constrained")
figN, axsN = plt.subplots(3, 2, figsize=(8, 10),layout="constrained")
totalParticleCnt=np.zeros((snapList.size,2))
COMoffset=np.zeros((snapList.size,2))
COMsepStore=np.zeros((snapList.size,1))

for i, snap in enumerate(snapList):
    infile1 = snapDir+'MW_'+str(snap).zfill(3)+'.txt'
    infile2 = snapDir+'M31_'+str(snap).zfill(3)+'.txt'
    concfile = snapDir+'CONC_'+str(snap).zfill(3)+'.txt'
    
    # Create a COM of object for M31 Dark Matter (particle type=1) Using Code from Homework 4
    COMDM31 = CenterOfMass(infile2,1)
    # Create a COM of object for MW Dark Matter (particle type=1) Using Code from Homework 4
    COMDMW = CenterOfMass(infile1,1)
    
    # Compute COM of M31 using disk particles
    COMPM31 = COMDM31.COM_P(0.1)
    COMVM31 = COMDM31.COM_V(COMPM31[0],COMPM31[1],COMPM31[2])
    # Compute COM of MW using disk particles
    COMPMW = COMDMW.COM_P(0.1)
    COMVMW = COMDMW.COM_V(COMPMW[0],COMPMW[1],COMPMW[2])
    
    # Calculate COM separation
    COMsep = np.sqrt(np.sum((COMPM31-COMPMW)**2))
    COMsepStore[i] = COMsep
    
    if False:
    # if COMsep.value < merged_distance_threshold:
        if ConcatFileNeededFlag:
            concatenate_snap_files(infile1, infile2, concfile)
        COMDSYS = CenterOfMass(concfile,1)
        COMPSYS = COMDSYS.COM_P(0.1)
        COMVSYS = COMDSYS.COM_V(COMPSYS[0],COMPSYS[1],COMPSYS[2])
    else:
        COMPSYS = (COMPM31 + COMPMW)/2
        COMVSYS = (COMVM31 + COMVMW)/2
        
        
    COMoffset[i,:] =np.array([np.sqrt(np.sum((COMPM31.value-COMPSYS.value)**2)),
                     np.sqrt(np.sum((COMPMW.value-COMPSYS.value)**2))])
    
    # Determine positions of disk particles relative to COM for M31
    xDM31 = COMDM31.x - COMPSYS[0].value 
    yDM31 = COMDM31.y - COMPSYS[1].value 
    zDM31 = COMDM31.z - COMPSYS[2].value 
    # Determine positions of disk particles relative to COM for MW
    xDMW = COMDMW.x - COMPSYS[0].value 
    yDMW = COMDMW.y - COMPSYS[1].value 
    zDMW = COMDMW.z - COMPSYS[2].value 
      
    # Determine positions of disk particles relative to COM for combined system
    xSYS = np.hstack((xDM31,xDMW))
    ySYS = np.hstack((yDM31,yDMW))
    zSYS = np.hstack((zDM31,zDMW))
    
# =============================================================================
# =============================================================================

    # total magnitude for M31
    rtotM31 = np.sqrt(xDM31**2 + yDM31**2 + zDM31**2)
    # total magnitude for MW
    rtotMW = np.sqrt(xDMW**2 + yDMW**2 + zDMW**2)
    # total magnitude for combined system
    rtotSYS = np.sqrt(xSYS**2 + ySYS**2 + zSYS**2)
    
    totalParticleCnt[i,:] = [rtotM31.size, rtotMW.size]
    
    # Determine velocities of disk particles relative to COM motion for M31
    vxDM31 = COMDM31.vx - COMVSYS[0].value 
    vyDM31 = COMDM31.vy - COMVSYS[1].value 
    vzDM31 = COMDM31.vz - COMVSYS[2].value 
    # Determine velocities of disk particles relative to COM motion for MW
    vxDMW = COMDMW.vx - COMVSYS[0].value 
    vyDMW = COMDMW.vy - COMVSYS[1].value 
    vzDMW = COMDMW.vz - COMVSYS[2].value 
    # Determine velocities of disk particles relative to COM motion for SYS
    
    vxSYS = np.hstack((vxDM31,vxDMW)) 
    vySYS = np.hstack((vyDM31,vyDMW))
    vzSYS = np.hstack((vzDM31,vzDMW))
    
    # total velocity for M31
    vtotM31 = np.sqrt(vxDM31**2 + vyDM31**2 + vzDM31**2)
    # total velocity for MW
    vtotMW = np.sqrt(vxDMW**2 + vyDMW**2 + vzDMW**2)
    # total velocity for MW
    vtotSYS = np.sqrt(vxSYS**2 + vySYS**2 + vzSYS**2)
    
    # Arrays for r and v for M31
    rM31 = np.array([xDM31,yDM31,zDM31]).T # transposed 
    vM31 = np.array([vxDM31,vyDM31,vzDM31]).T
    # Arrays for r and v for MW
    rMW = np.array([xDMW,yDMW,zDMW]).T # transposed 
    vMW = np.array([vxDMW,vyDMW,vzDMW]).T
    # Arrays for r and v for SYS
    rSYS = np.array([xSYS,ySYS,zSYS]).T # transposed 
    vSYS = np.array([vxSYS,vySYS,vzSYS]).T
    
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
        
        return pos, vel, R

# =============================================================================
# =======================================================================
# =======================================================================
    # compute the rotated position and velocity vectors
    rSYSn, vSYSn, R  = RotateFrame(rSYS,vSYS)
    
    # Rotate coordinate system
    rM31n = np.dot(R, rM31.T).T
    vM31n = np.dot(R, vM31.T).T
    # Rotate coordinate system
    rMWn = np.dot(R, rMW.T).T
    vMWn = np.dot(R, vMW.T).T
    
    
    #========================================================================
    #======================================================================== 
    #XY plane rotated SYS:    
    # Plot combined system (SYS)
    # fig, ax= plt.subplots(figsize=(12, 10))
    ax = axsS[i//2,i%2]
    ax.hist2d(rSYSn[:,0], rSYSn[:,1], bins=bincnt, range = ((-halfRange, halfRange), (-halfRange, halfRange)),
               norm=LogNorm(), cmap='magma')
    # plt.hist2d(rSYSn[:,0], rSYSn[:,1], bins=bincnt, range = ((-halfRange, halfRange), (-halfRange, halfRange)),
               # norm=LogNorm(), cmap='magma')
    # cbar = plt.colorbar()
    # cbar.set_label("Number of disk particle per bin for M31", fontsize=15)

    density_contour(rSYSn[:,0], rSYSn[:,1], 80, 80, ax=ax, colors=['red','blue','green','yellow'])
    
    # Add axis labels
    ax.set_xlabel('x (kpc) ', fontsize=8)
    ax.set_ylabel('y (kpc) ', fontsize=8)
    
    ax.set_title('System Particle Density Snap {}'.format(snap))
    
    #set axis limits
    # plt.ylim(-40,40)
    # plt.xlim(-40,40)
    
    #adjust tick label font size
    label_size = 8
    matplotlib.rcParams['xtick.labelsize'] = label_size 
    matplotlib.rcParams['ytick.labelsize'] = label_size
    
    # Save to a file
    # plt.savefig('ResearchAssignment3_CombinedXYHist.png')
    #========================================================================
    #========================================================================
    # Get histogram data for SYS
    Hsys, xedges, yedges = np.histogram2d(rSYSn[:,0], rSYSn[:,1], bins=bincnt,
                                          range = ((-halfRange, halfRange), (-halfRange, halfRange)),
                                          density=False)
    # Hsm = np.ma.masked_where(Hsys == 0, Hsys)
    X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
    
    # Get histogram data for M31
    # Use same edges as SYS for M31
    Hm31, xedges, yedges = np.histogram2d(rM31n[:,0], rM31n[:,1], bins=(xedges, yedges),
                                          density=False)
    
    # x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
    # y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))
    
    # Create array of count fractions M31/SYS
    Grat = Hm31/Hsys
    Gmix = 2*(0.5-np.abs(0.5-Hm31/Hsys))
    
    # Plot Grat fraction array
    # fig, ax = plt.subplots(figsize=(12, 10))
    ax = axsR[i//2,i%2]
    gplt = ax.pcolormesh(X, Y, Grat.T, cmap='rainbow')
    # # gplt = ax.pcolormesh(X, Y, Grat.T, cmap='rainbow')
    cbar = figR.colorbar(gplt,ax=ax)
    # cbar = fig.colorbar(gplt,ax=ax)
    cbar.set_label("Fraction of halo particle that are M31", fontsize=6)
    # Add axis labels
    ax.set_xlabel('x (kpc) ', fontsize=8)
    ax.set_ylabel('y (kpc) ', fontsize=8)
    ax.set_title('Mix Ratio - Snap {}'.format(snap))
    plt.show()
    
    # Plot Gmix fraction array
    # fig, ax = plt.subplots(figsize=(12, 10))
    ax = axsG[i//2,i%2]
    gplt = ax.pcolormesh(X, Y, Gmix.T, cmap='rainbow')
    # # gplt = ax.pcolormesh(X, Y, Grat.T, cmap='rainbow')
    cbar = figR.colorbar(gplt,ax=ax)
    # cbar = fig.colorbar(gplt,ax=ax)
    cbar.set_label("Mixing Number", fontsize=6)
    # Add axis labels
    ax.set_xlabel('x (kpc) ', fontsize=8)
    ax.set_ylabel('y (kpc) ', fontsize=8)
    ax.set_title('Mixing Number - Snap {}'.format(snap))
    plt.show()
    #========================================================================
    #========================================================================
    nbins_x,nbins_y = bincnt,bincnt
    x_H = X.reshape((1,nbins_y))*np.ones((nbins_y,1))
    y_H = Y.reshape((nbins_x,1))*np.ones((1,nbins_x))
    r_H = np.sqrt(x_H**2 + y_H**2)
    index = np.where(Gmix.reshape(bincnt**2,1) > 0.8)
    Rcnt = r_H.reshape(bincnt**2,1)[index]
    
    ax = axsH[i//2,i%2]
    # We can set the number of bins with the *bins* keyword argument.
    ax.hist(Rcnt, bins=20)
    ax.set_title('Unnormalized Mixing Ratio - Snap {}'.format(snap))
    
    hist, bin_edges = np.histogram(Rcnt, bins=20)
    rcnt, be = np.histogram(r_H.reshape(bincnt**2,1), bins=bin_edges)
    hist_n = hist/rcnt
    X = 0.5*(bin_edges[1:]+bin_edges[:-1])
    plt.style.use('_mpl-gallery')
    ax = axsN[i//2,i%2]
    # ax.plot(X,hist_n, 'o-')
    # ax.hist(r_H.reshape((bincnt**2,1)), bins=20)
    width = (bin_edges[1:] - bin_edges[:-1])
    ax.bar(X, hist_n, width=width, color='blue', edgecolor="white", linewidth=0.7)
    ax.set(xlim=(bin_edges[0],bin_edges[-1]), xticks=np.linspace(X[0],X[-1], 10),
       ylim=(0,1), yticks = np.linspace(0, 1, 9, endpoint=True))
    ax.text(0.5, 0.9, 'Mixing Number > 0.8', transform=ax.transAxes)
    ax.set_xlabel('radius (kpc) ', fontsize=8)
    ax.set_ylabel('Normalized Mixing Ratio', fontsize=8)
    ax.set_title('Normalized Mixing Ratio - Snap {}'.format(snap))
    
    #========================================================================
    threshNMN = 0.25
    
    delThresh = threshNMN - hist_n
    crossIdx = np.where(np.sign(delThresh[:-1]) == -np.sign(delThresh[1:])) 
    invSlope = (X[1:]-X[:-1])/(hist_n[1:]-hist_n[:-1])
    delX = invSlope * delThresh[:-1]
    xLoc = X[crossIdx] + delX[crossIdx]
    for xplt in xLoc:
        ax.plot([xplt,xplt],[0,1],'g:')
        ax.text(xplt, 0.6, 'rad: {:.0f} kpc'.format(xplt))
    plt.show()
    
figN.savefig('NormMixRatio')    
# figG.savefig('MixingNum')  
# figR.savefig('MixRatio')        
  