# -*- coding: utf-8 -*-
"""

@author: catherine
Code citation: used jupyter notebook by Himansh to get photoutils to work
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from ReadFile import Read
from CenterOfMass import CenterOfMass as COM
#you would need to install the packages below, uncomment the line below
#!pip install photutils
from photutils.isophote import EllipseGeometry
from photutils.aperture import EllipticalAperture
from photutils.isophote import Ellipse

#Loads Snapshot
# file = './MW_000.txt' #t = 0 particle data for the MW
file = "MW_610.txt"
_, _, data = Read(file) #reading particle data

#Extract Dark Matter Halo Particles
dm_idxs = np.where(data['type'] == 1) #extracting the indices corresponding to the DM particles
#extracting positions of the dark matter particles
x = data['x'][dm_idxs]
y = data['y'][dm_idxs]
z = data['z'][dm_idxs]

#Centering the Dark Mtter Halo
myMWHalo = COM(file, 1) #creating the center of mass object
MW_COM = myMWHalo.COM_P(0.1) #center of mass of the dark matter halo using the shrinking sphere
x_COM, y_COM, z_COM = MW_COM[0], MW_COM[1], MW_COM[2] #note: these are astropy quantities

#centering
x1 = x - x_COM.value
y1 = y - y_COM.value
z1 = z - z_COM.value

#Plot the X-Y Projection of the Dark Matter Halo
fig = plt.figure(figsize = (8, 6))
ax0 = fig.add_subplot(1, 1, 1)
_, _, _, plot = plt.hist2d(x1, y1, range = ((-200, 200), (-200, 200)), bins = 100, norm = LogNorm(), cmap = plt.cm.inferno)
fig.colorbar(plot, ax = ax0, label = 'DM particle counts')
ax0.set_xlabel('x [kpc]')
ax0.set_ylabel('y [kpc]')
plt.show()

#Construct a 2-D histogram of the dark matter particles distribution
counts, xedges, yedges = np.histogram2d(x1, y1, bins = 100, range = ((-200, 200), (-200, 200)))
#our bin size is:
x_binsize = xedges[1] - xedges[0]
y_binsize = yedges[1] - yedges[0]

#converting the center into pixel units
x_cen_pixel = -xedges[0]/x_binsize
y_cen_pixel = -yedges[0]/y_binsize

#defining the ellipse geometry
geometry = EllipseGeometry(x0 = x_cen_pixel, y0 = y_cen_pixel, sma = 50/x_binsize, eps = 0, pa = 0)
#sma is the semi-major axis, eps is the ellipticity and pa is the position angle
#constructing the ellipse
aper = EllipticalAperture((geometry.x0, geometry.y0), geometry.sma, geometry.sma * (1 - geometry.eps), 
                           geometry.pa) 

#plot the guess ellipse
#use imshow to plot the distribution of dark matter in pixel coordinates
fig = plt.figure(figsize = (6, 6))
ax0 = fig.add_subplot(1, 1, 1)
cm = plt.cm.inferno
im = plt.imshow(np.log10(counts).T, cmap = cm, origin = "lower") #imshow transposes the matrix, so transposing again to reverse
ax0.set_xlabel('x [pixel]')
ax0.set_ylabel('y [pixel]')
aper.plot(color = 'pink', linewidth = 3)
plt.show()

#Performing the fit
ellipse = Ellipse(counts.T, geometry) #initializing the ellipse object
#you need to transpose the counts matrix to reverse the effect of transposition that Ellipse does internally
isolist = ellipse.fit_image() 

#isolist stores a list of objects - each object corresponds to an elliptical iso-contour that is fitted to the image.

#Plotting the Fitted Ellipses
#plot the fitted elliptical contours for a range of semi-major axis, from 10 kpc to 150 kpc
fig = plt.figure(figsize = (8, 6))
ax0 = fig.add_subplot(1, 1, 1)
cm = plt.cm.inferno
im = plt.imshow(np.log10(counts).T, cmap = cm, origin = "lower")
ax0.set_xlabel('x [pixel]')
ax0.set_ylabel('y [pixel]')

#defining the semi-major axes in pixel coordinates
sma_min = 10/x_binsize
sma_max = 150/x_binsize
separation = 20/x_binsize
sma_array = np.arange(sma_min, sma_max, separation) #each contour will be separated by 20 kpc

for my_sma in sma_array:
    iso = isolist.get_closest(my_sma)
    #the get_closest method obtains the elliptical contour whose semi-major axis is closest to what you want
    x, y = iso.sampled_coordinates()
    ax0.plot(x, y, color = 'pink', linewidth = 3)

plt.show()

#Extracting the Properties of the Fitted Ellipses

# extract the ellipse with a semi-major axis of 100 kpc
sma_kpc = 100 #required semi-major axis in kpc
sma_pixel = 100/x_binsize #required semi-major axis in pixel units

#obtain the elliptical contour whose semi-major axis is closest to 100 kpc
iso_100 = isolist.get_closest(sma_pixel)

#iso_100 is itself an object, and it has various methods which correspond to various properties of this ellipse
eps_100 = iso_100.eps
print("Ellipticity of 100 kpc elliptical contour: ", np.round(eps_100, 2))
pa_100 = (iso_100.pa)*180/np.pi #converting to degrees
print("Position angle of 100 kpc elliptical contour in degrees: ", np.round(pa_100, 2))
#In extracting any property of the ellipse that has units (like semi-minor axis), make sure that you convert back to physical units from pixel units