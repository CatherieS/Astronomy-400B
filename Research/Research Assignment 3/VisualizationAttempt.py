# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:27:59 2025

@author: piguy
"""

import numpy as np

import matplotlib.pyplot as plt

from ReadFile import Read

colorChar = ['m.', 'b.','g.']
#colorChar = ['r.', 'b.']
starlim = [1500,25000,100000]

rng = np.random.default_rng()
fig = plt.figure()

ax = fig.add_subplot(projection='3d')

filenames = ['MW_000.txt','M31_000.txt']
time, total, data1 = Read(filenames[0])
time, total, data2 = Read(filenames[1])
data = np.hstack((data1,data2))
for idx in range(3):
    ptype = idx + 1.
    #ptype = idx
    index = np.where(data['type'] == ptype)
   
    p = data[['m','x','y','z']][index]
    psize = p.size
    shflidx = np.arange(psize)
    rng.shuffle(shflidx)
    x = p['x'][shflidx[0:starlim[idx]]] #store y pos
    y = p['y'][shflidx[0:starlim[idx]]] #store y pos
    z = p['z'][shflidx[0:starlim[idx]]] #store z pos
    ax.plot(x, y, z, colorChar[idx], markersize=1)
    # ax.plot(x, y, z, colorChar[idx], markersize=1)
    #ax.plot(y, z, colorChar[idx], markersize=1)    
   



ax.view_init(35,15, 0)
ax.set_xlim([-2000, 2000])
ax.set_ylim([-2000, 2000])
ax.set_zlim([-2000, 2000])
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()