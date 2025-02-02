#!/usr/bin/env python
# coding: utf-8

# In[148]:


import numpy as np
import astropy.units as u

def ParticleInfo(filename, particle_type, particle_number):
    """This function will calculate the distance, velocity, and mass of a particle)
        Inputs : filename : the file which is meant to read in a MW file
                 particle_type: specified in the file, possible particle types are Dark Matter,
                                Disk Stars, or Bulge Stars
                 particle_number: the number of the particle of that particular particle type

        Outputs : Distance_Mag : The magnitude of the distance of the particle from the 
                                 Center of Mass of the Milky Way in units of kpc
                  Velocity_Mag : The magnitude of the velocity of the particle             
     """
    from ReadFile import Read #imports the Read function from ReadFile so we can read and get the array data
    time,particles,data = Read(filename) #unpack the tuple from Read so we can get time,particles, and data sepeartely 
  
    particle_type_tuple = np.where(data['type'] == particle_type) #get a tuple index of all the row numbers of whichever particle type 
                                                                  #is passed into the function
    particle_type_index, = particle_type_tuple  #unpack the tuple index so we can use values in it
    
    x = data['x'][particle_type_index[particle_number]] #gets the x component of the particle's distance from the center of mass  
                                                        #of the Milky Way in kpc 
    y = data['y'][particle_type_index[particle_number]] #gets the y component of the particle's distance from the center of mass  
                                                        #of the Milky Way in kpc 
    z = data['z'][particle_type_index[particle_number]] #gets the z component of the particle's distance from the center of mass  
                                                        #of the Milky Way in kpc 
    vx = data['vx'][particle_type_index[particle_number]] #gets the x component of the particle's velocity  in km/s measured in a Cartesian
                                                          #coordinate system centered on the location of the Milky Way
    vy = data['vy'][particle_type_index[particle_number]] #gets the y component of the particle's velocity  in km/s measured in a Cartesian
                                                          #coordinate system centered on the location of the Milky Way
    vz = data['vz'][particle_type_index[particle_number]] #gets the x component of the particle's velocity  in km/s measured in a Cartesian
                                                          #coordinate system centered on the location of the Milky Way
    mass = data['m'][particle_type_index[particle_number]] #gets the mass of the particle in units of 10^10 Solar Mass

    Distance_Mag = np.around((np.sqrt((x**2 + y**2 + z**2))*u.kpc),3) #uses the 3 components of the distance and calculates the 
                                                                      #magnitude of the particle's distance in kpc
    
    Velocity_Mag = np.around((np.sqrt((vx**2 + vy**2 + vz**2))*u.km/u.s),3) #uses the 3 components of the velocity and calculates the 
                                                                            #magnitude of the particles's velocity in km/s

    Mass = (mass*10**10)*u.M_sun #converts the mass of the particle in solar mass
    
    
    return Distance_Mag, Velocity_Mag, Mass #returns the magnitude of the distance, the magnitue of the velocity, and the mass

#ParticleInfo('MW_000-Copy1.txt',2,99)  used this line to get answers to Part 5 of Homework 2
#Distance_Mag, Velocity_Mag, Mass = ParticleInfo('MW_000-Copy1.txt',2,99) used this line to unpack tuple
#print(np.around(Distance_Mag.to(u.lyr),3)) used this line to convert the distance of the particle from kpc to lightyears


# In[ ]:





# In[ ]:





