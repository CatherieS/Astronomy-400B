#!/usr/bin/env python
# coding: utf-8

# In[77]:


import numpy as np
import astropy.units as u


def Read(filename):
    """This function reads in a file, and gets the time of the data in units of Myr, the 
    number of particles in the file, and the rest of the data in an array
    
        Inputs: filename : the file which is meant to read in a MW file

        Outputs: time : the time of the data 
                 particles : the total number of particles in the file
                 data : the array of the rest of the data in the file after getting the time and particles    
    """

    file = open(filename, 'r') #Opens the file that we read in as filename

    line1 = file.readline() #Reads the first line of the file and makes it into a string
    label, value = line1.split() #Splits the string containing first line of file into a label and a value (time, actual time in a number)
    time = float(value)*u.Myr #converts the value in numbers from the file into a floating point value and puts it in units of Myr 
                              #(gives us the time of the file)

    
    line2 = file.readline() #Reads second line of file and makes it into a string
    label, value = line2.split() #Splits the string containing second line of file into a label and a value (time, actual time in a number)
    particles = float(value) #converts the value in numbers from the file into a floating point value(gives us total particles in file)

    file.close() #Closes file

    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3) #Creates an array to hold the rest of the data from the file
    #print(time)
    #print(particles)
    #print(data)
    

    return time,particles,data #returns the desired outputs as a tuple

time,particles,data = Read('MW_000-Copy1.txt')  #uncouples the tuple
#print(data[13]) test code

  


# In[ ]:





# In[ ]:




