

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMassMod import CenterOfMass




def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.
    inputs:
        galaxy = (float) name of the galaxy 
        star = (int) the number of the first snapshot to be read in
        end = (int) the number of the last snapshot to be read in
        n = (int) an interger indicating the integrals over which you will return the COM
        
          
    outputs: 
        fileout = text file that contains the time, COM position and velocity vectors of 
                  a given galaxy in each snapshot 
    """
    
    # compose the filename for output
    fileout = "Orbit_galaxyname.txt"
    fileout = "Orbit_" + galaxy + ".txt"
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    delta = 0.1
    VolDec = 2
    # for M33 that is stripped more, use different values for VolDec
    #if galaxy == "M33":
        #VolDec = 4
        
        
        
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
        
    snap_ids = np.arange(start, end, n, dtype=int)
    snap_length = len(snap_ids)
    if snap_length == 0: 
        raise RuntimeError("Invalid Snapshot ID Inputs")
        
        
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((snap_length,7))
    #print(orbit)
    
    dname = "./"+ galaxy +"/" #prepend snap directory to file name
    
    # a for loop 
    for  i, snap_id in enumerate(snap_ids):# loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000'+str(snap_id)
        ilbl = ilbl[-3:]
        filename = dname +"%s_"%(galaxy)+ilbl+'.txt'
        # Initialize an instance of CenterOfMass class, using disk particles
        
        COM = CenterOfMass(filename, 2)
        

        # Store the COM pos and vel. Remember that now COM_P required VolDec
        position_COM = COM.COM_P(delta, VolDec) #call COM_P
        x_COM = position_COM[0] #store x pos
        y_COM = position_COM[1] #store y pos
        z_COM = position_COM[2] #store z pos
        velocity_COM = COM.COM_V(x_COM, y_COM, z_COM) #call COM_V
        vx_COM = velocity_COM[0] #store vx vel
        vy_COM = velocity_COM[1] #store vy vel
        vz_COM = velocity_COM[2] #store vz vel
        time = (COM.time/1000).value #get time
       
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        orbit[i][0] = time #store ith time
        orbit[i][1] = x_COM.value #store ith x
        orbit[i][2]= y_COM.value #store ith y
        orbit[i][3] = z_COM.value #store ith z
        orbit[i][4] = vx_COM.value #store ith vx
        orbit[i][5] = vy_COM.value #store ith vy
        orbit[i][6] = vz_COM.value #store ith vz
        #print(orbit) this was test
        # note that you can store 
        # a[i] = var1, *tuple(array1)

        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
    return fileout #to use in future equations


#OrbitCOM("MW", 0, 26, 5) #test with small number files

# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first!
if True: #if statement used so that if I wanted to run plot multi times did not have to go through whole ORBITCOM of each galaxy over and over
    MW_Orbit = OrbitCOM("MW", 0, 801, 5) #make MW COM file
    M31_Orbit = OrbitCOM("M31", 0, 801, 5) #make M31 COM file
    M33_Orbit = OrbitCOM("M33", 0, 801, 5) #make M33 COM file
    print('MW_Orbit =',MW_Orbit) #print name of MW file(to use later in code if did not run it)
    print('M31_Orbit =',M31_Orbit) #print name of M31 file(to use later in code if did not run it)
    print('M33_Orbit =',M33_Orbit) #print name of M33 file(to use later in code if did not run it)


#MW_Orbit = Orbit_MW.txt test
#M31_Orbit = Orbit_M31.txt test
#M33_Orbit = Orbit_M33.txt test

# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt


MWdata = np.genfromtxt(MW_Orbit, dtype=None, names=True) #read in MW 
M31data = np.genfromtxt(M31_Orbit,dtype=None, names=True) #read in M31
M33data = np.genfromtxt(M33_Orbit,dtype=None, names=True) #read in M33
    


MW_t = MWdata['t'] #get time data
MW_x = MWdata["x"] #get x pos data
MW_y = MWdata["y"] #get y pos data
MW_z = MWdata["z"] #get z pos data
MW_vx = MWdata["vx"] #get vx vel data
MW_vy = MWdata["vy"] #get vy vel data
MW_vz = MWdata["vz"] #get vz vel data

M31_x = M31data["x"] #get x pos data
M31_y = M31data["y"] #get y pos data
M31_z = M31data["z"] #get z pos data
M31_vx = M31data["vx"] #get vx vel data
M31_vy = M31data["vy"] #get vy vel data
M31_vz = M31data["vz"] #get vz vel data

M33_x = M33data["x"] #get x pos data
M33_y = M33data["y"] #get y pos data
M33_z = M33data["z"] #get z pos  data
M33_vx = M33data["vx"] #get vx vel data
M33_vy = M33data["vy"] #get vy vel data
M33_vz = M33data["vz"] #get vz vel data



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def mag_dif(x1, x2,y1,y2,z1,z2): 
    x_dif = x1-x2
    y_dif = y1-y2
    z_dif = z1-z2
    mag = np.sqrt(x_dif**2+y_dif**2+z_dif**2)
    return mag


# Determine the magnitude of the relative position and velocities 

# of MW and M31

MW_M31_rel_pos = mag_dif(MW_x, M31_x, MW_y, M31_y, MW_z, M31_z)
MW_M31_rel_vel = mag_dif(MW_vx, M31_vx, MW_vy, M31_vy, MW_vz, M31_vz)
#print(MW_M31_rel_pos)
#print(MW_M31_rel_vel)



# of M33 and M31
M33_M31_rel_pos = mag_dif(M33_x, M31_x, M33_y, M31_y, M33_z, M31_z)
M33_M31_rel_vel = mag_dif(M33_vx, M31_vx, M33_vy, M31_vy, M33_vz, M31_vz)
#print(MW_M31_rel_pos)
#print(MW_M31_rel_vel)



fig = plt.figure() #make a 3d pos plot just for fun
ax = fig.add_subplot(projection='3d') #make it a 3d projection
ax.plot(MW_x, MW_y, MW_z, 'r', markersize=2) #plot MW x,y, and z
ax.plot(M31_x, M31_y, M31_z, 'g', markersize=2) #plot M31 x,y, and z


# Plot the Orbit of the galaxies 
#################################
#fig = plt.figure() test
fig,ax = plt.subplots() #make 2D plot
ax.plot( MW_t, MW_M31_rel_pos, 'ro-') #plot pos
ax.grid(visible=True) #make grid
ax.set_xlabel('Time (Gyr)') #label x axis
ax.set_ylabel('Relative Position Magnitude (kpc)') #label y axis
ax.set_title('MW and M31 Relative Position') #name plot
plt.show() #show it

#plot of MW and M31 pos with log y
fig,ax = plt.subplots() #make 2D plot
ax.semilogy( MW_t, MW_M31_rel_pos, 'ro-') #plot pos
ax.grid(visible=True) #make grid
ax.set_xlabel('Time (Gyr)') #label x axis
ax.set_ylabel('Relative Position Magnitude (kpc)') #label y axis
ax.set_title('MW and M31 Relative Position (with log y)') #name plot
plt.show() #show it

#fig = plt.figure()  test
fig,ax = plt.subplots() #make 2D plot
ax.plot( MW_t, MW_M31_rel_pos, 'ro-') #plot pos
ax.grid(visible=True) #make grid
ax.set_xlabel('Time (Gyr)') #label x axis
ax.set_ylabel('Relative Position Magnitude (kpc)') #label y axis
ax.set_title('M33 and M31 Relative Position') #name plot
plt.show() #show it



# Plot the orbital velocities of the galaxies 
#################################
#fig = plt.figure() test
fig,ax = plt.subplots() #make 2D plot
ax.plot( MW_t, MW_M31_rel_vel, 'ro-') #plot vel
ax.grid(visible=True) #make grid
ax.set_xlabel('Time (Gyr)') #lbel x axis
ax.set_ylabel('Relative Velocity Magnitude (kpc)') #label y axis
ax.set_title('MW and M31 Relative Velocity') #name plot
plt.show() #show plot

#fig = plt.figure()
fig,ax = plt.subplots() #make 2D plot
ax.plot( MW_t, MW_M31_rel_vel, 'ro-') #plot vel
ax.grid(visible=True) #make grid
ax.set_xlabel('Time (Gyr)') #label x axis
ax.set_ylabel('Relative Velocity Magnitude(kpc)') #label y axis
ax.set_title('M33 and M31 Relative Velocity') #name plot
plt.show() #show plot

