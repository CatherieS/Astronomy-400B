

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
from CenterOfMassFixMod import CenterOfMass




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
    if galaxy == "M33":
        VolDec = 4
        
        
        
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
        
    snap_ids = np.arange(start, end, n, dtype=int)
    snap_length = len(snap_ids)
    if snap_length <0: 
        print("Input ineligible")
        pass
        
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((snap_length,7))
    #print(orbit)
    
    #dname = ".\\"+ galaxy +"\\"
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
        position_COM = COM.COM_P(delta, VolDec)
        x_COM = position_COM[0]
        y_COM = position_COM[1]
        z_COM = position_COM[2]
        velocity_COM = COM.COM_V(x_COM, y_COM, z_COM)
        vx_COM = velocity_COM[0]
        vy_COM = velocity_COM[1]
        vz_COM = velocity_COM[2]
        time = (COM.time/1000).value
       
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        orbit[i][0] = time
        orbit[i][1] = x_COM.value
        orbit[i][2]= y_COM.value
        orbit[i][3] = z_COM.value
        orbit[i][4] = vx_COM.value
        orbit[i][5] = vy_COM.value
        orbit[i][6] = vz_COM.value
        #print(orbit)
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
    return fileout


#OrbitCOM("MW", 0, 26, 5) #test with small number files

# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first!
 
MW_Orbit = OrbitCOM("MW", 0, 801, 5)
M31_Orbit = OrbitCOM("M31", 0, 801, 5)
M33_Orbit = OrbitCOM("M33", 0, 801, 5)


# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt


MWdata = np.genfromtxt(MW_Orbit, dtype=None, names=True)
M31data = np.genfromtxt(M31_Orbit,dtype=None, names=True)
M33data = np.genfromtxt(M33_Orbit,dtype=None, names=True)

MW_t = MWdata["t"]  # ***************************** Added
MW_x = MWdata["x"]
MW_y = MWdata["y"]
MW_z = MWdata["z"]
MW_vx = MWdata["vx"]
MW_vy = MWdata["vy"]
MW_vz = MWdata["vz"]
#
MW_speed = np.sqrt(MW_vx**2 + MW_vy**2 +MW_vz**2)
#
M31_x = M31data["x"]
M31_y = M31data["y"]
M31_z = M31data["z"]
M31_vx = M31data["vx"]
M31_vy = M31data["vy"]
M31_vz = M31data["vz"]
#
M31_speed = np.sqrt(M31_vx**2 + M31_vy**2 +M31_vz**2)
#
M33_x = M33data["x"]
M33_y = M33data["y"]
M33_z = M33data["z"]
M33_vx = M33data["vx"]
M33_vy = M33data["vy"]
M33_vz = M33data["vz"]



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  
def mag_dif(x1, x2,y1,y2,z1,z2): #fix so that it can do it with an array
    x_dif = x1-x2
    y_dif = y1-y2
    z_dif = z1-z2
    mag = np.sqrt(x_dif**2+y_dif**2+z_dif**2)
    return mag


# Determine the magnitude of the relative position and velocities 

# of MW and M31

MW_M31_rel_pos = mag_dif(MW_x, M31_x, MW_y, M31_y, MW_z, M31_z)
MW_M31_rel_vel = mag_dif(MW_vx, M31_vx, MW_vy, M31_vy, MW_vz, M31_vz)
print(MW_M31_rel_pos)
print(MW_M31_rel_vel)



# of M33 and M31
M33_M31_rel_pos = mag_dif(M33_x, M31_x, M33_y, M31_y, M33_z, M31_z)
M33_M31_rel_vel = mag_dif(M33_vx, M31_vx, M33_vy, M31_vy, M33_vz, M31_vz)
print(MW_M31_rel_pos)
print(MW_M31_rel_vel)



fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(MW_x, MW_y, MW_z, 'r', markersize=2)
ax.plot(M31_x, M31_y, M31_z, 'g', markersize=2)
#
#
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,)
ax1.plot( MW_t, MW_x, 'ro-')
ax2.plot( MW_t, MW_y, 'go-')
ax3.plot( MW_t, MW_z, 'bo-')
ax1.grid(visible=True)
ax2.grid(visible=True)
ax3.grid(visible=True)
#
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,)
ax1.plot( MW_t, M31_x, 'ro-')
ax2.plot( MW_t, M31_y, 'go-')
ax3.plot( MW_t, M31_z, 'bo-')
ax1.grid(visible=True)
ax2.grid(visible=True)
ax3.grid(visible=True)
#
#
fig, (ax1, ax2) = plt.subplots(2, 1,)
ax1.plot( MW_t, MW_speed, 'ro-')
ax2.plot( MW_t, MW_speed, 'go-')
ax1.grid(visible=True)
ax2.grid(visible=True)
#
#
# Plot the Orbit of the galaxies 
#################################
fig, ax = plt.subplots()
ax.plot( MW_t, MW_M31_rel_pos, 'ro-')  # Plot some data on the Axes.
ax.grid(visible=True)
ax.set_xlabel('Time')
ax.set_ylabel('relative Position Magnitude')
ax.set_title('MW - M31 Relative Position')
ax.legend()
plt.show() 

#
fig, ax = plt.subplots()
ax.plot( MW_t, M33_M31_rel_pos, 'ro-')  # Plot some data on the Axes.
ax.grid(visible=True)
ax.set_xlabel('Time')
ax.set_ylabel('relative Position Magnitude')
ax.set_title('M33 - M31 Relative Position')
ax.legend()
plt.show() 
# Plot the orbital velocities of the galaxies 
#################################
fig, ax = plt.subplots()
ax.plot( MW_t, MW_M31_rel_vel, 'ro-')  # Plot some data on the Axes.
ax.grid(visible=True)
ax.set_xlabel('Time')
ax.set_ylabel('relative Velocity Magnitude')
ax.set_title('MW - M31 Relative Velocity')
ax.legend()
plt.show() 

#
fig, ax = plt.subplots()
ax.plot( MW_t, M33_M31_rel_vel, 'ro-')  # Plot some data on the Axes.
ax.grid(visible=True)
ax.set_xlabel('Time')
ax.set_ylabel('relative Velocity Magnitude')
ax.set_title('M33 - M31 Relative Velocity')
ax.legend()
plt.show() 
