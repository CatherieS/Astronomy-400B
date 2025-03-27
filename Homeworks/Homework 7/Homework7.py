
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component


# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """Function that takes in file that stores orbit and initializes quantities that will be used in future functions
        INPUTS
        ------
        filename: 'float'
            filename for the file in which you will store the integrated orbit
        """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = (const.G.to(u.kpc**3/u.Msun/u.Gyr**2)).value
        
        ### **** store the output file name
        self.filename = filename #"M33PredictedOrbit"
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        # Determine the relative position between M33 and M31 
        self.M31_COM = CenterOfMass("M31_000.txt", 2)
        self.M33_COM = CenterOfMass("M33_000.txt", 2) 
        self.M33_COM_p = self.M33_COM.COM_P(0.1)
        self.M33_COM_v = self.M33_COM.COM_V(self.M33_COM_p[0],self.M33_COM_p[1],self.M33_COM_p[2])   
        self.M31_COM_p = self.M31_COM.COM_P(0.1)
        self.M31_COM_v = self.M31_COM.COM_V(self.M31_COM_p[0],self.M31_COM_p[1],self.M31_COM_p[2])                                                             
        #M33_M31 = np.sqrt((M33_COM_p[0]-M31_COM_p[0])**2 + (M33_COM_p[1]-M31_COM_p[1])**2 + (M33_COM_p[2]-M31_COM_p[2])**2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        #self.r = M33_COM_p[0]-M31_COM_p[0], M33_COM_p[1]-M31_COM_p[1], M33_COM_p[2]-M31_COM_p[2]
        self.r = self.M33_COM_p.value - self.M31_COM_p.value
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        #vM33_M31 = np.sqrt((M33_COM_v[0]-M31_COM_v[0])**2 + (M33_COM_v[1]-M31_COM_v[1])**2 + (M33_COM_v[2]-M31_COM_v[2])**2)
        #self.v = M33_COM_v[0]-M31_COM_v[0], M33_COM_v[1]-M31_COM_v[1], M33_COM_v[2]-M31_COM_v[2]
        self.v = self.M33_COM_v.value - self.M31_COM_v.value
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 

        # **** store the position VECTOR of the M31 COM (.value to get rid of units)

        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = set scale length(no units)
        self.rdisk = 5
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = 0.12*1e12 
        ### bulge
        # ****  self.rbulge = set scale length (no units)
        self.rbulge = 1
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = 0.019*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = 1.921*1e12
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an input the position VECTOR 
        """Function returns the acceleration vector from a Hernquist potential for either the bulge or the halo
        
        INPUTS
        ------
        M: 'int'
            M is the total halo or bulge mass
        r_a: 'int'
            the corresponding scale length
        r: 'array'
            the magnitude of the relative position vector
        OUTPUTS
        -------
        Hern: 'array'
            acceleration vector from a Hernquist potential
        """
        
        ### **** Store the magnitude of the position vector
        #rmag = np.sqrt((self.M33_COM_p[0]-self.M31_COM_p[0])**2 + (self.M33_COM_p[1]-self.M31_COM_p[1])**2 + (self.M33_COM_p[2]-self.M31_COM_p[2])**2)
        rmag = np.sqrt(np.sum(r**2))
        ### *** Store the Acceleration
        Hern =  -self.G*M/(rmag*(r_a+rmag)**2)*r #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        #Hern = -self.G*M/(rmag*(r_a + rmag)**2) * r
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """"Function returns the acceleration vector from a Miyamoto-Nagai profile 
        
        INPUTS
        ------
        M: 'int'
            M is the total mass of the disk
        r_d: 'int'
            the corresponding scale length
        r: 'array'
            the magnitude of the relative position vector
        OUTPUTS
        -------
        MiyaNag: 'array'
            acceleration vector from a Miyamoto-Nagai profile
        """

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        z_d = self.rdisk/5.0 #define z_d
        B = r_d + np.sqrt(r[2]**2 + z_d**2) #define B
        R = np.sqrt(r[0]**2+r[1]**2) #define R
        ZSTUFF = B/np.sqrt(r[2]**2 + z_d**2) #define ZSTUFF
        MiyaNag = -self.G*M/(R**2 + B**2)**1.5*r*np.array([1,1,ZSTUFF])  #get acceleration vector 
        
       
        return MiyaNag
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """function takes as input the 3D position vector (x,y,z) and returns a 3D vector of the total acceleration
        
        INPUTS
        ------
        M: 'int'
            M is the total halo, bulge, or disk mass
        r_d:'int'
            the corresponding scale length of the disk
        r_a: 'int'
            the corresponding scale length of either the bugle or halo
        r: 'array'
            the magnitude of the relative position vector
        OUTPUTS
        -------
        accelsum: 'array'
            3D vector of the total of the sums of all acceleration vectors from each galaxy component
        """

        ### Call the previous functions for the halo, bulge and disk
        bulgeaccel = self.HernquistAccel(self.Mbulge, self.rbulge, self.r)
        haloaccel = self.HernquistAccel(self.Mhalo, self.rhalo, self.r)
        diskaccel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, self.r)
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        accelsum = bulgeaccel + haloaccel + diskaccel
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return accelsum
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """
        function that integrates the equation of motion using the LeapFrog method
        INPUTS
        ------
        dt: 'int'
            time interval for integration
        r: 'array'
            a starting position vector r for the M33 COM position relative to the M31
        v: 'array'
            a starting velocity vector v for the M33 relative to M31
        OUTPUTS
        -------
        rv_new: 'array'
            the new position and velocity vectors 
        """
        
        # predict the position at the next half timestep
        rhalf = r + v*dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf)*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew*dt/2
        rv_new = rnew, vnew
        
        return rv_new# **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """function that loops LeapFrog integrator and saves the orbit as a file
        INPUTS
        ------
        t0: 'int'
            starting time of orbit
        dt: 'int'
            time interval
        tmax: 'int'
            final time of orbit
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros(((int(tmax/dt)+2), 7))
        
        # initialize the first row of the orbit
        #orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        #orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        orbit[0] = t0, self.r[0], self.r[1], self.r[2], self.v[0], self.v[1], self.v[2]
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (i < int(tmax/dt)+2):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t = t +  dt
            # **** store the new time in the first column of the ith row
            orbit[i][0] = t
            
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            self.r, self.v = self.LeapFrog(dt, self.r, self.v)
         
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            orbit[i][1:4] = self.r
            orbit[i][4:7] = self.v
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            print(i, orbit[i])
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i = i + 1
            
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function
"""Plotting for Analysis Part 5"""
if __name__ == '__main__' : 

    IntegratedOrbit = M33AnalyticOrbit("M33PredictedOrbit") #cll the class
    PredictedOrbit = IntegratedOrbit.OrbitIntegration(0, 0.5, 10) #call function to make orbit file
    M33_prediction = np.genfromtxt('M33PredictedOrbit',dtype=None,names=True) #read in file
    M33_pre_pos = np.array([M33_prediction['x'], M33_prediction['y'], M33_prediction['z']]) #put position data into array
    M33_pre_pos_mag = np.sqrt(M33_pre_pos[0]**2 +M33_pre_pos[1]**2 +M33_pre_pos[2]**2) #get magnitude of the position vector
    M33_pre_vel = np.array([M33_prediction['vx'], M33_prediction['vy'], M33_prediction['vz']]) #put position data into array
    M33_pre_vel_mag = np.sqrt(M33_pre_vel[0]**2 + M33_pre_vel[1]**2 +M33_pre_vel[0]**2) #get magnitude of the velocity vector

    # Read in the data files for the orbits of each galaxy that you just created
    # headers:  t, x, y, z, vx, vy, vz
    # using np.genfromtxt
    M31_orbit = np.genfromtxt('Orbit_M31.txt',dtype=None,names=True) 
    M33_orbit = np.genfromtxt('Orbit_M33.txt',dtype=None,names=True) 
    
    # function to compute the magnitude of the difference between two vectors 
    # You can use this function to return both the relative position and relative velocity for two 
    # galaxies over the entire orbit  
    def relative_mag(a, b): 
        """
        Function that computes the magnitude of the difference between two vectors.
        Inputs with shape (3, n) will return n outputs
    
        PARAMETERS
        ----------
        a : `np.ndarray'
            first vector
        b : 'np.ndarray'
            second vector
    
        RETURNS
        -------
        mag : `float or np.ndarray`
            |a-b|
        """
        
        # compute the difference vector
        x = a[0] - b[0] 
        y = a[1] - b[1]
        z = a[2] - b[2]
    
        # return its magnitude
        return np.sqrt(x**2 + y**2 + z**2)
    
    # position vectors of each galaxy with respect to 0,0,0
    M31_pos = np.array([M31_orbit['x'], M31_orbit['y'], M31_orbit['z']])
    M33_pos = np.array([M33_orbit['x'], M33_orbit['y'], M33_orbit['z']])

    # velocity vectors of each galaxy with respect to 0,0,0
    M31_vel = np.array([M31_orbit['vx'], M31_orbit['vy'], M31_orbit['vz']])
    M33_vel = np.array([M33_orbit['vx'], M33_orbit['vy'], M33_orbit['vz']])
    
    # Determine the magnitude of the relative position and velocities of M33 and M31
    M33_M31_relpos = relative_mag(M31_pos, M33_pos)
    M33_M31_relvel = relative_mag(M31_vel, M33_vel)
    
    # set up plots
    fig, axes = plt.subplots(2, 1, figsize=(5, 8), sharex=True)
    
    # set fontsize 
    plt.rcParams['font.size'] = 15
    
    # Plot the Orbit of the galaxies 
    #################################
    ax = axes[0]
    # Plot the separtion of M31 and MW
    ax.plot(M33_prediction['t'],M33_pre_pos_mag.T, c='b', lw=3, label='Predicted M31')
    # Plot the separtion of M33 and M31
    ax.plot(M31_orbit['t'], M33_M31_relpos, c='r', lw=3, ls="-.", label='M33 and M31')
    ax.set(ylabel='COM Separation [kpc]')
    ax.legend()
    
    # Set the title
    ax.set_title("Orbits Comparison", fontsize=15)
    
    # Plot the orbital velocities of the galaxies 
    #################################
    ax = axes[1]
    ax.plot(M33_prediction['t'],M33_pre_vel_mag , c='b',  lw=3, label='Predicted Velocity')
    ax.plot(M31_orbit['t'], M33_M31_relvel, c='r', lw=3, ls="-.", label='M33 and M31')
    ax.set(ylabel='Relative COM Velocity [km/s]', xlabel='Time [Gyr]')
    #adjust tick label font size
    #label_size = 15
    #matplotlib.rcParams['xtick.labelsize'] = label_size 
    #matplotlib.rcParams['ytick.labelsize'] = label_size


    plt.tight_layout()
    plt.savefig('Homework7_orbits.png', bbox_inches='tight')
    
    """
    #########################
    # # Answering Questions
    #########################
    
    1. How do the plots compare?
    The predicted position plot appears to show the distance from M31 increasing almost linearlly where as the 
    orbit from HW 6 stays pretty almost oscillates and decreases over time
    
    2.What missing physics could make the difference?
    The effect of MW as well as the change in direction
    
    The MW is missing in these calculations. How might you include its effects?
    You ccould include its effects by counting its impact into the M31 COM information
    """

