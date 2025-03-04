
# Homework 4
# Center of Mass Position and Velocity
# Solutions: G.Besla, R. Li, H. Foote


# remember this is just a template, you don't need to follow every step
# if you have your own method to solve the homework, it is totally fine



# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read




class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot

    def __init__(self, filename, ptype):
        ''' Class to calculate the 6-D phase-space position of a galaxy's center of mass using
        a specified particle type. 
            
            PARAMETERS
            ----------
            filename : `str`
                snapshot file
            ptype : `int; 1, 2, or 3`
                particle type to use for COM calculations
        '''
     
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        self.x = self.data['x'][self.index] #store x pos
        self.y = self.data['y'][self.index] #store y pos
        self.z = self.data['z'][self.index] #store z pos
        self.vx = self.data['vx'][self.index] #store vx velocity
        self.vy = self.data['vy'][self.index] #store vy velocity
        self.vz = self.data['vz'][self.index] #store vz velocity


    def COMdefine(self,a,b,c,m):
        ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
        
        PARAMETERS
        ----------
        a : `float or np.ndarray of floats`
            first vector component
        b : `float or np.ndarray of floats`
            second vector component
        c : `float or np.ndarray of floats`
            third vector component
        m : `float or np.ndarray of floats`
            particle masses
        
        RETURNS
        -------
        a_com : `float`
            first component on the COM vector
        b_com : `float`
            second component on the COM vector
        c_com : `float`
            third component on the COM vector
        '''
        sum_m = sum(m)
        # write your own code to compute the generic COM 
        #using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        a_com = sum(a*m)/sum_m
        # ycomponent Center of mass
        b_com = sum(b*m)/sum_m
        # zcomponent Center of mass
        c_com = sum(c*m)/sum_m
        #print(a_com, b_com, c_com, np.sqrt(a_com**2 + b_com**2 + c_com**2))
        # return the 3 components separately
        return a_com, b_com, c_com
    
    def COM_P(self, delta, volDec):
        '''Method to compute the position of the center of mass of the galaxy 
        using the shrinking-sphere method.

        PARAMETERS
        ----------
        delta : `float, optional`
            error tolerance in kpc. Default is 0.1 kpc
        
        RETURNS
        ----------
        p_COM : `np.ndarray of astropy.Quantity'
            3-D position of the center of mass in kpc
        '''                                                                     

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        x_COM, y_COM, z_COM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        r_COM = np.sqrt(x_COM**2 + y_COM**2 + z_COM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # write your own code below
        x_new = self.x - x_COM
        y_new = self.y - y_COM
        z_new = self.z - z_COM
        r_new = np.sqrt(x_new**2 + y_new**2 + z_new**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/volDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            #print ("CHANGE1 = ", change) 
            #print(r_max)
            # select all particles within the reduced radius (starting from original x,y,z, m)
            # write your own code below (hints, use np.where)
            index2 = np.where(r_new <= r_max) #index of particles in range
            x2 = self.x[index2] #index x values in range
            y2 = self.y[index2] #index y values in range
            z2 = self.z[index2] #index z values in range
            m2 = self.m[index2] #index m values in range

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            x_COM2, y_COM2, z_COM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position
            # write your own code below
            r_COM2 = np.sqrt(x_COM2**2+y_COM2**2+z_COM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(r_COM - r_COM2)
            # uncomment the following line if you want to check this                                                                                               
            #print ("CHANGE2 = ", change)                                                                                     

            # Before loop continues, reset : r_max, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max /= volDec
            # check this.                                                                                              
            #print ("maxR", r_max)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            x_new = self.x - x_COM2  # Bug Fix
            y_new = self.y - y_COM2  # Bug Fix
            z_new = self.z - z_COM2  # Bug Fix
            r_new = np.sqrt(x_new**2+y_new**2+z_new**2)

            # set the center of mass positions to the refined values                                                   
            x_COM = x_COM2
            y_COM = y_COM2
            z_COM = z_COM2
            r_COM = r_COM2

            # create an array (np.array) to store the COM position                                                                                                                                                       
            p_COM = np.array([x_COM, y_COM, z_COM])

        # set the correct units using astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        p_COM = np.round(p_COM,3)*u.kpc
        return p_COM
        
        
    def COM_V(self, x_COM, y_COM, z_COM):
        ''' Method to compute the center of mass velocity based on the center of mass
        position.

        PARAMETERS
        ----------
        x_COM : 'astropy quantity'
            The x component of the center of mass in kpc
        y_COM : 'astropy quantity'
            The y component of the center of mass in kpc
        z_COM : 'astropy quantity'
            The z component of the center of mass in kpc
            
        RETURNS
        -------
        v_COM : `np.ndarray of astropy.Quantity'
            3-D velocity of the center of mass in km/s
        '''
        
        # the max distance from the center that we will use to determine 
        #the center of mass velocity                   
        rv_max = 15.0*u.kpc
        
        x_COM = x_COM/u.kpc #taking out the units for computing purposes
        y_COM = y_COM/u.kpc #taking out the units for computing purposes
        z_COM = z_COM/u.kpc #taking out the units for computing purposes

        # determine the position of all particles relative to the center of mass position (x_COM, y_COM, z_COM)
        # write your own code below
        # Note that x_COM, y_COM, z_COM are astropy quantities and you can only subtract one astropy quantity from another
        # So, when determining the relative positions, assign the appropriate units to self.x
        xV = self.x - x_COM
        yV = self.y - y_COM
        zV = self.z - z_COM
        rV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(rV <= (rv_max/u.kpc))
        
        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vx_new = self.vx[indexV]  # Bug Fix
        vy_new = self.vy[indexV]  # Bug Fix
        vz_new = self.vz[indexV]  # Bug Fix
        m_new =  self.m[indexV]  # Bug Fix
        
        # compute the center of mass velocity using those particles
        # write your own code below
        vx_COM, vy_COM, vz_COM = self.COMdefine(vx_new, vy_new, vz_new, m_new)
        
        # create an array to store the COM velocity
        # write your own code below
        v_COM = np.array([vx_COM, vy_COM, vz_COM])
        v_COM = np.round(v_COM,3)*u.km/u.s
      

        # return the COM vector
        # set the correct units usint astropy
        # round all values  
        return v_COM
                                                                                      
     
    

# ANSWERING QUESTIONS
#######################
if __name__ == '__main__' : 

    # Create a Center of mass object for the MW, M31 and M33
    # below is an example of using the class for MW
    MW_COM = CenterOfMass("C:\\Users\\piguy\\Catherine's HW\\Astr400B\\MW_000.txt", 2)
    M31_COM = CenterOfMass("C:\\Users\\piguy\\Catherine's HW\\Astr400B\\M31_000.txt", 2)
    M33_COM = CenterOfMass("C:\\Users\\piguy\\Catherine's HW\\Astr400B\\M33_000.txt", 2)


    # below gives you an example of calling the class's functions
    # MW:   store the position and velocity COM
    
    #Answer to question 1 of part 6
    MW_COM_p = MW_COM.COM_P(0.1)
    print("COM position of MW:")
    print(MW_COM_p)
    MW_COM_v = MW_COM.COM_V(MW_COM_p[0], MW_COM_p[1], MW_COM_p[2])
    print("COM velocity of MW:")
    print(MW_COM_v)

    # now write your own code to answer questions
    M31_COM_p = M31_COM.COM_P(0.1)
    print("COM position of M31:")
    print(M31_COM_p)
    M31_COM_v = M31_COM.COM_V(M31_COM_p[0], M31_COM_p[1], M31_COM_p[2])
    print("COM velocity of M31:")
    print(M31_COM_v)
    
    M33_COM_p = M33_COM.COM_P(0.1)
    print("COM position of M33:")
    print(M33_COM_p)
    M33_COM_v = M33_COM.COM_V(M33_COM_p[0], M33_COM_p[1], M33_COM_p[2])
    print("COM velocity of M33:")
    print(M33_COM_v)
    
    #Answer to question 2 of part 6
    MW_COM_pos_vector = np.array(MW_COM_p) #turning tuple into array/vector
    M31_COM_pos_vector = np.array(M31_COM_p) #turning tuple into array/vector
    MW_M31_pos_vector = MW_COM_pos_vector-M31_COM_pos_vector #separation vector
    MW_M31_pos_mag = np.round(np.sqrt(sum(MW_M31_pos_vector**2)),3)*u.kpc #mag of vector
    print("Magnitude of current separation of MW and M31:")
    print(MW_M31_pos_mag)
    
    MW_COM_vel_vector = np.array(MW_COM_v) #turning tuple into array/vector
    M31_COM_vel_vector = np.array(M31_COM_v) #turning tuple into array/vector
    MW_M31_vel_vector = MW_COM_vel_vector-M31_COM_vel_vector #separation vector
    MW_M31_vel_mag = np.round(np.sqrt(sum(MW_M31_vel_vector**2)),3)*u.kpc #mag of vector
    print("Magnitude of current separatiion velocity of MW and M31:")
    print(MW_M31_vel_mag)
    
    #Answer to question 3 of part 6 
    M33_COM_pos_vector = np.array(M33_COM_p) #turning tuple into array/vector
    M31_COM_pos_vector = np.array(M31_COM_p) #turning tuple into array/vector
    M33_M31_pos_vector = M33_COM_pos_vector-M31_COM_pos_vector #separation vector
    M33_M31_pos_mag = np.round(np.sqrt(sum(M33_M31_pos_vector**2)),3)*u.kpc #mag of vector
    print("Magnitude of current separation of M33 and M31:")
    print(M33_M31_pos_mag)
    
    M33_COM_vel_vector = np.array(M33_COM_v) #turning tuple into array/vector
    M31_COM_vel_vector = np.array(M31_COM_v) #turning tuple into array/vector
    M33_M31_vel_vector = M33_COM_vel_vector-M31_COM_vel_vector #separation vector
    M33_M31_vel_mag = np.round(np.sqrt(sum(M33_M31_vel_vector**2)),3)*u.kpc #mag of vector
    print("Magnitude of current separatiion velocity of M33 and M31:")
    print(M33_M31_vel_mag)
    
    #Answer to question 4 of part 6
    print("The iterative process is important because with the galaxies merging")
    print("the shape of the galaxies will be distorted and the COM may be skewed")



