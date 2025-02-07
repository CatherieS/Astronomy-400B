import numpy as np
import astropy.units as u

def ComponentMass(filename, particle_type):
""" 
This function will calculate the mass of a galaxy component
Inputs: (float) filename : the name of the file that is contains the galaxy particles
        (integer) particle_type : the particle type that makes up the component we want the mass of
Outputs:(integer) component_mass : the mass of all the particles that make up that galaxy component
"""
    from ReadFile import Read #imports Read so we can read file in
    time,particles,data = Read(filename) #unpack the tuple from Read so we can get time,particles, and data sepeartely 
  
    particle_type_tuple = np.where(data['type'] == particle_type) #get a tuple index of all the row numbers of whichever particle type 
                                                                  #is passed into the function
                                                                  
    particle_type_index, = particle_type_tuple #unpack the tuple so we can use its values
    #print(particle_type_index) this was test
    
    
    mass_array = data['m'][particle_type_index] #create an index of mass of each of the particle type we desire
    
    #print(mass_array) this was test
    
    
    sum_of_mass = sum((mass_array)) #sum up the masses we indexed
    #print(sum_of_mass) this was test
    
    
    component_mass = np.round((sum_of_mass/100),3) #convert to correct units and round
    
    
    #print(component_mass) this was a test
    
    
    return component_mass #return our output


#ComponentMass('MW_000.txt', 2) this was test


name_list = ("MW_000.txt", 'M31_000.txt', 'M33_000.txt') #creates list of the files we want to read

#local_group_components = [[0.]*5]*3
local_group_components = np.zeros((3,5),dtype=float) #create an empty array of we will fill in with component masses for table


#this for loop goes through the three files and each particle type to fill in array for table
for galaxy_index in range(len(name_list)):#runs through list of files
        for type_number in range(1,4):  #runs through possible particle types
            typer_index = type_number-1 #accounts for the fact index starts at zero
            fn = name_list[galaxy_index] #value for running through loop
            element_sum = ComponentMass(fn,type_number) #calls on function to calculate component mass
            local_group_components[galaxy_index][typer_index] = element_sum #adds mass to list
        list_galazy_mass = local_group_components[galaxy_index][0:3] # makes a list of all the components of a galaxy
        total_galaxy_mass = sum(list_galazy_mass) #sums list to get a total mass of the galaxy
        local_group_components[galaxy_index][3] = total_galaxy_mass #puts total mass into array
        baryon_fraction = (total_galaxy_mass-local_group_components[galaxy_index][0])/total_galaxy_mass #calaculates baryon fraction of galaxy
        local_group_components[galaxy_index][4] = baryon_fraction #puts the baryon fraction of the galaxy into array
        
local_group_total_mass = sum(local_group_components[:,3]) #sums the total masses of each galaxy to get mass of Local Group   
local_group_dark_matter_sum = sum(local_group_components[:,0]) #sums the dark matter of each galaxy to get total local group dark matter
local_group_baryon_fraction = (local_group_total_mass - local_group_dark_matter_sum)/local_group_total_mass  #calculates baryon fraction of local group

#print(local_group_components) for testing purposes
#print(local_group_total_mass) for testing purposes
#print(local_group_baryon_fraction) for testing purpose


with open('Homework_3_Galaxy_Table.txt','w', encoding="utf-8") as f: #this with creates a file that will store information for table to be used in Latex
    f.write('This is the output table\n\n')
    f.write(str(local_group_components))
    f.write('\n\n This is the Local Group total mass\n\n')
    f.write(str(local_group_total_mass))
    f.write("\n\n This is the Local Group baryon fraction\n\n")
    f.write(str(local_group_baryon_fraction))
print(f.closed)


            


