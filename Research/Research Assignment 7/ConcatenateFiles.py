# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 16:30:54 2025


"""


def concatenate_snap_files(infile1, infile2, outfile):

    with open(infile1, 'r') as f1in, open(infile2, 'r') as f2in:
        
        f1line1 = f1in.readline()
        f2line1 = f2in.readline()
        if f1line1 != f2line1:
            raise RuntimeError("Snap dates don't match")
            
        with open(outfile, 'w') as fout:
            fout.write(f1line1) 
            
            # read and store total number of particles
            f1line2 = f1in.readline()
            label, value = f1line2.split()
            print('file 1 count: ', value)
            total1 = int(value)
            
            f2line2 = f2in.readline()
            label, value = f2line2.split()
            print('file 2 count: ', value)
            total2 = int(value)
            total = total1 +total2
            line2 = f1line2[:6] + str(int(total)).rjust(11) +'\n'
            fout.write(line2) 
            
            f1line = f1in.readline()
            f2line = f2in.readline()
            fout.write(f1line) 
            
            f1line = f1in.readline()
            f2line = f2in.readline()
            fout.write(f1line)    
            
            for line in f1in:
                fout.write(line)  
                
            for line in f2in:
                fout.write(line)  
                
    return

#######################
if __name__ == '__main__' : 
    
    infile1 ='MW_801.txt'
    infile2 ='M31_801.txt'
    outfile = 'CONC_801.txt'
    
    concatenate_snap_files(infile1, infile2, outfile)
    
    