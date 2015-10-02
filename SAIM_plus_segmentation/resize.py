#!/usr/local/share/bin/EPD/python -tt

'''
Created on Jul 24, 2012

@author: gou


takes input, 1D heights csvs and repositions fields to create a corresponding 512x512 csv (going across then down)

'''

import sys
import csv
import numpy as np
import os
import math
import re

def try_int(s):
    try: return int(s)
    
    except: return s
    
def natural_key(s):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', s)]


def readcsv(home_dir, filename):
    print 'file: ' + home_dir + filename
    os.chdir(home_dir)
    f=open(filename,'r')
    
    #grabs lines
    lines = f.readlines()
    #size =[4,4]
    size =[512,512]
    
    #allocates space
    list_from_text= np.zeros(size[1]*size[0])
    
    for line in lines:
        lineagain = line.rstrip()
            #uses ',' as delimiter to separate out the relevant information
        line_split = lineagain.split(',')
        list_from_text = line_split     
        
    return list_from_text

def populate_array (home_dir, filename):
    
    orig_array = readcsv(home_dir,filename)
    length=len(orig_array)
    size = [math.sqrt(length),math.sqrt(length)]
    
    #allocates space for new 512x512 array
    array = np.zeros(size)
    
    for i in range(length):
        
        #if the number is nonzero, then we go to that index in the array, determined by how the index divides 512, and put it there
        #otherwise, it's just zero, which we've already assigned as the default value
        if float(orig_array[i])>1:
            row_num = (i)/512 
            col_num = math.fmod(i,512)
#            print 'original index:{0} repositioned to row {1} column {2}'.format(i,row_num+1,col_num+1)
            array[row_num,col_num] = orig_array[i]
            
        else:
            continue
    result_array = array
    return result_array


def writeoutput(home_dir, output_dir,filename,outputname):
    #get the 512x512 array
    array = populate_array (home_dir, filename)
#    print array
    print 'output file: ' + output_dir+outputname
    if os.path.exists(output_dir):
        print 'already exists'
    else: 
        os.makedirs(output_dir)
        print 'made new directory'
    #goes and writes the new csv
    os.chdir(output_dir)
    f = open(outputname, 'w')
    #iterate through the rows of the array and writes each row of the csv
    writer = csv.writer(f)
    for i in range(len(array)):
        writer.writerow(array[i,:])
    
    return 1

def main():
    if len(sys.argv) != 3:
        print 'Provide right inputs'
        sys.exit(1)
    home_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    #if there are NO folders in the home directory, use the below body of code
    
    #'''
    filenames = os.listdir(home_dir)
    filenames.sort(key=natural_key)
    count = 0
    for filename in filenames:
        if '.csv' in filename:
            count +=1
            print natural_key(filename)
            output_name = 'cell' + str(count)+'_c_H.csv'
            writeoutput(home_dir,output_dir, filename,filename)
        else: continue
    
    #'''
    
    
    #if there are folders in the home direcotory, use the below body of code
    # to do so, uncomment the " ''' " 
    
    
    '''
    
    folders = os.listdir(home_dir)
    for folder in folders:
        count = 0
        
        new_home_dir = home_dir+'/'+folder+'/'
        new_output_dir = output_dir + '/'+folder+'/'
        
        filenames = os.listdir(new_home_dir)
        os.chdir(new_home_dir)       
        
        
        #sorts the filenames by cell number (does natural sort)
        #sorts by the last number in the file name
        # for instance, if AAAA##BBBB## sort by the ## after BBB
        filenames.sort(key=natural_key)
        
        
        #iterates through each csv in the folder
        #therefore only have files in here that are in this format
        #if it's a csv, it will get resized like this
        #and renamed cell+#
        for filename in filenames:
            if '.csv' in filename:
                
                #keeps track of how many actually cells are in the data
                #this matches the naming scheme used in the UNC / segmentation code
                count +=1
                print filename
                output_name = 'cell' + str(count)+'_c_H.csv'
                writeoutput(new_home_dir,new_output_dir, filename,output_name)
            else: continue
    
    
    '''
    
    
        
    #if writeoutput(home_dir,output_dir) == 1:
        #print 'done reading/writing csvs/txt files'
    #else:
        #print 'unable to complete command'
        #sys.exit(1)
if __name__ == '__main__':
    main()

