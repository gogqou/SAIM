#!/usr/local/share/bin/EPD/python -tt

'''
Created on Jul 20, 2012

@author: gou
'''

import sys
import csv
import numpy as np
import os

def readcsv(home_dir, filename):
    #this function is here so I can play with any parameters that need to be changed in case of a csv vs. a txt file
    #07242012
    
    print 'file: ' + home_dir + filename
    os.chdir(home_dir)
    f=open(filename,'r')
    lines = f.readlines()
    #size =[4,4]
    size =[512,512]
    list_from_text= np.zeros(size)
    #enumerate numbers the lines, the 1 specifies that n starts at 1 instead of 0   
    for n, line in enumerate(lines, 1):
        #strips away the \n delimiter at end of each line that makes a new line
        lineagain = line.rstrip()
        #uses ',' as delimiter to separate out the relevant information
        line_split = lineagain.split(',')
        list_from_text [n-1,:] = line_split        
    return list_from_text

def readinput(home_dir, filename):
    print 'file: ' + home_dir + filename
    os.chdir(home_dir)
    f=open(filename,'r')
    lines = f.readlines()
    #size =[4,4]
    size =[512,512]
    
    list_from_text= np.zeros(size)
        
    #enumerate numbers the lines, the 1 specifies that n starts at 1 instead of 0   
    for n, line in enumerate(lines, 1):
        #strips away the \n delimiter at end of each line that makes a new line
        lineagain = line.rstrip()
        
        #uses ',' as delimiter to separate out the relevant information
        line_split = lineagain.split(',')
        list_from_text [n-1,:] = line_split        
 

    return list_from_text
def populate_array (home_dir, filename):
    array = readinput(home_dir,filename)
    #array = readcsv(home_dir,filename)
    result_array = array
    return result_array


def writeoutput(home_dir, filename,outputname):
    array = populate_array (home_dir, filename)
    print array
    print 'output file: ' + home_dir+outputname
    os.chdir(home_dir)
    f = open(outputname, 'w')
    writer = csv.writer(f)
    for i in range(len(array)):
        writer.writerow(array[i,:])
    
    return 1


def readinput_alternative(home_dir, filename,line_num,index):
    
    print 'file: ' + home_dir + filename
    os.chdir(home_dir)
    f=open(filename)
    lines = f.readlines()
    line_num = int(line_num)
    index = int(index)
    size =[512,512]
    #size =[4,4]
    list_from_text= np.zeros(size)
    for n, line in enumerate(lines, 1):
        if n <= line_num:
            #standin = str(n) + '. ' + line
            print 'line number:' + str(n) + '\nwhole line is: ' + line  
            #strips away the \n delimiter at end of each line that makes a new line
            lineagain = line.rstrip()
            #uses ',' as delimiter to separate out the relevant information
            line_split = lineagain.split(',')
            list_from_text [n-1,:] = line_split        
            value = line_split[index]
            print 'position: '+str(index+1)+ '\nvalue: '+value
        else: break
    return [value,list_from_text]






def main():
    if len(sys.argv) != 6:
        print 'Provide something'
        sys.exit(1)
    home_dir = sys.argv[1]
    filename = sys.argv[2]
    line_num = sys.argv[3]
    index = sys.argv[4]
    outputname =sys.argv[5]
    if writeoutput(home_dir,filename,outputname) == 1:
        print 'done reading/writing csvs/txt files'
    else:
        print 'unable to complete command'
        sys.exit(1)
if __name__ == '__main__':
    main()
