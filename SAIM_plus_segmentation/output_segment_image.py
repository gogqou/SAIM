#!/usr/local/share/bin/EPD/python -tt
#!C:\EPD\python -tt
'''
Created in October, 2012

@author: gogqou, au


inputs:
1) directory with segment_adhesions.csvs are (these are the csvs where each pixel is either 0 or has an assigned adhesion number to it
they are generated in the middle of the M. Berginski code and are saved in each individual pictures folder, outside of the raw_data folder
with the actual images
2) this csv's file name (usually adhesions_segment.csv)
3) directory with the processed segment+heights csv where the adhesion numbers for each pixel have been replaced with their calculated height from the 
SAIM image processing 
4) this above csv's file name, usually under the original data files, under segmented_csvs folder, as cell##t##_cH.csv
5) output directory
6) output name

order of operations:

1a)take the UNC program segmented csv, which is a 512x512 grid with either 0s or a number corresponding to adhesion number in each field
1b)pulls out the values and puts it into a matrix

2) reads through each line of the matrix:
3) creates a dummy bitmap image that can then be drawn later, after its pixel values have been assigned
4) for each pixel (entry in the matrix), reads the adhesion number (or zero)
5) if non zero, records its position
6) adds the adhesion number to the value_map
7) adds the position to the position map (these two are a ways to keep real time track of the center of this adhesion
8) draws a box in the new image based on the scaling factor from the old image (if scaling factor is 5, then draws a box of 5x pixels )
9) uses the map of adhesion number --> average position to write the adhesion number over that adhesion
10)draw picture


'''

import sys
import Image, ImageDraw
import random
import csv
import os
import numpy as np
#import surfacePlot

def cycle_through_files(home_dir,filename,heights_dir,folders,naming_scheme,end_dir):
    size = [512, 512]
    row = size[0]
    column = size[1]
    #gets lists of the files from the unc program and corresponding heights_csv, using the get_file_list function
    filenames= get_file_list (home_dir, filename,heights_dir,folders,naming_scheme,end_dir)
    #assign the results from that function to indicate which list is which
    segmentfiles = filenames[0] #unc csvs--512 x512, with pixels with adhesions indicated by the number (#) adhesion it belongs to
    heightsfiles=filenames[1] #heights csvs--512x512, with values at each pixel corresponding with the SAIM calculated heights
    outputfiles = filenames[2]
    outputfiles2 = filenames[3]
    #for each segment file, get the csv data into an array
    for i in range(len(segmentfiles)):
        segment_array = readcsv(segmentfiles[i])
        #for each heights csv, get the data into an array
        heights_array = readcsv(heightsfiles[i])
        #allocate space for the new segmented_heights array
        segmented_heights_array = np.zeros(size)
        adhesion_array = np.zeros([800,4])        
        writeoutput(segment_array, end_dir, 'segment_array_orig_'+str(i+1)+'.csv')
        writeoutput(heights_array, end_dir, 'height_array_orig_'+str(i+1)+'.csv')
        #cycle through the segment array and if the value is not zero, go to the heights csv at the same position, and populate the segmented_heights array with that number
        adh_row=0
        for k in range(row):
            for j in range(column):
                rownum=k
                colnum=j
                
                if segment_array[rownum,colnum]>0 and heights_array[rownum,colnum]>0:
                    segmented_heights_array[rownum,colnum] = heights_array [rownum,colnum]
                    #print 'row '+str(rownum+1) +' column ' + str(colnum+1)
                    
                    #also populate the 4 column array with adhesion number, row number, column number
                    adhesion_array[adh_row,0] = segment_array[rownum,colnum]
                    adhesion_array[adh_row,1] = rownum+1
                    adhesion_array[adh_row,2] = colnum+1
                    adhesion_array[adh_row,3] = segmented_heights_array[rownum,colnum]
                    #increment the row of the adhesion_array that we will be writing to
                    adh_row+=1
                else:
                    continue  
        adhesion_array = adhesion_array[0:adh_row]
        writeoutput(segmented_heights_array,end_dir,outputfiles[i])
        writeoutput(adhesion_array,end_dir,outputfiles2[i])
       
def draw():
    colormap = {0:(0, 0, 0)};

    positions = {}
    val_n = {}

    def readFile(filename):
        f = open(filename)
        lines = f.readlines()
        f.close()
    
        matrix = []
    
        for line in lines:
            nums = line[:-1].split(',')
            lineMatrix = [];
            for n in nums:
                lineMatrix += [int(n)]
    
            matrix += [lineMatrix]
    
        #print matrix
        return matrix
    
    def makeImage(imgName, matrix):
        height = len(matrix)
        width = len(matrix[0])
    
        factor = 5
    
        img = Image.new('RGB', (factor*width, factor*height), 'white')
        draw = ImageDraw.Draw(img)
    
        for x in xrange(height):
            for y in xrange(width):
    
                val = matrix[x][y]
    
                if val in colormap:
                    color = colormap[val]
                else:
                    color = (128, random.randint(0, 256), random.randint(0, 256))
                    colormap[val] = color
    
                if val!=0:
                    pos = (x, y)
                    
                    if val in val_n:
                        val_n[val] += 1
                        oldPos = positions[val]
                        newPos = ( (pos[0]+oldPos[0]), (pos[1]+oldPos[1]) )
                        positions[val] = newPos;
                    else:
                        val_n[val] = 1
                        positions[val] = pos
    
                for i in xrange(factor*x, factor*(x+1)):
                    for j in xrange(factor*y, factor*(y+1)):
                        #print (x, y), val
                        img.putpixel((j, i), color)   
    
    
        for val in val_n.keys():
            n = val_n[val]
            pos = positions[val]
            avg_pos = ( factor*pos[0]/n, factor*pos[1]/n)
            draw.text( (avg_pos[1], avg_pos[0]), str(val))
    
        img.save(imgName, 'BMP')
    
    matrix = readFile(sys.argv[1])
    makeImage(sys.argv[2], matrix)
return 1


def get_file_list(home_dir, unc_csv_name,heights_dir,subfolders,naming_scheme,end_dir):
    subfolders= subfolders.split(',')

    #BECAUSE I'M LAZY: JUST CODING IN THESE NUMBERS
    subfoldernum1 = 1 #cell number
    subfoldernum2 = 33 #time points
    
    file_list = []
    heights_list = []
    output_list =[]
    output_list2 =[]
    #iterates through each subfolder in the main home directory
    for i in range(subfoldernum1):
        for j in range(subfoldernum2):
            num1 = i+1
            num2 = j+1
            #necessary because there are 0s before the single digit numbers
            #be aware that this may change for data sets with even more timepoints (like >100)
            if num2 <10:
                csv_dir = home_dir+'0'+str(num2)+'/'+unc_csv_name 
            else:
                csv_dir = home_dir+str(num2)+'/'+unc_csv_name
            #home directory, plus subfolder, plus the actual name of the file (cell1t1_cH.csv, for example)
            heights_file =heights_dir+subfolders[0]+ str(num1) + '/'+subfolders[1]+str(num2)+'/'+subfolders[0]+str(num1)+subfolders[1]+str(num2)+naming_scheme
            output_file = end_dir+subfolders[0]+str(num1)+subfolders[1]+str(num2)+naming_scheme
            output_file2=end_dir+subfolders[0]+str(num1)+subfolders[1]+str(num2)+'_adh_array.csv'
            file_list = file_list+[csv_dir]
            heights_list = heights_list+[heights_file]
            output_list = output_list + [output_file]
            output_list2=output_list2+[output_file2]
            
    return file_list,heights_list,output_list,output_list2


def readcsv(filename):
    print filename
    f=open(filename,'rb')
    lines = f.readlines()
    size =[512,512]
    csv_data= np.zeros(size)
    #enumerate numbers the lines, the 1 specifies that n starts at 1 instead of 0   
    for n, line in enumerate(lines, 1):
        #strips away the \n delimiter at end of each line that makes a new line
        lineagain = line.rstrip()
        #uses ',' as delimiter to separate out the relevant information
        line_split = lineagain.split(',')
        csv_data [n-1,:] = line_split        
    return csv_data

def writeoutput(array,end_dir,outputname):
    os.chdir(end_dir)
    f = open(outputname, 'wb')
    writer = csv.writer(f)
    for i in range(len(array)):
        writer.writerow(array[i,:])
    return 1


def main():
    if len(sys.argv) > 8:
        print 'Provide 1)directory,  \n 2)filename (of the csv file) \n 3) subfolders in the SAIM directory \n 4)directory where SAIM csvs are kept \n 5)namings scheme of SAIM csvs \n 6)final directory'
        sys.exit(1)
    home_dir = sys.argv[1]
    filename = sys.argv[2]
    folders = sys.argv[3]
    heights_dir = sys.argv[4]
    naming_scheme = sys.argv[5]
    #end_dir = "C:\Users\gogqou\Documents\Research\Weaver lab\coding\Data\3ExampleCells\segmented_csvs_2"
    end_dir = sys.argv[6]
    if cycle_through_files(home_dir,filename,heights_dir,folders,naming_scheme,end_dir) == 1:
        print 'done segmenting and writing csvs'
    else:
        print 'unable to complete command'
        sys.exit(1)
if __name__ == '__main__':
    main()