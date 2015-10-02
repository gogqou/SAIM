#!/usr/local/share/bin/EPD/python -tt
#!C:\EPD\python -tt
'''
Created on Jul 20, 2012

@author: gou

inputs: 
    1) results directory for UNC program
    2) UNC csv filename
    3) subfolder names in the SAIM directory, and also how they'll be parsed
    4) directory where SAIM heights csvs are
    5) naming scheme of SAIM csvs; for instance: cell###t###_c_H.csv
    6) target directory for where to put new, segmented heights data

libraries needed: numpy, PIL   

order of operations:

take the UNC program csv, which is a 512x512 grid with either 0s or a number corresponding to adhesion number in each field

go through and for only the fields with a non-zero value, grab the corresponding value (based on location in the 512 by 512 grid)
at the same time, also assign the adhesion number to that particular pixel--that is, add that position, and the height, to the "segment" defining that particular adhesion
--so for instance, let's say we write a four column vector: 1) adhesion number 2) row number 3) column number 4) SAIM height
--if we have an entire list of this, we can then look at average height of adhesion, etc etc



KEEP
PREVIOUS ARGUMENTS LIST 10-19-2012

"/media/Windows Storage/Linux/FAAS_1 (copy).5/results/rub/sample_data/individual_pictures/"
adhesions_segment.csv
cell,t
"/media/Windows Storage/Linux/Data/3ExampleCells/Green-vinculin/csvs/"
_cH.csv
"/media/Windows Storage/Linux/Data/3ExampleCells/Green-vinculin/segmented_csvs/"

KEEP





once we have these two sets of data: 1 512x512 rewritten heights heatmap with just the adhesions that 
passed the UNC program's parameters and the four column vector

--we write the first to a txt file/csv
--we make a heatmap using that csv
--and we write the four column vector to another txt file/csv that contains: adhesion number, x,y, and height for each relevant pixel


we then populate a list with only unique adhesion numbers
sort it by ascending adhesion number
iterate through these adhesion numbers to find the average height for each adhesion
and output a histogram

'''

import sys
import csv
import os
import numpy as np
import matplotlib.pyplot as plt
import pylab
import matplotlib.path as path
import matplotlib.patches as patches

def cycle_through_files(home_dir,filename,heights_dir,folders,naming_scheme,end_dir,minAdh_size=5):
    size = [512, 512]
    row = size[0]
    column = size[1]
    
    
    #gets lists of the files from the unc program and corresponding heights_csv, using the get_file_list function
    filenames= get_file_list (home_dir, filename,heights_dir,folders,naming_scheme,end_dir)
    
    
    #assign the results from that function to indicate which list is which
        
    #unc csvs--512 x512, with pixels with adhesions indicated by the number (#) adhesion it belongs to
    segmentfiles = filenames[0] 
    
    #heights csvs--512x512, with values at each pixel corresponding with the SAIM calculated heights
    heightsfiles=filenames[1] 
    
    outputfiles = filenames[2]
    outputfiles2 = filenames[3]
    
    #cell #
    #subfoldernum1=filenames[4] 
    
    #time points #
    #subfoldernum2=filenames[5] 
    
    
    #for each segment file, get the csv data into an array
    for i in range(len(segmentfiles)):
        segment_array = readcsv(segmentfiles[i])
        
        
        #for each heights csv, get the data into an array
        heights_array = readcsv(heightsfiles[i])
        
        
        #allocate space for the new segmented_heights array
        segmented_heights_array = np.zeros(size)
        adhesion_array = np.zeros([180000,4])        
        writeoutput(segment_array, end_dir, 'segment_array_orig_'+str(i+1)+'.csv')
        writeoutput(heights_array, end_dir, 'height_array_orig_'+str(i+1)+'.csv')
        
        
        #cycle through the segment array and if the value is not zero, go to the 
        #heights csv at the same position, and populate the segmented_heights array with that number
        adh_row=0
       
        for k in range(row):
            for j in range(column):
                rownum=k
                colnum=j
                
                if segment_array[rownum,colnum]>0 and heights_array[rownum,colnum]>0:
                    segmented_heights_array[rownum,colnum] = heights_array [rownum,colnum]
                    #print 'row '+str(rownum+1) +' column ' + str(colnum+1)
                    
                    #also populate the 4 column array with adhesion number, row number, column number
                    print rownum, colnum
                    adhesion_array[adh_row,0] = segment_array[rownum,colnum]
                    adhesion_array[adh_row,1] = rownum+1
                    adhesion_array[adh_row,2] = colnum+1
                    adhesion_array[adh_row,3] = segmented_heights_array[rownum,colnum]
                    
                    #increment the row of the adhesion_array that we will be writing to
                    adh_row+=1
                else:
                    continue  
        adhesion_array = adhesion_array[0:adh_row]
        
        #using the writeoutput function to write these arrays to csv files
        writeoutput(segmented_heights_array,end_dir,outputfiles[i])
        writeoutput(adhesion_array,end_dir,outputfiles2[i])
        
        #pull out the unique adhesions by number using the separate_adhesions function
        adh_nums= separate_adhesions(adhesion_array,0)
        
        #creates key for adhesion number to number of pixels in adhesion
        adh_pixel_count = {}
        
        #creates key for adhesion number to sum of all pixel heights
        adh_total_height = {}
        
        #creates key for adhesion number to average height
        adh_avg_heights = {}
        
        row_count = -1
        
        #allocate space that will later be populated with the average heights of each adhesion
        #there are two lists/arrays because the avg_heights list will be used to generate the histogram
        #the array will be used to output a csv, and is the repository for the more complete information
        #it includes pixel count and average height for each adhesion number
        adh_avg_heights_array = np.zeros([len(adh_nums),3])
        avg_heights =np.zeros([len(adh_nums),1])
        
        
        #cycle through each adhesion number in the unique list we created
        for adh in adh_nums:
                        
            for m in range(adh_row):
                
                #to deal with the way python treats arrays, we have to pull out each row separately
                #in order to access individual entries within the row
                #so we cycle through each row individually
                
                adhesion_row = adhesion_array[m]
                
                #first we compare the first value to our adhesion number of interest
                #if they're not the same, we move on
                if adhesion_row[0]!=adh:
                    continue
                
                #if the adhesion number is the same as the adhesion number of the row (first entry)
                else:
                    
                    # if there's already a mapping for this adhesion number
                    if adh in adh_pixel_count:
                        
                        #increment pixel count for the adhesion
                        adh_pixel_count[adh]+=1
                        
                        #pull old height number
                        adh_height_total_old = adh_total_height[adh]
                        adh_height_total_new = adh_height_total_old+ adhesion_row[3]
                        adh_total_height[adh]=adh_height_total_new
                    
                    #if there isn't a mapping for this adhesion number to pixel count, which means
                    #it's the first time we're encountering it
                    
                    else: 
                        
                        #we initialize each mapping
                        adh_pixel_count[adh]=1
                        #the total height is just the height of this one adhesion, which is the fourth entry in the row
                        #(recall that the row goes: adhesion number, x,y, height)
                        adh_total_height[adh]=adhesion_row[3]
                        
            #we assign the pixel count and total height values so we can operate on them later
            pixel_count = adh_pixel_count[adh]
            total_height = adh_total_height[adh]
            
            #use pixel count and heights sum to calculate average height
            adh_avg_heights[adh]=total_height/pixel_count
            
            #we set a minimum adhesion size, based on pixel count
            #this is defined in the input of this function, and is hard coded in right now
            if adh_pixel_count[adh]>=minAdh_size:
                #we're keeping track of how many adhesions pass this cutoff
                row_count +=1
                
                #populating the average heights lists/arrays
                
                #each array row goes: adhesion number, pixel count, average height of pixels in adhesion
                adh_avg_heights_array [row_count]= [adh, adh_pixel_count[adh],adh_avg_heights[adh]]
                #the list is just the average adhesion heights in ascending order of the adhesion number
                avg_heights[row_count]=[adh_avg_heights[adh]]
                
                
        adh_avg_heights_array=adh_avg_heights_array[0:row_count]
        avg_heights = avg_heights[0:row_count]
        print adh_avg_heights_array
        writeoutput(adh_avg_heights_array, end_dir, 'adh_avg_heights_'+str(i+1)+'.csv')

        #plot histogram
        fig=plt.figure()
        ax = fig.add_subplot(111)
        
        #we assign the data source (n) and number of bins 
        n,bins = np.histogram(avg_heights,50)
        left = np.array(bins[:-1])
        right = np.array(bins[1:])
        bottom = np.zeros(len(left))
        top = bottom + n
        XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

        # get the Path object
        barpath = path.Path.make_compound_path_from_polys(XY)
        
        # make a patch out of it
        patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
        ax.add_patch(patch)
        
        # update the view limits
        ax.set_xlim(left[0], right[-1])
        ax.set_ylim(bottom.min(), top.max())
        
        #save the image
        pylab.savefig(end_dir+'histogram'+str(i+1)+'.png')

    return 1

def separate_adhesions(array,index):
    #need to hard code the pattern
    
    #pulling out all the unique adhesion numbers within the adhesion array
    s=set([e[index] for e in array])
    adh_nums= []
    
    for adh_num in s:
        adh_nums=adh_nums + [adh_num]
    #sorts from smallest to largest
    adh_nums.sort()
       
    return adh_nums


def get_file_list(home_dir, unc_csv_name,heights_dir,subfolders,naming_scheme,end_dir):
    subfolders= subfolders.split(',')

    #BECAUSE I'M LAZY: JUST CODING IN THESE NUMBERS
    subfoldernum1 = 1 #cell number (the number of cells, not which number to start on, that's later)
    subfoldernum2 = 33 #time points
    
    file_list = []
    heights_list = []
    output_list =[]
    output_list2 =[]
    
    #iterates through each subfolder in the main home directory
    for i in range(subfoldernum1):
        for j in range(subfoldernum2):
            
            #this is unique to this set of data and may need to be changed
            #right now i is the cell number
            #j is the time, because that is how the data was saved (cell#t#)
            num1 = i+3
            num2 = j+1
            
            #necessary because there are 0s before the single digit numbers
            #be aware that this may change for data sets with even more timepoints (like >100)
            if num2 <10:
                csv_dir = home_dir+'0'+str(num2)+'/'+unc_csv_name 
            else:
                csv_dir = home_dir+str(num2)+'/'+unc_csv_name
            
            #home directory, plus subfolder, plus the actual name of the file (cell1t1_cH.csv, for example)
            
            heights_subfolder_name = subfolders[0]+ str(num1) + '/'+subfolders[1]+str(num2)+'/'+subfolders[0]+str(num1)+subfolders[1]+str(num2)
            heights_file =heights_dir+heights_subfolder_name+naming_scheme
            
            output_file = end_dir+subfolders[0]+str(num1)+subfolders[1]+str(num2)+naming_scheme
            output_file2=end_dir+subfolders[0]+str(num1)+subfolders[1]+str(num2)+'_adh_array.csv'
            
            file_list = file_list+[csv_dir]
            
            heights_list = heights_list+[heights_file]
            
            output_list = output_list + [output_file]
            output_list2=output_list2+[output_file2]
            
    return file_list,heights_list,output_list,output_list2,subfoldernum1, subfoldernum2


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
    
    #change the directory
    os.chdir(end_dir)
    #open the file 
    f = open(outputname, 'wb')
    #create a writer object
    writer = csv.writer(f)
    
    #iterate through each row of the array and populate a matrix 
    #corresponding to the data in the array
    for i in range(len(array)):
        writer.writerow(array[i,:])
    return 1


def main():
    if len(sys.argv) != 7:
        print 'Provide 1)directory,  \n 2)filename (of the csv file) \n 3) subfolders in the SAIM directory \n 4)directory where SAIM csvs are kept \n 5)namings scheme of SAIM csvs \n 6)final directory'
        sys.exit(1)
    print sys.argv
    
    home_dir = sys.argv[1]
    filename = sys.argv[2]
    folders = sys.argv[3]
    heights_dir = sys.argv[4]
    naming_scheme = sys.argv[5]
    end_dir = sys.argv[6]
    
    
    if cycle_through_files(home_dir,filename,heights_dir,folders,naming_scheme,end_dir) == 1:
        print 'done segmenting and writing csvs'
    else:
        print 'unable to complete command'
        sys.exit(1)
if __name__ == '__main__':
    main()