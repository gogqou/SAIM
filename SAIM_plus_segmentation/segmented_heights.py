#!/usr/local/share/bin/EPD/python -tt

'''
Created on Jul 17, 2012

@author: gou
'''

import sys
import csv

def get_segmentation(home_dir,filename):
    
    segmentation_list =home_dir+filename
    pixelslist = csv.reader(open(segmentation_list,'rb'))
    return pixelslist

def write_heights(pixels_list, heights_csv,naming_scheme,end_dir):
    ad_num = 0
    
    for row in pixels_list:
        ad_num +=1
        
        for i in range(len(row)):
            print ad_num
    print str(ad_num)+' adhesions processed'
    return 1

def cycle_through_files(home_dir,pixels_list, heights_csv, naming_scheme, end_dir):
    #set up a list of each time point in the experiment--and if necessary, each experiment, if processing more than one set of data
    
    #iterate through each time point and experiment and correlate the segmented pixels with calculated heights from
    #the SAIM analysis
    
    data_files = []
    for data_file in data_files:
        a=write_heights(pixels_list,heights_csv,naming_scheme,end_dir)
    
    
    
    return 1


def main():
    if len(sys.argv) != 6:
        print 'Provide 1)directory,  \n 2)filename (of the csv file) \n 3)directory where SAIM csvs are kept \n 4)namings scheme of SAIM csvs \n 5)final directory'
        sys.exit(1)
    home_dir = sys.argv[1]
    filename = sys.argv[2]
    heights_dir = sys.argv[3]
    naming_scheme = sys.argv[4]
    end_dir = sys.argv[5]
    #if make_image(home_dir, subfolders,searchterm,newdir)== 1:
        #print 'Adding maps done' 
    #else: 
        #print 'Unable to complete command'
    pixels_list = get_segmentation(home_dir,filename) 
    if write_heights(pixels_list, heights_dir,naming_scheme,end_dir) == 1:
        print 'done segmenting and writing csvs'
    else:
        print 'unable to complete command'
        sys.exit(1)
if __name__ == '__main__':
    main()