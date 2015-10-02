#!/usr/local/share/bin/EPD/python -tt

'''
Created on Jul 17, 2012

@author: gou
'''

import sys
import csv

def draw_heat_map():

    heights_map = 1


    return heights_map


def get_segmentation(dir,filename):
    
    segmentation_list =dir+filename
    pixelslist = csv.reader(open(segmentation_list,'rb'))
    return pixelslist

def make_image(heights_map, other_image):
    
    return 1

def main():
    #if len(sys.argv) != 3:
        #print 'Provide 1)directory,  \n 2)filename (of the csv file)'
        #sys.exit(1)
    home_dir = sys.argv[1]
    filename = sys.argv[2]
    #if make_image(home_dir, subfolders,searchterm,newdir)== 1:
        #print 'Adding maps done' 
    #else: 
        #print 'Unable to complete command'
    list = get_segmentation(home_dir,filename) 
    for row in list:
        print row
    
    sys.exit(1)
if __name__ == '__main__':
    main()