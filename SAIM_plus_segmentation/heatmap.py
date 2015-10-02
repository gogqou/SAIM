
#!/usr/local/share/EPD/bin/python -tt

'''
Created on Jun 23, 2012

@author: paszekm
'''

'''
Created on 2010-09-08

@author: matt
'''

import numpy as np
import pylab as py
import Image
import csv
from scipy.misc import imsave
import os as os

#convert image into array according to im.mode - Still in development
def image2arrayFast(im):
    print im.mode
    if im.mode in ("L",):
        return np.reshape(np.fromstring(im.tostring(), np.uint8), (im.size[1],im.size[0]))
    elif im.mode in ("I;16",):
        return np.reshape(np.fromstring(im.tostring(), np.uint16), (im.size[1],im.size[0]), order='A')
    elif im.mode in ("I;16B",):
        Temp = np.reshape(np.fromstring(im.tostring(), np.uint16), (im.size[1],im.size[0]))
        Temp2 = np.float64(Temp)
        return Temp2
        #return np.reshape(np.fromstring(im.tostring(), np.uint16), (im.size[1],im.size[0]))
    else:
        print "here"

#convert image into array
def image2array(im):
    pix = im.load()
    a = np.zeros((im.size[1], im.size[0]))
    for i in range(im.size[1]):
        for j in range(im.size[0]):
            a[i,j] = pix[j,i]
            if pix[j,i] < 0:
                print pix[j,i]
    return a
            

#convert array into image according to a.dtype
def array2image(a):
    if a.dtype == np.uint8: 
        mode = "L"
    elif a.dtype == np.uint16: 
        mode = "I;16"
    else:
        raise ValueError, "unsupported image mode"
    return Image.fromstring(mode, (a.shape[1],a.shape[0]), a.tostring())


#The file path to where the raw data files are stored
#path = 'C:\\Users\\paszekM\Desktop\\SAIMforSAM\\output\\'
path = 'E:\\062812 JRT3\\AnalysisGreen\\output\\'
#path = 'E:\\062812 JRT3\\AnalysisRed\\output\\'
     
#The pixel size in nanometers
xyspacing = 110  #was 0.45

#modify this parameter if you wish to change the z-scaling (height) relative to the xy-scaling
zscale = 0 #40 #80 for microtubule

#The range of the colormap for the z heights
#plotRange = [70,350] #adhesion [40,120]  #microtubule [70,350]
plotRange = [30, 100]

#The range of z heights to include in the graphing
#hRange = [10,430]#adhesion [20,140] #microtubule [20,450]
hRange = [20,100]

#Color map file path
lut_file_path = "C:\\Users\\paszekM\\Desktop\\SAIMforSAM\\LUT14.txt"

#Set to true if you want a colorbar on your image.  Otherwise set to False
colorbar_on = True

#View figure - set this option to true if you would like to view and modify the surface plot and then manually save it
view = False

#Save Figure? - Set this option to false if you plan on viewing the surface plot, manually adjusting it, and saving it
save = True

#The desired file extension for the save surface plot - can be .jpg, .png, .tif, etc.  The reconstruction will be saved as this file type
surface_ext = '.png'

#The desired number of labels on the color scalebar
nlabels = 8

#The list of files that you would like to open and create surface plots of 
#name_list=['gluter_control001_c','gluter_control003_c','gluter_control005_c','gluter_control007_c','gluter_control009_c','gluter_control011_c','gluter_control013_c','gluter_control015_c','gluter_control017_c','gluter_control019_c','gluter_control021_c','gluter_control023_c','gluter_control025_c','gluter_control027_c','gluter_control029_c','gluter_control031_c','gluter_control033_c','gluter_control035_c','gluter_control037_c','gluter_control039_c','gluter_control041_c','gluter_control043_c','gluter_control045_c']
#name_list=['gluter_clasp1003_c'] 
#name_list=['gluter_control001_c', 'gluter_control003_c','gluter_control005_c','gluter_control007_c','gluter_control009_c','gluter_control011_c','gluter_control013_c','gluter_control015_c','gluter_control017_c','gluter_control019_c','gluter_control021_c','gluter_control023_c','gluter_control025_c','gluter_control027_c','gluter_control029_c','gluter_control031_c','gluter_control033_c','gluter_control035_c','gluter_control037_c','gluter_control039_c','gluter_control041_c','gluter_control043_c','gluter_control045_c','gluter_control047_c','gluter_control049_c']
#name_list=['gluter_clasp1001_c','gluter_clasp1003_c','gluter_clasp1005_c','gluter_clasp1007_c','gluter_clasp1009_c','gluter_clasp1011_c','gluter_clasp1013_c','gluter_clasp1015_c','gluter_clasp1017_c','gluter_clasp1019_c','gluter_clasp1021_c','gluter_clasp1023_c','gluter_clasp1025_c','gluter_clasp1027_c','gluter_clasp1029_c','gluter_clasp1031_c','gluter_clasp1033_c','gluter_clasp1035_c','gluter_clasp1037_c','gluter_clasp1039_c','gluter_clasp1041_c','gluter_clasp1043_c','gluter_clasp1045_c','gluter_clasp1047_c']
#name_list=['gluter_clasp2003_c','gluter_clasp2005_c','gluter_clasp2007_c','gluter_clasp2009_c','gluter_clasp2011_c','gluter_clasp2013_c','gluter_clasp2015_c','gluter_clasp2017_c','gluter_clasp2019_c','gluter_clasp2021_c','gluter_clasp2023_c','gluter_clasp2025_c','gluter_clasp2027_c','gluter_clasp2029_c','gluter_clasp2031_c','gluter_clasp2033_c','gluter_clasp2035_c','gluter_clasp2037_c','gluter_clasp2039_c','gluter_clasp2041_c','gluter_clasp2043_c','gluter_clasp2045_c']
name_list = ['c']

#The number of columns in each of the image files above
ncols=512

#Make any necessary folders
file_path = path + "maps"
if not os.path.isdir(file_path):
    os.makedirs(file_path)
    
#Open the LUT
LUT = np.genfromtxt(lut_file_path)

#Make the colorbar
fig = py.figure(figsize=(8, 0.8))
fig.set_facecolor('black')   
ax1 = fig.add_axes([0.05, 0.1, 0.9, 1])
x = np.arange(255)
Bar = np.zeros((15,255,3))
for i in x:
    Bar[:,i,0] = LUT[i,0]*255
    Bar[:,i,1] = LUT[i,1]*255
    Bar[:,i,2] = LUT[i,2]*255

file_path = path + 'maps//Colorbar.png'
imsave(file_path, Bar)
colorbar = py.imread(file_path)
ax1.imshow(colorbar)
ax1.axis("off")
inc_scale = (plotRange[1]-plotRange[0])/(nlabels-1)
inc_pos = 0.9/(nlabels-1)
for i in range(nlabels):
    fig.text(0.035+i*inc_pos, 0.07, str(np.round(plotRange[0]+i*inc_scale)), fontsize=16, color='white')

file_path = path + 'maps//Scalebar.png'
py.savefig(file_path, dpi=600, facecolor='black')
fig.clear()

#Make the maps of the images
for name in name_list:
    
    #Open the csv file containing the best-fit heights output from the flic fitting software
    file_path = path + name + '_H.csv'
    reader = csv.reader(open(file_path, 'rb'), quoting=csv.QUOTE_NONNUMERIC)
    
    #Convert the first row of the file to a numpy array.  Note. the csv file should only have a single row
    for row in reader:
        H = np.array(row)

    #The best-fit heights in the output file were stored in row-major format. Reshape as a two-dimensional matrix             
    Hmap = H.reshape(-1,ncols)
    nrows, ncols = Hmap.shape
    
    
    print 'Processing: ', name

    #Exclude out of range data
    for i in range(nrows):
        for j in range(ncols):
            
            if Hmap[i,j] < plotRange[0]:
                if Hmap[i,j] >= hRange[0]:
                    Hmap[i,j] = plotRange[0]
                else:
                    Hmap[i,j] = -1
            
            if Hmap[i,j] > plotRange[1]:
                if Hmap[i,j] <= hRange[1]:
                    Hmap[i,j] = plotRange[1]
                else:
                    Hmap[i,j] = -1
        
    #print Hmap.shape
    
    #Elimate any negative numbers in the array for exporting
    HmapNoZero = Hmap
    for i in range(nrows):
        for j in range(ncols):
            if HmapNoZero[i,j] < 0:
                HmapNoZero[i,j] = 0
    
    #Convert Hmap into an RGB image
    h_unit = (plotRange[1] - plotRange[0])/255.0
    RGB = np.zeros((nrows,ncols,3))
    for i in range(nrows):
        for j in range(ncols):
            #print Hmap[i,j], h_unit
            if Hmap[i,j] == plotRange[0]:
                index = 1
            elif (Hmap[i,j] - plotRange[0]) > 0: 
                index = np.round((Hmap[i,j]-plotRange[0])/h_unit)
            else:
                index = 0
            RGB[i,j,0] = LUT[index,0]*255
            RGB[i,j,1] = LUT[index,1]*255
            RGB[i,j,2] = LUT[index,2]*255
  
    
    #Save the mapped image
    RGB8 = np.uint8(RGB)
    file_path = path + 'maps//' + name + '_mapIS2' + surface_ext
    imsave(file_path, RGB8)
    #file_path = path + 'maps2/' + name + '_map' + surface_ext
    
    #Multiply the Map image times the intensity image
    #Open the intensity image
    file_path = path + name + 'T.png'
    im = Image.open(file_path)
    Intensity = image2array(im)
    I_max = np.max(np.max(Intensity))
    Intensity = Intensity/I_max
    
    for i in range(nrows):
        for j in range(ncols):
            RGB[i,j,0] = RGB[i,j,0]*Intensity[i,j]
            RGB[i,j,1] = RGB[i,j,1]*Intensity[i,j]
            RGB[i,j,2] = RGB[i,j,2]*Intensity[i,j]
    
    #Save the intensity mapped image
    RGB8 = np.uint8(RGB)
    file_path = path + 'maps//' + name + '_mapIS2' + surface_ext
    imsave(file_path, RGB8)
       
    print "finished"