#!/usr/local/share/bin/EPD/python -tt
#!C:\EPD\python -tt
'''
Created on Oct 21, 2012

@author: gou
'''
# getAdhesionImage.py

import Image, ImageDraw
import sys
import random
import numpy as np

#Import a CSV File
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

def readFile(filename):
# Open the file, read everything, and then close the file
    f = open(filename)
    lines = f.readlines()
    f.close()
    print 'opened'
    
    # Construct a matrix from the file
    
    matrix = []
    
    # Each line gives us a column in the matrix
    for line in lines:
    
        # Get the column vector of strings by splitting the line according to commas
        rowStr = line[:-1].split(',')
        
        # Go through the array of number strings, create int's from them, and put the
        # new row in the matrix
        row = [];
        
        for n in rowStr:
            # Generate the row by successively adding the numbers to an array
            row += [int(n)]
            
            # Add the row to the matrix
            matrix += [row]
    
    return matrix



# Create an image from a matrix
# scalingFactor is the scaling factor for the image. For example, 2 will 
# create an image twice as big
# minClumpSize is the minimum size of a group of pixels to draw
def makeImage(imgName, matrix, scalingFactor=2, minClumpSize=5):

    # Get the dimensions of the image from the matrix
    height = len(matrix)
    width = len(matrix[0])
    print height
    print width
    
    # Initalize the image with a black background
    img = Image.new('RGB', (scalingFactor*width, scalingFactor*height), 'black')
    
    # This is the object we'll use to draw onto the image
    draw = ImageDraw.Draw(img)
    
    
    # This is a map of all the colors we use
    # For now, we know that 0 (no adhesion) shouldn't be colored
    # We'll fill in this map as we go
    colormap = {0:(0, 0, 0)};
    
    # This is a map that holds that sum of the x-y position of the pixels, 
    # by adhesion value
    # We'll use this later to find the average position of the adhesion
    positions = {}
    
    # This is a map that simply keeps track of how many pixels were in
    # the adhesion.
    # This is also useful for computing the average position, and also
    # deciding whether or not to draw the adhesion.
    val_n = {}
    
    
    
    # Let's first figure out what clumps are there
    
    
    # Iterate over every entry in the matrix
    for x in xrange(height):
        for y in xrange(width):
    
            # This is the value in the image
            val = matrix[x][y]
                       
            
            # We don't care about 0's, so if we have a 0, ignore it and move on
            #if val is 0:
            if val==0:
                #print x,y
                continue
            
            # If we haven't yet assigned a color for this value, do so now
            if not val in colormap:
                # Construct a random RGB color
                # I'm using a baseline of 100 so that we don't draw black
                # pixels, which would be hard to see on a black background
                
                R = 100+random.randint(0, 156)
                G = 100+random.randint(0, 156)
                B = 100+random.randint(0, 156)
                color = (R, G, B)
                
                # Store the color in our color for use later
                colormap[val] = color
            
            # We need to keep track of the positions of all the pixels, so we can know
            # where to draw the text later
            pos = (x, y)
            if val in val_n:
                # If we already came across this value, increment the count
                # and update the sum of the x and y coordinates accordingly
                val_n[val] += 1
    
                # Grab the old position from the map, add in the new position
                # by x and y component, and stick it back in the map
                oldPos = positions[val]
                newPos = ( (pos[0]+oldPos[0]), (pos[1]+oldPos[1]) )
                positions[val] = newPos;
                
        
            else:
                # Otherwise, create the intial entry
                
                # This is the first one, so we only have 1 of this value    
                val_n[val] = 1
                # This is the only position in the map so far, so just stick it in
                positions[val] = pos
            
    
    # Now let's go through and actually draw in things
    for x in xrange(height):
        for y in xrange(width):
    
            # Grab the value at this position    
            val = matrix[x][y]
            
            # We'll only draw something if it's nonzero (so, actually an adhesion)
            # and it's bigger than the minimum size we're interested in
            if val!=0 and val_n[val] > minClumpSize:
            
                # Grab the color we had previously assigned
                color = colormap[val]
                
                # Generate the box, according to the scaling factor
                # ( (top left corner y x), (bottom right corner y x) )
                # y and x are swapped, because that's how they're drawn
                topLeft = (scalingFactor*y, scalingFactor*x) 
                bottomRight = (scalingFactor*(y+1), scalingFactor*(x+1))
                box = (topLeft, bottomRight)
                
                # Now draw the box that represents this matrix entry
                draw.rectangle(box, fill=color)
    
    
    # Now let's actually label things
    for val in val_n.keys():
    
        # How many pixels of this value were there?
        n = val_n[val]
        
        # Only label if we had at least the minimum size
        if n > minClumpSize:
            print val
            # Compute the average position of all of the pixels we care about
            pos = positions[val]
            avg_pos = ( scalingFactor*pos[0]/n, scalingFactor*pos[1]/n)
            
            # Draw the label, with the text just being the adhesion number, and the
            # position being in the middle of the adhesion itself
            draw.text( (avg_pos[1], avg_pos[0]), str(int(val)))
                    
    # Save the image to the given filename
    img.save(imgName, 'BMP')
    return 1

def main():
    if len(sys.argv)!=3:
        print 'getAdhesionImage.py'
        print ''
        print 'Usage:'
        print ''
        print 'python getAdhesionImage.py csvFilename imageFilename'
        print ''
        print 'Generates an image of labeled adhesions from an input CSV.'
        print 'Note that any existing image will be overwritten!'
        sys.exit()

    # Read in the CSV file
    #matrix = readFile(sys.argv[1])
    matrix = readcsv(sys.argv[1])
    
    # Create the image
    if makeImage(sys.argv[2], matrix)==1:
        print 'complete'
if __name__ == '__main__':
    main()