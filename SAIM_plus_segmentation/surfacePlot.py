#!/usr/local/share/bin/EPD/python -tt

'''
Created on 2010-11-23

@author: weaver
'''

import numpy as np
import mayavi.mlab as mlab


def surfacePlot(Hmap, nrows, ncols, xyspacing, zscale, name, hRange, file_path, lutfromfile, lut, lut_file_path, colorbar_on, save, show):
    
    #Create a grid of the x and y coordinates corresponding to each pixel in the height matrix
    x, y = np.mgrid[0:ncols*xyspacing:xyspacing, 0:nrows*xyspacing:xyspacing]
    
    #Create a new figure
    mlab.figure(size=(1000, 1000))  
    
    #Set the background color if desired
    #bgcolor=(0.16, 0.28, 0.46)

    #Create the surface plot of the reconstructed data
    plot = mlab.surf(x, y, Hmap, warp_scale=zscale, vmin=hRange[0], vmax=hRange[1], colormap=lut)
    
    #Import the LUT from a file if necessary
    if lutfromfile:
        plot.module_manager.scalar_lut_manager.load_lut_from_file(lut_file_path)
        

    #Draw the figure with the new LUT if necessary
    mlab.draw()
    
    #Zoom in to fill the entire window
    f = mlab.gcf()
    f.scene.camera.zoom(1.05)

    #Change the view to a top-down perspective
    mlab.view(270,0)
    
    #Add a colorbar if indicated by colorbar_on (=True)
    if colorbar_on:
        #mlab.colorbar(title='Height (nm)', orientation='vertical')
        mlab.colorbar(orientation='vertical')

    
    #Save the figure if indicated by save (=True)
    if save:
        mlab.savefig(file_path, size=(1000,1000))
        if show == False:
            mlab.close()
    
    #Keep the figure open if indicated by show (=True)        
    if show:
        mlab.show()
        