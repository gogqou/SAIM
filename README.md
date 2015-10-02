# SAIM
Analysis code for SAIM (scanning angle interference microscopy) analysis

MultiAngleFLIC--source code for executable that uses Intel libraries to execute a lot of the behind the scenes math and optics. Base authored by Matthew Paszek and published in Nature Methods. My work is built on top of his framework to implement the ability to handle more models, including the addition of a material of variable refractive index on top of the silicon imaging surface.

Mapping: post processing to output visualization of SAIM data. 


Original paper:
http://www.nature.com/nmeth/journal/v9/n8/full/nmeth.2077.html

This method is build upon multi angle FLIC (fluorescence interference contrast microscopy). The imaging method takes advantage of interference between light bouncing off of a silicon oxide surface and a fluorophore resting <1000nm above the surface to calculate its height from the surface. The MultiAngleFLIC code takes  microscope settings and material properties to perform this calculation from a series of images collected at different angles of incident light. It outputs estimates (model fits) for height of the fluorophore, which is then taken by Mapping to draw a height map on top of the input images. 
