Estimating the Distance to Star Clusters using Automated Photometry

The aim of this package is to use automated photometry on images on star clusters and determine the distance to them using the colour-magnitude diagram.

Input:
   1. B and V filter CCD images in FITS format of the star cluster whose distance you want to measure
   2. Instrumental Zero Point Magnitude (either directly specified or estimated using Aladin or other such tools)

Output:
   Distance measuremnt in parsecs (3.26 light years)

This works in two stages:

1. The first stage deals with accurate measurement of the B and V magnitudes of the stars in a cluster. 
    1. reduction.py: reduce the raw CCD images
    2. image_red.py: generate SExtractor source catalogue files
    3. HR_plot.py: plot the colour-magnitude diagrams

2. The second stage deals with distance esimation using main sequence fitting technique. 
I have used a robust method for distance computation which uses a nearby stars catalogue for main sequence fitting.
-> use classdistance.py methods to compute distance using:
   1. Statistical binning
	 2. Quadratic fitting
	 3. Clustering and statistical binning

  
