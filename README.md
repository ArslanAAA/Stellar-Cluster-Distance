Estimating the Distance to Star Clusters using Automated Photometry

1. Abstract:
   Images are needed as a first step in classifying celestial objects based on their morphology. 
   To go beyond this initial stage of investigation and to set constraints on the models of stellar structure and evolution, we need quantitative information or measurements of the properties of the objects, collectively referred to as Photometric Analysis. 
   The aim of this project has been to use automated photometry on images of star clusters and determine the distance to them using the colour-magnitude diagram. This project was completed in two stages. 
   The first stage dealt with accurate measurement of the magnitudes of the stars in a cluster, either open or globular, in two different filters, visible (λcenter = 550 nm) and blue (λcenter = 440 nm). 
   This required acquisition of B and V filter CCD images in FITS format (Flexible Image Transport System), of the open cluster NGC 2420 and the globular cluster M80, reduction of the raw instrumental data followed by automated point spread function (PSF) photometry to calculate magnitudes tied to the standard UBV system. 
   The second stage involved estimating the distance to the two clusters using main sequence fitting technique. 
   Algorithms were designed, implemented and tested to automate the task of computation of the distance modulus of each, of which unsupervised learning for identifying main sequence data coupled with statistical data binning produced better results than those involving model fitting. 
   The values and their corresponding error bars were calculated. The standard values were found to be well within the error bars.

2. Package usage:
   Although the aim has been to completely automate the pipeline, there are certain intermediary steps which need to done.
   These involve the usage of certain tools as mentioned below.

3. Input:
    1. B and V filter CCD images in FITS format of the star cluster whose distance you want to measure
    2. Instrumental Zero Point Magnitude (either directly specified or estimated using Aladin or other such tools)

4. Output:
    Distance measuremnt in parsecs (3.26 light years)
   
5. Requirements:
    1. SExtractor version 2.19.5 for source catalogue extraction from given CCD images
    2. Topcat for matching two source catalogues with respect to RA, DEC values
    3. Aladin and Aperture Photometry Tool for calculation of zero point magnitude(optional)
    4. Xmgrace for plotting the colour magnitude and HR diagrams and for calculations(optional)

This works in two stages:

1. The first stage deals with accurate measurement of the B and V magnitudes of the stars in a cluster. 
    1. CCDreduction.py: reduce the raw CCD images
    2. photometry.py: generate SExtractor source catalogue files
    3. HR_plot.py: plot the colour-magnitude diagrams

2. The second stage deals with distance esimation using main sequence fitting technique. 
I have used a robust method for distance computation which uses a nearby stars catalogue for main sequence fitting.
-> use classdistance.py methods to compute distance using:
    1. Statistical binning
    2. Quadratic fitting
    3. Clustering and statistical binning

Procedure:

1. The input CCD images are first reduced. This step involves reduction for:
	1. Bias
	2. Dark current 
	3. Flat fielding
    
    This is done using CCDreduction.py if bias, dark and flat field frames are available. This also involves checking if the CCD images have the keywords which are relevant to reduction and whether they are set or not.
    
2. The reduced CCD images are then used to extract source catalogues using SExtractor. This is done using photometry.py, by utilising all the relevant header information for extraction. This requires specifying the zero point magnitude (MAG_ZERO), for converting instrumental magnitudes to standard magnitudes.
   We assume uniform seeing throughout the image field of view and thus the instrumental magnitudes can be converted to standard magnitudes using the zero point magnitude calculated using Aladin and Aperture Photometry Tool. 
   This requires calculating the instrumental magnitude using the APT sky algorithm, sky-annulus median subtraction for a few known stars in the cluster(use RA, DEC to identify a star) and subtracting that from the standard magnitude of those stars obtained using Aladin. This will be done for both the filters, to calculate zero point magnitudes for both the filter images.

3. Then, we need to match the two source catalogues, which are extracted from images of the same cluster, at two different times, in two different clusters. This is done by matching the stars in the cluster according to the RA(right ascension), DEC(declination) values. This is done using Topcat which is a very powerful tool for matching catalogues with respect to RA, DEC or coordinates in general. 
   An object is defined as being the same one in both tables if the co-ordinates in both rows are "similar", for instance if the difference between the positions indicated by RA and Dec columns differ by no more than a specified angle on the sky. Matching rows to produce the join requires you to specify the criteria for rows in both tables to refer to the same object and what to do when one is found.
   Thus, we use the Topcat Sky Algorithm with a max error of 3*Seeing(arcsecs) as an error measure.

4. The next step is to use the matched catalogue for plotting the colour-magnitude diagram, followed by distance modulus computation.

Distance in parsecs = 10**(distance_modulus/5 + 1)
