import sys
import numpy as np
from astropy.io import fits
import subprocess
import sewpy
import os

def configuration(head):
	gain=0.0
	pix_scale=1.0
	seeing_arc=1.2
	saturate=50000.0
	try:
		gain=head['GAIN']
	except KeyError:
		print "GAIN Key not found! using defaults!"
	try:
		pix_scale=head['PIXSCALE']
	except KeyError:
		print "PIXSCALE Key not found! using defaults!"
	try:
		seeing_arc=head['SEEING']*pix_scale
	except KeyError:
		print "SEEING Key not found! using defaults!"
	try:
		if(head['L1SAT']==False):
				saturate=head['L1COUNTS']+head['BACKGRD']+head['STDDEV']
	except KeyError:
		print "L1 Keys not found! using defaults!"			
	return gain,pix_scale,seeing_arc,saturate



science=list()
header=list()
if(len(sys.argv)<3):
	print "Input B and V filter images as argument!"
	exit()
try:
	img,hdr=fits.getdata("Fits/"+sys.argv[1],header=True)
	science.append(img)
	header.append(hdr)
	img,hdr=fits.getdata("Fits/"+sys.argv[2],header=True)
	science.append(img)
	header.append(hdr)
except:
	print "Files could not be opened!"


cat=list()
i=0
while(i<2):
	ret=configuration(header[i])
	sew = sewpy.SEW(params=["NUMBER","ALPHA_SKY","DELTA_SKY","X_IMAGE","Y_IMAGE","MAG_AUTO","MAG_MODEL","MAGERR_AUTO","MAGERR_MODEL"],config={"GAIN":ret[0],"SATUR_LEVEL":ret[3],"MAG_ZEROPOINT":float(sys.argv[i+3]),"PIXEL_SCALE":ret[1],"SEEING_FWHM":ret[2]})
	out = sew("Fits/"+sys.argv[i+1])
	print out["table"]
	cat.append(out["catfilepath"])
	print cat[i]
	i+=1

#p1=subprocess.Popen(["sextractor",""],stdin=PIPE,stdout=PIPE,stderr=PIPE)
#p2=subprocess.Popen(["sextractor",""],stdin=PIPE,stdout=PIPE,stderr=PIPE)

