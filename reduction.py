import numpy as *
from astropy.io import fits

biases=array([fits.getdata("bias-%d.fit" % n ) for n in range(1,bias)])
bias=median(biases,axis=0)
darks=array([fits.getdata("dark-%d.fit" % n ) for n in range(1,dark)])
dark=median(darks,axis=0)
flatsb=array([fits.getdata("flat-b-%d.fit" % n ) for n in range(1,flat_b)])
flatb=median(flatsb,axis=0)
flatsv=array([fits.getdata("flat-v-%d.fit" % n ) for n in range(1,flat_v)])
flatv=np.median(flatsv,axis=0)

scienceb=array([fits.getdata("" %n ) for n in range(1,sci_b)])

finalb = ((scienceb-dark)/(flatb-bias))*np.mean(flatb-bias)
B=median(finalb,axis=0)


sciencev=array([fits.getdata("" %n ) for n in range(1,sci_v)])

finalv = ((sciencev-dark)/(flatv-bias))*np.mean(flatv-bias)
V=median(finalv,axis=0)

fits.write('Reduction.fits')
