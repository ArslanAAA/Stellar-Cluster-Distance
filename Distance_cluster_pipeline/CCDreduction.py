import numpy as np
from astropy.io import fits

biases = np.array([fits.getdata("bias-%d.fit" % n ) for n in range(1, bias)])
bias = np.median(biases,axis=0)
darks = np.array([fits.getdata("dark-%d.fit" % n ) for n in range(1, dark)])
dark = np.median(darks,axis=0)
flatsb = np.array([fits.getdata("flat-b-%d.fit" % n ) for n in range(1, flat_b)])
flatb = np.median(flatsb,axis=0)
flatsv = np.array([fits.getdata("flat-v-%d.fit" % n ) for n in range(1, flat_v)])
flatv = np.median(flatsv,axis=0)

scienceb = np.array([fits.getdata("science-b-%d" %n , header=True ) for n in range(1, sci_b)])
scienceb_data = [ x[0] for x in scienceb ]
scienceb_header = [ x[1] for x in scienceb ]

finalb = ((scienceb_data - dark)/ (flatb-bias)) * np.mean(flatb - bias)
B=median(finalb, axis=0)

fits.writeto("Reduced_B.fits", B)

sciencev=array([fits.getdata("science-v-%d" %n , header=True ) for n in range(1, sci_v)])
sciencev_data = [ x[0] for x in sciencebv ]
sciencev_header = [ x[1] for x in sciencev ]

finalv = ((sciencev_data - dark)/(flatv - bias)) * np.mean(flatv - bias)
V=median(finalv, axis=0)

fits.writeto("Reduced_V.fits", V)
#[fits.writeto('ReducedV-%d.fits' %n , finalv[n], sciencev_header[n]) for n in range(0, len(finalv))]
