import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.modeling import models,fitting
#nearby=raw_input('Enter the filename/URL containing nearby stars catalogue : ')
#read data into an astropy ascii table
dat_n = ascii.read('NearbyStars.tsv')

#cluster=raw_input('Enter the filename/URL containing cluster catalogue : ')
#read data into an astropy ascii table
dat_c = ascii.read('NGC2420.dat')

color_n = dat_n['B-V']       #array of values with column B-V
Mv_n = dat_n['Mv']           #array of values with column Mv

nearby=[]

for i in range(len(color_n)):
	if color_n.mask[i]==False:    
    		nearby.append((np.array(color_n)[i],np.array(Mv_n)[i]))

nearby_struct=np.array(nearby,dtype=[('color',float),('Mag',float)])

B = dat_c['MAG_AUTO_1']
V = dat_c['MAG_AUTO_2']
color_NGC2420= B - V

sigmaB = dat_c['MAGERR_AUTO_1']
sigmaV = dat_c['MAGERR_AUTO_2']

sigma= np.sqrt(sigmaB**2 + sigmaV**2)

#plt.errorbar(color,V,xerr=sigma,yerr=sigmaV,linewidth=1,color='black',fmt="none")
#plt.scatter(color,V,marker='.')
NGC2420=[]

for i in range(len(color_NGC2420)):
	NGC2420.append((np.array(color_NGC2420)[i],np.array(V)[i]))

NGC2420_struct=np.array(NGC2420,dtype=[('color',float),('Mag',float)])

def split(arr,cond):
	return [arr[cond],arr[~cond]]

#NGC2420V=V.tolist()
#NGC2420color=color_NGC2420.tolist()
#NGC2420=np.vstack((NGC2420color,NGC2420V)).T

#Mv=Mv_n.tolist()
#color=color_n.tolist()
#nearby=np.vstack((color,Mv)).T
min_NGC2420=np.amin(NGC2420_struct['color'])
min_nearby=np.amin(nearby_struct['color'])
a=[min_nearby,min_NGC2420]

max_NGC2420=np.amax(NGC2420_struct['color'])
max_nearby=np.amax(nearby_struct['color'])
b=[max_nearby,max_NGC2420]

if max(a)==min_NGC2420:
	nearby_struct=split(nearby_struct,nearby_struct['color']<min_NGC2420)[1]
else:
	NGC2420_struct=split(NGC2420_struct,NGC2420_struct['color']<min_nearby)[1]

if min(b)==max_NGC2420:
	nearby_struct=split(nearby_struct,nearby_struct['color']>max_NGC2420)[1]
else:
	NGC2420_struct=split(NGC2420_struct,NGC2420_struct['color']>max_nearby)[1]

print NGC2420_struct
print nearby_struct

#min=max(min(NGC2420color),min(color))
#max=min(max(NGC2420color),max(color))
#print min,max

#okay=[tuple(p) for p in NGC2420.tolist()]
#final=np.array(okay,dtype=[('color',float),('Mag',float)])
sorted_NGC2420=np.sort(NGC2420_struct,order='color')

#okay_n=[tuple(p) for p in nearby.tolist()]
#final_n=np.array(okay_n,dtype=[('color',float),('Mag',float)])
sorted_nearby=np.sort(nearby_struct,order='color')

delta=(min(b)-max(a))/10

list=[]

for i in range(1,10):
    list.append(max(a)+i*delta)

#s=sorted.view(np.float64).reshape(sorted.shape+(-1,))

split_at_NGC2420=sorted_NGC2420['color'].searchsorted(list)
NGC2420_split=np.split(sorted_NGC2420,split_at_NGC2420)

split_at_nearby=sorted_nearby['color'].searchsorted(list)
nearby_split=np.split(sorted_nearby,split_at_nearby)

#counts,edges=np.histogram(sorted_NGC2420['color'],bins=list)

#bin_middles= (edges[:-1]+edges[1:]) /2

#weights = np.array(range(len(counts)) / sum(range(len(counts))
weights_NGC2420=[]
for i in range(10):
	weights_NGC2420.append(NGC2420_split[i].size)
	
weights_nearby=[]
for i in range(10):
	weights_nearby.append(nearby_split[i].size)

y_NGC2420=[]
for i in range(10):
	y_NGC2420.append(np.mean(NGC2420_split[i]['Mag']))

y_nearby=[]
for i in range(10):
	y_nearby.append(np.mean(nearby_split[i]['Mag']))

y=np.asarray(y_NGC2420)-np.asarray(y_nearby)

y_avgNGC2420=np.average(y,weights=weights_NGC2420)
y_avgnearby=np.average(y,weights=weights_nearby)

print y_avgNGC2420,y_avgnearby

p1=models.Polynomial1D(2)
pfit=fitting.LinearLSQFitter()

model_NGC2420=pfit(p1,NGC2420_struct['color'],NGC2420_struct['Mag'])
model_nearby=pfit(p1,nearby_struct['color'],nearby_struct['Mag'])

distance=model_NGC2420.c0.value-model_nearby.c0.value

print distance

plt.figure()
plt.gca().invert_yaxis()
plt.scatter(color_n,Mv_n,marker='o')

plt.scatter(NGC2420_struct['color'],NGC2420_struct['Mag'],marker='o')
plt.scatter(nearby_struct['color'],nearby_struct['Mag'],marker='.')
plt.show()




