import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.modeling import models,fitting
from math import sqrt
import random,sys
import time

class distance():
	
	def __init__(self,nearby,cluster):			#read the cluster and nearby stars data into ascii tables
		self.nearby= ascii.read(nearby)
		self.cluster= ascii.read(cluster)
		self.flag=False
		self.call()
		
	def data(self):								#read in required data
		self.color_n = self.nearby['B-V']       #array of values with column B-V
		self.Mv_n = self.nearby['Mv']           #array of values with column Mv
		B = self.cluster['MAG_AUTO_1']			
		self.V = self.cluster['MAG_AUTO_2']		#array of cluster magnitude in V band 
		self.color_c= B - self.V				#array of cluster color i.e. B-V
		sigmaB = self.cluster['MAGERR_AUTO_1']
		self.sigma_y = self.cluster['MAGERR_AUTO_2']		#standard deviation of cluster magnitude in V band i.e. y axis error 
		self.sigma_x= np.sqrt(sigmaB**2 + self.sigma_y**2)	#std dev of B-V i.e. x axis error
	
	def coordinates(self,data_x,data_y):		#tuples (x,y) :coordinates numpy structure ('color','Mag')
		list=[]
		flag_x=0
		flag_y=0
		if type(data_x).__name__ =='MaskedColumn':	#masked values are to be ignored
				flag_x=1
		if type(data_y).__name__=='MaskedColumn':
				flag_y=1
		try:
				for i in range(len(data_x)):
					if flag_x!=1 and flag_y!=1:
							list.append((np.array(data_x)[i],np.array(data_y)[i]))
					elif data_x.mask[i]==False and data_y.mask[i]==False:
    						list.append((np.array(data_x)[i],np.array(data_y)[i]))
		except AttributeError:
    			if flag_x!=1:
    					if data_y.mask[i]==False:
    							list.append((np.array(data_x)[i],np.array(data_y)[i]))
    			else:
    					if data_x.mask[i]==False:
    							list.append((np.array(data_x)[i],np.array(data_y)[i]))
		list_struct=np.array(list,dtype=[('color',float),('Mag',float)])
		return list_struct
    		   	
	def split(self,arr,cond):			#split an array at a given condition
		return [arr[cond],arr[~cond]]

	def min(self,nearby_struct,cluster_struct):		#find minimum of the two datasets and split
		min_n=np.amin(nearby_struct['color'])
		min_c=np.amin(cluster_struct['color'])
		a=[min_n,min_c]
		
		if max(a)==min_c:
			nearby_struct=self.split(nearby_struct,nearby_struct['color']<min_c)[1]
		else:
			cluster_struct=self.split(cluster_struct,cluster_struct['color']<min_n)[1]
			
		return max(a),nearby_struct,cluster_struct
			
	def max(self,nearby_struct,cluster_struct):		#find maximum of two datasets and split 
		max_n=np.amax(nearby_struct['color'])
		max_c=np.amax(cluster_struct['color'])
		a=[max_n,max_c]
		
		if min(a)==max_c:
			nearby_struct=self.split(nearby_struct,nearby_struct['color']>max_c)[1]
		else:
			cluster_struct=self.split(cluster_struct,cluster_struct['color']>max_n)[1]
		return min(a),nearby_struct,cluster_struct
	
	def binning(self,max,min,bins=5):
		delta=(max-min)/bins
		list=[]
		for i in range(1,bins):
			list.append(min+i*delta)
    		
		split_at_c=self.inter_sorted_cluster['color'].searchsorted(list)
		cluster_split=np.split(self.inter_sorted_cluster,split_at_c)
		split_at_n=self.inter_sorted_nearby['color'].searchsorted(list)
		nearby_split=np.split(self.inter_sorted_nearby,split_at_n)
		weights_c=[]
		weights_n=[]
		for i in range(bins):
				weights_c.append(cluster_split[i].size)
				weights_n.append(nearby_split[i].size)
		return cluster_split,nearby_split,weights_c,weights_n
		  
	def intersection(self,nearby_struct,cluster_struct):
		minima,nearby_struct,cluster_struct=self.min(nearby_struct,cluster_struct)
		maxima,nearby_struct,cluster_struct=self.max(nearby_struct,cluster_struct)
		return (nearby_struct,cluster_struct,minima,maxima)
 		
	def distance_bins(self,bins=5,plot=False):
		nearby_struct,cluster_struct,minima,maxima=self.intersection(self.nearby_struct,self.cluster_struct)
		self.inter_sorted_cluster=np.sort(cluster_struct,order='color')
		self.inter_sorted_nearby=np.sort(nearby_struct,order='color')
		cluster_split,nearby_split,weights_c,weights_n=self.binning(maxima,minima,bins=bins)
		
		if plot==True:
			plt.clf()
			plt.gca().invert_yaxis()
			x_min=minima*np.ones(50)
			x_max=maxima*np.ones(50)
			y=np.linspace(-5,20,50,endpoint=True)
			plt.plot(x_min,y)
			plt.plot(x_max,y)
			plt.scatter(cluster_struct['color'],cluster_struct['Mag'],marker='.',color='blue')
  			plt.scatter(nearby_struct['color'],nearby_struct['Mag'],marker='o',color='black')
  			plt.pause(5)
		
		y_c=[]
		for i in range(bins):
			y_c.append(np.mean(cluster_split[i]['Mag']))
		
		y_n=[]
		for i in range(bins):
			y_n.append(np.mean(nearby_split[i]['Mag']))
		
		y=np.asarray(y_c)-np.asarray(y_n)
		self.y_mean=np.mean(y)
		self.y_std=np.std(y)
		self.y_avg_c=np.average(y,weights=weights_c)
		self.y_avg_n=np.average(y,weights=weights_n)
		return self.y_mean,self.y_std,self.y_avg_c,self.y_avg_n  #self.y_std/(5**0.5)

	def distance_quadraticfit(self):
		nearby_struct,cluster_struct,minima,maxima=self.intersection(self.nearby_struct,self.cluster_struct)
		p1=models.Polynomial1D(2)
		pfit=fitting.LinearLSQFitter()
		
		model_c=pfit(p1,cluster_struct['color'],cluster_struct['Mag'])
		model_n=pfit(p1,nearby_struct['color'],nearby_struct['Mag'])

		self.fit_distance=model_c.c0.value-model_n.c0.value
		return self.fit_distance
	
	def cluster_points(self,X, mu):
		clusters  = {}
		for x in X:
				bestmukey = min([(i[0], np.linalg.norm(x-mu[i[0]])) for i in enumerate(mu)], key=lambda t:t[1])[0]
				try:
						clusters[bestmukey].append(x)
				except KeyError:
						clusters[bestmukey] = [x]
		return clusters
			
	def reevaluate_centers(self,mu, clusters):
		newmu = []
		keys = sorted(clusters.keys())
		for k in keys:
			newmu.append(np.mean(clusters[k], axis = 0))
		return newmu
		
	def has_converged(self,mu,oldmu):
		return (set([tuple(a) for a in mu])== set([tuple(a) for a in oldmu]))
		
	def find_center(self,X,K):
		#Initialise to K random centers
		oldmu = random.sample(X,K)
		mu=random.sample(X,K)
		while not self.has_converged(mu,oldmu):
			oldmu=mu
			clusters =self.cluster_points(X,mu)
			mu = self.reevaluate_centers(oldmu,clusters)
		return mu,clusters
  		
  	def clustering(self,k):
  		nd_cluster=self.cluster_struct.view(dtype=np.float64).reshape(-1,2)
  		mu_c,cluster=self.find_center(nd_cluster,k)
  		nd_nearby=self.nearby_struct.view(dtype=np.float64).reshape(-1,2)
  		mu_n,nearby=self.find_center(nd_nearby,k)
  		cluster_list=[]
  		nearby_list=[]
  		for i in range(k):
  				cluster_list.append(self.listing(cluster[i]))
  				nearby_list.append(self.listing(nearby[i]))
  		#self.intersection()
		#print self.minima,self.maxima,self.nearby_struct,self.cluster_struct
  		return cluster_list,nearby_list
  		
  	def listing(self,list):
  		l=[]
  		for i in range(len(list)):
			l.append((list[i][0],list[i][1]))
		l_struct=np.array(l,dtype=[('color',float),('Mag',float)])
		return l_struct
  	
  	def call(self):
  		self.data()
		self.cluster_struct=self.coordinates(self.color_c,self.V)
		self.nearby_struct=self.coordinates(self.color_n,self.Mv_n)
		
	def update_struct(self,nearby,cluster):
		self.cluster_original=self.cluster_struct
		self.nearby_original=self.nearby_struct
		self.cluster_struct=cluster
		self.nearby_struct=nearby
		self.flag=True
		
	def original_struct(self):
		if self.flag==True:
			self.cluster_struct=self.cluster_original
			self.nearby_struct=self.nearby_original
			self.flag=False
		
if __name__=='__main__':
	#nearby=raw_input('Enter the filename/URL containing nearby stars catalogue : ')
	#cluster=raw_input('Enter the filename/URL containing cluster catalogue : ')
	nearby='NearbyStars.tsv'
	cluster='NGC2420.dat'
	dist=distance(nearby,cluster)
	
	while(True):
		choice=int(raw_input('''Distance to cluster using:
		1.Statistical binning
		2.Quadratic fitting
		3.Clustering and statistical binning
		4.Exit
		'''))
		if choice==1:
			y_mean,y_std,y_avg_c,y_avg_n=dist.distance_bins(bins=5,plot=True)
			print 'Distance modulus using statistical binning is : ',y_mean,y_std,y_avg_c,y_avg_n
			r=10**(y_mean/5+1)
			print 'Distance in parsecs is: ',r
		elif choice==2:
			distance_mod_fit=dist.distance_quadraticfit()
			print 'Distance modulus using quadratic fitting is : ',distance_mod_fit
			r=10**(distance_mod_fit/5+1)
			print 'Distance in parsecs is: ',r
		elif choice==3:
			cluster,nearby=dist.clustering(2)
		
			plt.ion()
			plt.figure()
 			plt.gca().invert_yaxis()
 			plt.scatter(cluster[0]['color'],cluster[0]['Mag'],marker='.',color='blue')
  			plt.scatter(cluster[1]['color'],cluster[1]['Mag'],marker='o',color='black')
 			plt.scatter(nearby[0]['color'],nearby[0]['Mag'],marker='.',color='red')
  			plt.scatter(nearby[1]['color'],nearby[1]['Mag'],marker='o',color='green')
 			plt.pause(1)
 	
 			c=int(raw_input('Enter cluster number: { 0:blue, 1:black }'))
 			n=int(raw_input('Enter nearby number: { 0:red, 1:green }'))
 	
 			plt.clf()
 			plt.gca().invert_yaxis()
 			plt.scatter(cluster[c]['color'],cluster[c]['Mag'],marker='.',color='blue')
  			plt.scatter(nearby[n]['color'],nearby[n]['Mag'],marker='o',color='black')
 			plt.pause(2)
 	
			dist.update_struct(nearby[n],cluster[c])
	
			y_mean,y_std,y_avg_c,y_avg_n=dist.distance_bins(bins=5,plot=True)
			print 'Distance modulus using clustering and binning is : ',y_mean,y_std,y_avg_c,y_avg_n
			r=10**(y_mean/5+1)
			print 'Distance in parsecs is: ',r
			dist.original_struct()
		elif choice==4:
			break
	
	#r=10**(distance/5+1)
	#print 'Distance in parsecs is: ',r