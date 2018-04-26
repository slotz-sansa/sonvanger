# SONVANGER
# 2018/04/20
# Original Author:	SI Lotz
# NWU NASSP Hons
# Student: CH Tleane

import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import sys
import lulu

# CONSTANTS
fitsDir = '../SunspotData_Lotz/MAG_ROTATE/'
#fitsDir = '../SunspotData_Lotz/CONT_ROTATE/'
#fitsDir = '../SunspotData_Lotz/SMART_MASK_MAG/'
imgDir = './img/'
imgExt = '.png'

def loadimg(infile):

	hdulist = pf.open(infile)
	dat = hdulist[0].data

	return dat
	
def test0(dat):
	""" input is KxK fits image matrix
	output is adapted KxK matrix """
	
	R,K = dat.shape
	dK = 50; dR = 50
	
	m = 30
	
	datP = np.zeros_like(dat)
	ipos = np.where(dat>0)[0]
	print(ipos)
	datP[ipos] = dat[ipos]
	
	datN = np.zeros_like(dat)
	ineg = np.where(dat<0)[0]
	datN[ineg] = dat[ineg]
	
	datu = np.zeros_like(dat); datl = datu.copy()
	
	#for i in range(R):
		#datu[i,:] = lulu.U(dat[i,:],m)
		#datl[i,:] = lulu.L(dat[i,:],m)
	
	for i in range(dK,K-dK):
		datu[:,i] = lulu.U(dat[:,i],m)
		datl[:,i] = lulu.L(dat[:,i],m)
	
	plt.figure()
	plt.imshow(datu,cmap='gray')
	plt.title('U(X,m), m = %d' % m); plt.grid()
	
	plt.figure()
	plt.title('L(X,m), m = %d' % m)
	plt.imshow(datl,cmap='gray'); plt.grid()
	
	#plt.figure()
	
	return 999
	
def test1(dat):
	""" input is KxK fits image matrix
	output is adapted KxK matrix """
	
	R,K = dat.shape
	dK = 50; dR = 50
	
	m = 30
	
	datR = np.zeros_like(dat); datK = datR.copy()
	
	for i in range(dR,R-dR):
		datR[i,:] = lulu.L(dat[i,:],m)
	
	for j in range(dK,K-dK):
		datK[:,j] = lulu.L(datR[:,j],m)
	
	plt.figure()
	plt.imshow(dat,cmap='gray')
	#plt.title('U(X,m), m = %d' % m)
	plt.grid()
	
	plt.figure()
	#plt.title('L(X,m), m = %d' % m)
	plt.imshow(datK,cmap='gray'); plt.grid()
	
	#plt.figure()
	
	return 999

def test2(dat=None):
	from skimage import data, io, filters
	
	if dat is None:
		fn = 'fd_M_96m_01d.1285.0000.fits'
		#fn = 'smart_mask_19960709_0004.fits'
		infile = fitdir + fn
		dat = loadimg(infile)
	
	#image = skimage.data.coins()
	# ... or any other NumPy array!
	edges = filters.sobel(np.abs(dat))
	plt.figure()
	plt.imshow(edges,cmap='gray')
	#io.show()
	
	
	return 999

def dnn(i,j,d):
	from matplotlib.cbook import flatten
	
	#i = 5; j = 5; d = 3
	
	zz = np.arange(-d,d+1)
	
	rows=[]
	cols=[]
	
	for z in zz:
		#print(z)
		row = z+i
		#print(row)
		kc = d-np.abs(z)
		cols.append(list(np.arange(j-kc,j+kc+1)))
		
		rows.append(list(np.ones(np.abs(2*kc+1),dtype=int)+row-1))
	
	R = np.array(list(flatten(rows)))
	K = np.array(list(flatten(cols)))
	
	return R,K
	

def main(decSet=None,saveflag=False):
	
	from lulu import connected_region_handler as crh

	#fn = 'fd_Ic_6h_01d.1285.0000.fits'
	fileName = 'fd_M_96m_01d.1285.0000.fits'
	#fn = 'smart_mask_19960709_0004.fits'

	# read FITS image
	infile = fitsDir + fileName
	dat0 = loadimg(infile)

	dat = dat0.copy()

	## plot image using imshow() function from matplotlib
	fig = plt.figure()
	plt.imshow(dat,cmap='gray')
	
	# save image as png file
	if saveflag:
		saveimg(fig,filehdr=fn[:-5])
		plt.close(fig)
		
	# use 2D lulu to decompose image
	## make all nan's zero -- XXX not correct, just a plug for now
	dat[np.where(np.isnan(dat))] = 0
	
	# cast array to int
	#dat = dat.astype(int)
	dat = np.array(dat/10,dtype=int)*10
	# decompose and cast from dict to list
	#decSet = list(lulu.decompose(dat).values())
	if decSet is None:
		decSet = lulu.decompose(dat)
	Nr = len(decSet)
	dsAttr = np.zeros((Nr,2))
	
	i=0
	for a in decSet:
		dsAttr[i,0] = a
		dsAttr[i,1] = len(decSet[a])
		i += 1
		
	plt.figure()
	plt.plot(dsAttr[:,0],dsAttr[:,1],'.')
	plt.xlabel('Area')
	plt.ylabel('# regions')
	
	#plt.figure()
	#plt.imshow(dat)
	
	area0 = int(500); area1 = int(1000)
	
	iii = np.where((dsAttr[:,0]>area0) & (dsAttr[:,0]<area1))[0]
	
	print('area %d - %d: %d' % (area0,area1,len(iii)))
	
	plt.figure()
	plt.pcolor(dat0)

	for i in iii:
		print('%d: %d' % (i,len(decSet[dsAttr[i,0]])))
		c = crh.todense(decSet[dsAttr[i,0]][0])
		y,x = findedge(c)
		#plt.pcolor(c)
		plt.plot(x,y,'.')
		
	return dat0,decSet
	#return dsAttr

def uselulu(dat0):
	"""
	Input:	dat: array from .fits image
	Output:	x
	"""
	
	import lulu
	from lulu import connected_region_handler as crh
	
	# make all nan's zero -- XXX not correct, just a plug for now
	dat0[np.where(np.isnan(dat0))] = 0
	
	# cast array to int
	dat = dat0.astype(int)
	#dat = np.array(dat0,dtype=int)
	
	
	# first decompose magnetogram in to connected regions
	# cast to a list from dict
	decSet = list(lulu.decompose(dat).values())
	
	# print out the number of regions per resolution level
	#for l in range(len(decSet)):
	for l in range(10):
		print('level: %d\t regions: %d' % (l,len(decSet[l])))
		
		# get sum of all images in current level
		for i in range(len(decSet[l])):
			d = crh.todense(decSet[l][i])
			print(d.shape)
	
	# identify some level of interest
	
	
	return dat

def findedge(x):

	"""
	Input is lulu.connectedregionhandler.todense() array with 0/1
	indicating region of interest.
	This just finds the edges of that.
	"""
	
	nr,nk = np.shape(x)
	
	Ire = []; Ice = []
	#Ire1 = []; Ice1 = []
	#Ire2 = []; Ice2 = []
	
	for i in range(nr):
		jj = np.where(np.diff(x[i,:]!=0))[0]
		for j in jj:
			Ire.append(i)
			Ice.append(j)
	
	#for i in range(nk):
		#jj = np.where(np.diff(x[:,i]!=0))[0]
		#for j in jj:
			#Ire2.append(i)
			#Ice2.append(j)
	
	#N1 = len(Ice1); N2 = len(Ice2); N = np.max([nr,nk])
	#Xr = np.zeros(N1)
	
	#for i in range(N1):
		#eks = -int(1+np.log10(N))
		#Xr[i] = a + b*10**(eks)
	
	return Ire,Ice
	
def getbb(x,y):
	
	bb = np.zeros((4,2))
	
	bb[0,:] = [np.min(y),np.max(x)] # ll
	bb[1,:] = [np.min(y),np.min(x)] # ul
	bb[2,:] = [np.max(y),np.min(x)] # ur
	bb[3,:] = [np.max(y),np.max(x)] # lr
	
	return bb

def main2(dat0 = None, fileName = 'fd_M_96m_01d.1285.0000.fits'):
	
	if dat0 is None:
		# read data
		#fn = 'fd_Ic_6h_01d.1285.0000.fits'
		#fileName = 'fd_M_96m_01d.1285.0000.fits'
		#fn = 'smart_mask_19960709_0004.fits'

		## read FITS image, store in array
		infile = fitsDir + fileName
		dat0 = loadimg(infile)
	
	# define area size of interest
	## max area
	areaMax = int(1000)
	## min area
	areaMin = int(300)
	
	# define sensitivity level
	sensLevel = 100
	
	# pre-process image
	dat1 = dat0.copy()
	## make all nan's zero -- XXX not correct, just a plug for now
	dat1[np.where(np.isnan(dat1))] = 0
	## flatten image to desired sensitivity level
	dat2 = np.array(dat1/sensLevel,dtype=int)*sensLevel
	## split image in to pos / neg
	dat2p = dat2.copy()
	dat2p[dat2p<0] = 0 # make all negatives zero
	dat2n = dat2.copy()
	dat2n[dat2n>0] = 0 # make all positives zero
	## make abs value image
	dat2a = np.abs(dat2.copy())
	
	# decompose processed images
	## decompose positive
	dcDat2p = lulu.decompose(dat2p)
	## decompose negative
	dcDat2n = lulu.decompose(dat2n)
	## decompose absolute
	dcDat2a = lulu.decompose(dat2a)
	
	regList, regAttr, bb = manageRegions(dcDat2a,areaMin)
	
	# plot image for each step of processing
	## define colormap to use
	colMapList = ['gray','jet','cool','autumn','terrain','spectral']
	cm = 0 # pick 'gray' out of the list
	
	plt.figure()
	plt.imshow(dat0,cmap=colMapList[cm])
	plt.title(fileName + '\nOriginal')
	plt.grid()
	
	plt.figure()
	plt.imshow(dat1,cmap=colMapList[cm])
	plt.title(fileName + '\nNaN -> 0\'s')
	plt.grid()

	plt.figure()
	plt.imshow(dat2,cmap=colMapList[cm])
	plt.title(fileName + '\nSensitivity p > %d' % sensLevel)
	plt.grid()
	
	plt.figure()
	plt.imshow(dat2p,cmap=colMapList[cm])
	plt.title(fileName + '\nOnly p > 0')
	plt.grid()

	plt.figure()
	plt.imshow(dat2n,cmap=colMapList[cm])
	plt.title(fileName + '\nOnly p < 0')
	plt.grid()

	plt.figure()
	plt.imshow(dat2a,cmap=colMapList[cm])
	plt.title(fileName + '\nAbsolute values')
	plt.grid()
	
	return regList,regAttr,bb
	
def manageRegions(decSet,areaMin,areaMax=None):
	
	from lulu import connected_region_handler as crh
	
	regAttr = np.zeros((len(decSet),2))
	
	i=0
	for s in decSet:
		regAttr[i,0] = s
		regAttr[i,1] = len(decSet[s])
		i += 1

	if areaMax is None:
		iii = np.where((regAttr[:,0]>areaMin))[0]
	else:
		iii = np.where((regAttr[:,0]>areaMin) &\
						(regAttr[:,0]<areaMax))[0]
	
	regList = []
	bbMatrix = np.array([-1,-1,-1,-1])
	for i in iii:
		print('%d: %d' % (i,len(decSet[regAttr[i,0]])))
		for region in decSet[regAttr[i,0]]:
			c = crh.todense(region)
			regList.append(c)
			bbMatrix = np.vstack((bbMatrix,crh.bounding_box(region)))

			
		
	return regList,regAttr,bbMatrix

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
