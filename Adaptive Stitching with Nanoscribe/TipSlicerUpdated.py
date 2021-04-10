#!/usr/bin/env python
import numpy as np
import numpy.random as npRand
import matplotlib.pyplot as plt
import matplotlib.colors as colcol
import matplotlib.cm as cmx
import cv2
import tables
import scipy.ndimage as im
from scipy.spatial import Voronoi
import os
import glob
import re
import scipy.signal as signal

def squareSlice(X,Y,Z,extraParams):
	amp,period,heightOff=extraParams
	ResMat=(amp*np.cos(X*2.0*np.pi/period)+amp*np.cos(Y*2.0*np.pi/period)+amp+heightOff)>Z
	return ResMat

#Z=sin(X)+sin(-X/2+sqrt(3)/2*Y) + sin(-X/2-sqrt(3)/2*Y)
# minimum Z = heightOff
# max Z = heightOff + 4.5*amp
def hexagonalSlice(X,Y,Z,extraParams):
	amp,period,heightOff=extraParams
	ResMat=(amp*np.cos(X*2.0*np.pi/period)+amp*np.cos((-0.5*X+np.sqrt(3.0)/2.0*Y)*2.0*np.pi/period)+amp*np.cos((-X/2.0-np.sqrt(3.0)/2.0*Y)*2.0*np.pi/period)+1.5*amp+heightOff)>Z
	return ResMat
	
def polarSineSlice(X,Y,Z,extraParams):
	xc,yc,period,amp,offset=extraParams
	D=np.sqrt((X-xc)**2+(Y-yc)**2)*2.0*np.pi/period
	ResMat=amp*np.cos(D)+offset>Z
	return ResMat
	
# height falls off as the sum of the distances to the two focal points of an ellipse
# smallAxis determines the scale and is the distance from centre to edge - perpendicular to separating line, 
# maxHeight scales the height
def ellipticalConeSlice(X,Y,Z,extraParams):
	centreX,centreY,sepDist,sepAngle,smallAxis,maxHeight=extraParams
	centre=np.array([centreX,centreY])
	cosang=np.array([np.cos(np.deg2rad(sepAngle)),np.sin(np.deg2rad(sepAngle))])
	f1=centre+sepDist/2.0 * cosang 
	f2=centre-sepDist/2.0 * cosang
	D=np.sqrt((f1[0]-X)**2+(f1[1]-Y)**2)+np.sqrt((f2[0]-X)**2+(f2[1]-Y)**2)
	heightOff=sepDist/np.cos(np.arctan(2*smallAxis/sepDist))
	ResMat=maxHeight*(heightOff-D)/(heightOff-sepDist)>Z
	return ResMat

# will create a X,Y & Z matrix based on the given stlMesh
def getMeshXYZ(stlMesh,hashing,slicing):
	# now determine the extent of the figure, in order to create X,Y and Z axes
	xext=stlMesh[:,:,0].max()-stlMesh[:,:,0].min()
	yext=stlMesh[:,:,1].max()-stlMesh[:,:,1].min()
	zext=stlMesh[:,:,2].max()-stlMesh[:,:,2].min()
	maxd=max(xext,yext)
	maxd=max(maxd,zext)
	
	Xmins=stlMesh[:,:,0].min(axis=1)
	XMI=Xmins.min()
	Xmaxs=stlMesh[:,:,0].max(axis=1)
	XMA=Xmaxs.max()
	Ymins=stlMesh[:,:,1].min(axis=1)
	YMI=Ymins.min()
	Ymaxs=stlMesh[:,:,1].max(axis=1)
	YMA=Ymaxs.max()
	Zmins=stlMesh[:,:,2].min(axis=1)
	ZMI=Zmins.min()
	Zmaxs=stlMesh[:,:,2].max(axis=1)
	ZMA=Zmaxs.max()
	
	offsetXY=5*hashing
	offsetZ=2*slicing
	
	# should set the range to match the hashing resolution! otherwise discrepancy with grid used later
	# assuming the grid goes through zero => each border should be a multiple of het hashing
	nrxmin=np.int(np.round((XMI-offsetXY)/hashing))
	nrxmax=np.int(np.round((XMA+offsetXY)/hashing))+1
	nrymin=np.int(np.round((YMI-offsetXY)/hashing))
	nrymax=np.int(np.round((YMA+offsetXY)/hashing))+1
	
	X,Y=np.meshgrid(np.arange(nrxmin*hashing,nrxmax*hashing,hashing),np.arange(nrymin*hashing,nrymax*hashing,hashing))
	Z=np.arange(ZMI,ZMA+offsetZ,slicing) # should start from zero
	
#	print('Maximum extent: {} with {} slices'.format(maxd,len(Z)))
	if (maxd<1):
		print('Warning, max extent is <1: {}'.format(maxd))
	
	return X,Y,Z
	
def loadStl(stlName,mulDimension=1.,centred=True,zToZero=True,flipZ=False):
	from stl import mesh
	mm = mesh.Mesh.from_file(stlName)
	mm2=mm.vectors*mulDimension
	
	if (centred):
		xCentre=(mm2[:,:,0].min()+mm2[:,:,0].max())/2.0
		yCentre=(mm2[:,:,1].min()+mm2[:,:,1].max())/2.0
		mm2[:,:,0]=mm2[:,:,0]-xCentre
		mm2[:,:,1]=mm2[:,:,1]-yCentre
	if (flipZ):
		mm2[:,:,2]=mm2[:,:,2].max()-mm2[:,:,2]
	if (zToZero):
		mm2[:,:,2]=mm2[:,:,2]-mm2[:,:,2].min()
		
	Zmins=mm2[:,:,2].min(axis=1)
	Zmaxs=mm2[:,:,2].max(axis=1)
	
	return mm2,Zmins,Zmaxs


# complete slicing into a hdf5 file 
# first load the stl with loadStl, this also gives the Zmins and Zmaxs values
# when giving X & Y smaller than total size STL => not possible to get a good infill => use inFillMaster (hdf5-stack) with smaller resolution possibly
def stlToStack(mesh,h5name,X,Y,Z,Zmins,Zmaxs,inFillMaster='',writingMask=np.ones((1,1)),epsilon=0.001,onlyShell=False):
	zslices=len(Z)
	height=X.shape[0]
	width=X.shape[1]
	
	if (writingMask.shape[0]==1):
		writingMask=np.ones(X.shape)
	

	compFilt=tables.Filters(complevel=3)
	with tables.open_file(h5name,filters=compFilt,mode="w") as h5file:
		Array3D=h5file.create_carray("/",'data',tables.BoolAtom(),(zslices,height,width),filters=compFilt)
		ArrayX=h5file.create_carray("/",'X',tables.Float32Atom(),(height,width))
		ArrayY=h5file.create_carray("/",'Y',tables.Float32Atom(),(height,width))
		ArrayZ=h5file.create_carray("/",'Z',tables.Float32Atom(),(zslices,))
		ArrayX[:]=X.copy()
		ArrayY[:]=Y.copy()
		ArrayZ[:]=Z.copy()
		h5file.root.data.attrs.GalvoCentre=[X[0].mean(),Y[:,0].mean(),Z.min()]
		h5file.root.data.attrs.PiezoCentre=[X[0].min(),Y[:,0].min(),Z.min()]
		
		if (inFillMaster!=''):
			with tables.open_file(inFillMaster,mode='r') as infill:
				Xglob=infill.root.X[:]
				Yglob=infill.root.Y[:]
				Zglob=infill.root.Z[:]
				Mglob=infill.root.data
				
				doCarving=False
				dXg=Xglob[0,1]-Xglob[0,0]
				dX=X[0,1]-X[0,0]
				dZg=Zglob[1]-Zglob[0]
				dZ=Z[1]-Z[0]
				
				if ((dX==dXg) and (dZ==dZg)):
					doCarving=True
			
				for i in range(len(Z)):
				
					if (doCarving):
						zind=np.searchsorted(Zglob,Z[i])  # check if it shouldn't be -1
						if (zind>=len(Zglob)):
							zind=len(Zglob)-1
						tempS=carvePiece(X,Y,Xglob,Yglob,Mglob[zind])
						Array3D[i,:,:]=tempS*writingMask
					else:
			#           print(i)
						try:
							tempS=singleSliceStl(mesh,Zmins,Zmaxs,Z[i],X,Y,epsilon,onlyShell=True)
							zind=np.searchsorted(Zglob,Z[i])
							if (zind>=len(Zglob)):
								zind=len(Zglob)-1
							G=Mglob[zind][:]
							Ger=im.binary_erosion(G,iterations=3)   # perhaps keep on eroding until number connected components changes? Or extract skeleton? => should reduce in-fill time?
							nrLabs,Labs=cv2.connectedComponents(Ger.astype(np.uint8))
#							print('First erosion: ',nrLabs)
							oldNrLabs=nrLabs
							if (nrLabs>1):
								while (oldNrLabs==nrLabs):
									Ger=im.binary_erosion(Ger,iterations=3)   
									nrLabs,Labs=cv2.connectedComponents(Ger.astype(np.uint8))  
#									print('Next erosion: ',nrLabs)
							tempS=fillInFromTemplate(tempS,X,Y,Ger,Xglob,Yglob)
							Array3D[i,:,:]=tempS*writingMask
						except:
							print('Problem with slice {} at position {}'.format(i,[X[0].mean(),Y[:,0].mean(),Z.min()]))
							# now just take the slice below the current one.
							if (i>0):
								Array3D[i,:,:]=Array3D[i-1,:,:]
		else:
			for i in range(len(Z)):
	#           print(i)
				#try:
				Array3D[i,:,:]=writingMask*singleSliceStl(mesh,Zmins,Zmaxs,Z[i],X,Y,epsilon,onlyShell=onlyShell)
				#except:
					# now just take the slice below the current one. I could perhaps take the average of above vs below
			#		print('Problem with slice {}'.format(i))
					#Array3D[i,:,:]=singleSliceStl(mm2,mm.normals*mulDimension,Zmins,Zmaxs,Z[i-1],X,Y,epsilon) # Not a perfect solution, but better than nothing!!
	return


# will start with the simplest mask and complexify when needed (=> to function pointer probably)
def applyStackMask(stackName,multMatrix):
	with tables.open_file(stackName,mode='a') as H:
		D=H.root.data
		for i in range(len(D)):
			D[i,:multMatrix.shape[0],:multMatrix.shape[1]]=D[i,:multMatrix.shape[0],:multMatrix.shape[1]]*multMatrix
	return


# give as input already the loaded stl file with the mesh and the normals separately    
def singleSliceStl(mesh,Zmins,Zmaxs,Zslice,X,Y,epsilon=0.001,onlyShell=False,debug=False):
	
	lines=[]
	gI=np.argwhere(np.logical_and((Zmaxs>=Zslice),(Zmins<=Zslice)))[:,0]
	triangs=mesh[gI]

	T=np.zeros((X.shape[0],X.shape[1]),dtype=np.uint8) # this will be the filled in image

	for triang in triangs:
		dists=triang[:,2]-Zslice
		copoints=np.argwhere(np.abs(dists)<epsilon)

		if ((len(copoints)>1)):
#			print('In copoints situation',copoints.shape)
#			if (len(copoints)==2):
			tArray=triang[copoints[:,0]]

			if (len(copoints)!=2):
				tArray=np.row_stack((tArray,triang[copoints[0,0]]))
			#print(tArray)

			lines.append(tArray)
				
			if ((len(tArray)<3)):
				pp=triang[copoints[:,0]]

			else: # for in-plane triangles, add the centre point as a seed point
				vv=tArray[:3,:2]
				yy=convertToPixels(X,Y,vv)
				cv2.fillConvexPoly(T,yy.reshape(-1,1,2),255)
				
		else:
			ps=[]
			if (len(copoints)==1):
				poi=triang[copoints[:,0]]
				ps.append([poi[0,0],poi[0,1],Zslice])	
			
			
			if (dists[0]*dists[1]<0):
#				print('Case 1')
				x1,y1,z1=triang[0,:]
				x2,y2,z2=triang[1,:]
				ts=(Zslice-z1)/(z2-z1)
				xs=(x2-x1)*ts+x1
				ys=(y2-y1)*ts+y1
				ps.append([xs,ys,Zslice])

			if (dists[1]*dists[2]<0):
#				print('Case 2')
				x1,y1,z1=triang[1,:]
				x2,y2,z2=triang[2,:]
				ts=(Zslice-z1)/(z2-z1)
				xs=(x2-x1)*ts+x1
				ys=(y2-y1)*ts+y1
				ps.append([xs,ys,Zslice])

			if (dists[0]*dists[2]<0):
#				print('Case 3')
				x1,y1,z1=triang[0,:]
				x2,y2,z2=triang[2,:]
				ts=(Zslice-z1)/(z2-z1)
				xs=(x2-x1)*ts+x1
				ys=(y2-y1)*ts+y1
				ps.append([xs,ys,Zslice])
				
			if (len(ps)==0):
				print('Freaky!!! Not case 1, 2 or 3')

			if (len(ps)>0):
				pp=np.asarray(ps)
			#	print(pp)
				lines.append(pp)
				
				
		
	li=[(convertToPixels(X,Y,lines[ind][:,:2])).reshape(-1,1,2) for ind in range(len(lines))]
	li3=[ np.round(litem).astype(np.int32) for litem in li]
	res=cv2.polylines(T,li3,False,255,lineType=8)
	
	if (onlyShell):
		if (debug):
			return T,gI,lines
		else:
			return T
	
	if (T.sum()!=0):
		T=colourStlShell(T)
	return T
	
def colourStlShell(TT):
	nrLabs,Labs=cv2.connectedComponents(TT)
	im,cnts,hier = cv2.findContours(TT, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
	hier=hier[0]

	bigTree=[]
	gen=np.argwhere(hier[:,3]==-1)[:,0]

	thisGen,nextGen=checkGeneration(gen,hier,cnts,Labs)
	bigTree.append(thisGen)
	while (len(nextGen)>0):
		thisGen,nextGen=checkGeneration(nextGen,hier,cnts,Labs)
		bigTree.append(thisGen)
		
	N=np.zeros((TT.shape))
	genNr=0
	for gen in bigTree:
		if (genNr%2==0): # should be white
#            print('Should paint white: {}'.format(gen))
			for gg in gen:
				cv2.drawContours(N,cnts,gg,1,-1)
			genNr=genNr+1
		else:
#            print('Should paint black: {}'.format(gen))
			for gg in gen:
				cv2.drawContours(N,cnts,gg,0,-1)
			for gg in gen:
				cv2.drawContours(N,cnts,gg,1,1)
			genNr=genNr+1
	return N

def checkGeneration(gen,hier,cnts,Labs):
	thisGen=gen
	nextGen=[]

	findex=0
	while (findex<len(thisGen)):
		ff=thisGen[findex]
		# find all the clones and all the children of these clones
		clones,childs=checkParent(ff,hier,cnts,Labs)
		#print('output: ',clones,childs)
		for child in childs:
			nextGen.append(child)
		for clone in clones[1:]:
			thisGen=np.append(thisGen,clone)
		findex=findex+1
	return thisGen,np.asarray(nextGen)
	
def checkParent(number,hier,cnts,Labs):
	clones=[]
	trueChildren=[]
	clones.append(number)
	#print('Checking first parent: {}'.format(number))
	fc=hier[number,2]
	if (fc==-1):
		#print('No children')
		pass
	else:
		# first check if we should collapse with some children
		P=np.zeros(Labs.shape,np.uint8)
		cv2.drawContours(P,cnts,number,1,1)
		nrLabsp,Labsp=cv2.connectedComponents(P)  # should only have a single label

		# now get all children
		cc=[fc]
		currPoint=fc
		while (hier[currPoint,0]!=-1):
			currPoint=hier[currPoint,0]
			cc.append(currPoint)
		#print('Children: {}'.format(cc))
		
		for chi in cc:
			chiCont=cnts[chi]
			labelChi=np.median(Labsp[chiCont[:,0,1],chiCont[:,0,0]])
			if (labelChi==1): 
				clones.append(chi)
			else:
				trueChildren.append(chi)
	return clones,trueChildren


# this function will extract a part of Mslice corresponding to the coordinates Xloc,Yloc 
def carvePiece(Xloc,Yloc,Xglob,Yglob,Mslice):
	coord=np.array([[Xloc[0,0],Yloc[0,0]],[Xloc[0,-1],Yloc[-1,0]]])
	ncoord=convertToPixels(Xglob,Yglob,coord)

	carveX=ncoord[:,0]
	startOffX=0
	endOffX=0
	if (carveX[0]<0):
		startOffX=-carveX[0]
		carveX[0]=0
	if (carveX[1]>Xglob.shape[1]-1):
		endOffX=carveX[1]-Xglob.shape[1]+1
		carveX[1]=Xglob.shape[1]-1
	#print(carveX,startOffX,endOffX)


	carveY=ncoord[:,1]
	startOffY=0
	endOffY=0
	if (carveY[0]<0):
		startOffY=-carveY[0]
		carveY[0]=0
	if (carveY[1]>Xglob.shape[0]-1):
		endOffY=carveY[1]-Xglob.shape[0]+1
		carveY[1]=Xglob.shape[0]-1
	#print(carveY,startOffY,endOffY)

	M=np.zeros(Xloc.shape)
	M[startOffY:M.shape[0]-endOffY,startOffX:M.shape[1]-endOffX]=Mslice[carveY[0]:carveY[1]+1,carveX[0]:carveX[1]+1]
	return M

	
def fillInFromTemplate(I,Xloc,Yloc,Template,X,Y):
	Mask=np.zeros((I.shape[0]+2,I.shape[1]+2),dtype=np.uint8)
	Mask[1:-1,1:-1]=I[:]
	TestIm=np.zeros(I.shape,np.uint8)
	gI=np.where(Template>0)
	Xnon=X[gI]
	Ynon=Y[gI]
	ggI=np.where(np.logical_and(Ynon<Yloc.max(),np.logical_and(Ynon>Yloc.min(),np.logical_and(Xnon>Xloc.min(),Xnon<Xloc.max()))))
	if (len(ggI)==0):
		print('No fill in found!')
	seedXY=np.column_stack((Xnon[ggI],Ynon[ggI]))
	seedXYloc=convertToPixels(Xloc,Yloc,seedXY)

	for seed in seedXYloc:
		rr,TestIm,Mask,rect=cv2.floodFill(TestIm,Mask,(seed[0].astype(np.int32),seed[1].astype(np.int32)),1)
	return TestIm
	
def getProjection(stackName,Zmin,Zmax):
	with tables.open_file(stackName,mode='r') as H:
		X=H.root.X[:]
		Y=H.root.Y[:]
		Z=H.root.Z[:]
		D=H.root.data
		bottomI=np.searchsorted(Z,Zmin)
		if (Z[bottomI-1]==Zmin):
			bottomI+=1
		topI=np.searchsorted(Z,Zmax)+1
		P=D[bottomI:topI].max(axis=0)
	return P,X,Y
		

# Options should be approximation used ('NONE','SIMPLE','L1' or 'COS') and number of layers
def getShell(RR,X,Y,z,approx='SIMPLE',random=False):
	ordList=[]
	if (approx=='NONE'):
		conts=cv2.findContours(RR,cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
	elif (approx=='SIMPLE'):
		conts=cv2.findContours(RR,cv2.RETR_LIST,cv2.CHAIN_APPROX_SIMPLE)
	elif (approx=='L1'):
		conts=cv2.findContours(RR,cv2.RETR_LIST,cv2.CHAIN_APPROX_TC89_L1)
	elif (approx=='COS'):
		conts=cv2.findContours(RR,cv2.RETR_LIST,cv2.CHAIN_APPROX_TC89_KCOS)
	Ctour=conts[1]
	for CC in Ctour:
		CC=CC.reshape((-1,2))
		xh=X[CC[:,1],CC[:,0]] 
		yh=Y[CC[:,1],CC[:,0]]
#		xh=np.append(xh,xh[0]) # in order to make the shell fully closed		
#		yh=np.append(yh,yh[0])
		zh=np.ones(len(xh))*z
		fullMat=np.column_stack((xh,yh,zh))
		if (random):
			newMat=np.zeros((fullMat.shape[0]+1,fullMat.shape[1]))
			newStart=np.random.randint(0,len(fullMat))
			newMat[:len(fullMat)-newStart]=fullMat[newStart:]
			newMat[len(fullMat)-newStart:-1]=fullMat[:newStart]
			newMat[-1]=fullMat[newStart]
			ordList.append(newMat)
		else:
			ordList.append(fullMat)
	return ordList

def getHatching(sli,X,Y,z,angle=0,subdiv=1,crossed=True,debug=False):
	Ir=im.rotate(sli,angle,order=0)[::subdiv,:] # only subdivide in the rows after the rotation
	Xr=im.rotate(X,angle)[::subdiv,:]
	Yr=im.rotate(Y,angle)[::subdiv,:]
	if (Ir.sum()==0): # After the subdivide there can be nothing left..
		return [],[]
	
	tesBig=np.zeros([Ir.shape[0],Ir.shape[1]+2],dtype='bool')
	tesBig[:,1:-1]=Ir
	edgeXStart=(tesBig[:,1:]>tesBig[:,:-1])
	edgeXEnd=(tesBig[:,1:]<tesBig[:,:-1])
	SS=np.argwhere(edgeXStart)
	EE=np.argwhere(edgeXEnd)
	EE[:,1]=EE[:,1]-1
	stCoX=Xr[SS[:,0],SS[:,1]]
	stCoY=Yr[SS[:,0],SS[:,1]]
	endCoX=Xr[EE[:,0],EE[:,1]]
	endCoY=Yr[EE[:,0],EE[:,1]]
	ordList=[]
	distanceList=[]                  # keep track of the written length for time estimation
	if (debug):
		return SS,EE,stCoX,stCoY,endCoX,endCoY
	if (not crossed):
		for i in np.arange(len(stCoX)):
			ordList.append(np.array([[stCoX[i],stCoY[i],z],[endCoX[i],endCoY[i],z]]))
			distanceList.append(np.abs(endCoX[i]-stCoX[i])) # normally Y and Z stays constant
	else:
		unSS=np.unique(SS[:,0])
		lineBreaks=np.asarray([np.argwhere(SS[:,0]==unVal)[0][0] for unVal in unSS])
		ordInd=0
		for lb in np.arange(len(lineBreaks)-1):
			indRange=[lineBreaks[lb],lineBreaks[lb+1]]
			if (lb%2==0):
				for ii in np.arange(indRange[0],indRange[1]):
					ordList.append(np.array([[stCoX[ii],stCoY[ii],z],[endCoX[ii],endCoY[ii],z]]))
					distanceList.append(np.sqrt((endCoX[ii]-stCoX[ii])**2+(endCoY[ii]-stCoY[ii])**2))
#					if (np.abs(endCoX[ii]-stCoX[ii])==0):
#					print('Found a zero length :1',endCoX[ii],stCoX[ii],endCoY[ii],stCoY[ii])
			else:
				for ii in np.arange(indRange[1]-1,indRange[0]-1,-1):
					ordList.append(np.array([[endCoX[ii],endCoY[ii],z],[stCoX[ii],stCoY[ii],z]]))
					#distanceList.append(np.abs(endCoX[ii]-stCoX[ii]))
					distanceList.append(np.sqrt((endCoX[ii]-stCoX[ii])**2+(endCoY[ii]-stCoY[ii])**2))
#					if (np.abs(endCoX[ii]-stCoX[ii])==0):
					#print('Found a zero length :2',endCoX[ii],stCoX[ii],endCoY[ii],stCoY[ii])
		# now do last line
		for ii in np.arange(lineBreaks[-1],len(stCoX)):
					ordList.append(np.array([[stCoX[ii],stCoY[ii],z],[endCoX[ii],endCoY[ii],z]]))
					#distanceList.append(np.abs(endCoX[ii]-stCoX[ii]))
					distanceList.append(np.sqrt((endCoX[ii]-stCoX[ii])**2+(endCoY[ii]-stCoY[ii])**2))
#					if (np.abs(endCoX[ii]-stCoX[ii])==0):
					#print('Found a zero length :3',endCoX[ii],stCoX[ii],endCoY[ii],stCoY[ii])
					
	return ordList,np.sum(distanceList)

	
# amShells is the amount of shells. The distance between shells is determined by hatchStep * divider
# scaffoldStep is the distance between the scaffolds (multiplied with divider)
# contApprox is a list of maximum 2 approx methods. First one is for outer shell, rest is done with second one if present
# divider: After the first shell is done with the resolution of the original X,Y. A coarser image is obtain by subsampling
# iterations: is amount of erosions after divider
# atAngles: for every layer, a new 'direction' is chosen in which the main lines are drawn. -1 implies a random direction
def getSliceOfShellScaffold(D,X,Y,Z,sliceIndList,sliceInd,amShells,scaffoldStep,hatchStep=1,atAngles=[0,90,180,270],contApprox=['L1'],divider=1,randomShell=True,doTopBottom=True):
	sliceList=[]
	sliceCodes=[]
	sliceZList=[]
	sliceDist=[]
	slicInd=sliceIndList[sliceInd]
	zslice=Z[slicInd]
	RR=D[slicInd][:]
	if (RR.sum()==0):
		return sliceList,sliceCodes,sliceZList,sliceDist
	
	if (amShells>0):
		startInd=max(0,sliceInd-amShells)
		stopInd=min(len(sliceIndList),sliceInd+amShells+1)
		extInds=sliceIndList[startInd:stopInd]
		Rstack=np.zeros((len(extInds),D.shape[1],D.shape[2]),dtype=bool)
		for jj in np.arange(len(extInds)):
			Rstack[jj]=D[extInds[jj]][:]
		
		extraStart=0
		extend=False
		if (sliceInd<amShells):
			extraStart=-(sliceInd-amShells)
			extend=True
		if (sliceInd+amShells+1>len(sliceIndList)-1):
			extend=True
		if (extend):
			#print(startInd,stopInd,extraStart,extraStart+len(Rstack),len(Rstack))
			Rstack2=np.zeros((2*amShells+1,Rstack.shape[1],Rstack.shape[2]),dtype=bool)
			Rstack2[extraStart:extraStart+len(Rstack)]=Rstack
			Rstack=Rstack2
		
		RR=Rstack[amShells]
		if (RR.sum()==0):
			return sliceList,sliceCodes,sliceZList,sliceDist
		
		Xr=X.copy()
		Yr=Y.copy()
		erosionKernel=disk(hatchStep)
		
		# First do in-plane shells
		for i in range(amShells):
			if (RR.sum()>0):
				if (len(contApprox)==1):
					gSlice=getShell(RR.astype(np.uint8),Xr,Yr,zslice,contApprox[0],random=randomShell)
				else:
					if (i==0):
						gSlice=getShell(RR.astype(np.uint8),Xr,Yr,zslice,contApprox[0],random=randomShell)
					else:
						gSlice=getShell(RR.astype(np.uint8),Xr,Yr,zslice,contApprox[1],random=randomShell)
				sliceList.append(gSlice)
				sliceCodes.append(1) # 1 is for shell
				sliceZList.append(zslice)
				nDist=0
				for gArr in gSlice:
					nDist=nDist+np.sum(np.sqrt(np.sum((gArr[1:]-gArr[:-1])**2,axis=1)))
				sliceDist.append(nDist)
				if ((i==0) and (divider>1)): # should reduce resolution at this stage if necessary
					Rstack=Rstack[:,::divider,::divider]
					RR=Rstack[amShells]
					Xr=X[::divider,::divider]
					Yr=Y[::divider,::divider]
				#RR=im.binary_erosion(RR,iterations=iterations)
				RR=im.binary_erosion(RR,structure=erosionKernel)
		if (RR.sum()==0):
			return sliceList,sliceCodes,sliceZList,sliceDist
		
		# now check if some floors/ceilings need to be printed.
	#   Rstack=Rstack[:,::iterations,::iterations]
	#   RR=RR[::iterations,::iterations]
	#   Xr=Xr[::iterations,::iterations]
	#   Yr=Yr[::iterations,::iterations]
		
		rotAngle=-1
		if (len(np.atleast_1d(atAngles))>1):
			rotAngle=atAngles[sliceInd%len(atAngles)]
		else:
			if (atAngles==-1):
				rotAngle=np.random.randint(0,360)
			else:
				rotAngle=atAngles
		
		
		if (doTopBottom): # if you are doing shell + solid => better to switch of top-bottom.
			inds=np.arange(len(Rstack))
			inds=np.delete(inds,amShells)
			Cmin=Rstack[inds].min(axis=0)
			
			
			FillIn=np.logical_and(RR,np.logical_not(Cmin))
			if (FillIn.sum()>0):
				# now hatch these ceilings/floors without coarsening..
				gHatch,gDist=getHatching(FillIn.astype(np.uint8),Xr,Yr,zslice,angle=rotAngle,subdiv=hatchStep,crossed=True)
				sliceList.append(gHatch)
				sliceCodes.append(3) # 3 is for solid (straight lines with overlap), this is used for bottom/top hatching as well.
				sliceZList.append(zslice)
				sliceDist.append(gDist)
				
			RR=np.logical_and(RR,np.logical_not(FillIn))
	
	else: # no shells, only do scaffolds => Solid mode...
		Xr=X.copy()
		Yr=Y.copy()
		rotAngle=-1
		if (len(np.atleast_1d(atAngles))>1):
			rotAngle=atAngles[slicInd%len(atAngles)]
		else:
			if (atAngles==-1):
				rotAngle=np.random.randint(0,360)
			else:
				rotAngle=atAngles

	#print('Will try to hatch')
	# Finally, hatch in a coarse way the remaining material => Scaffolding
	if (scaffoldStep>0):
		if (RR.sum()>0):
			gHatch,gDist=getHatching(RR.astype(np.uint8),Xr,Yr,zslice,angle=rotAngle,subdiv=scaffoldStep)
			if (len(gHatch)>0):
				sliceList.append(gHatch)
				sliceCodes.append(2) # 2 is for scaffolding (straight lines with no hatching, doubles slicing)
				sliceZList.append(zslice)
				sliceDist.append(gDist)
		
	return sliceList,sliceCodes,sliceZList,sliceDist

#~ def getFullShellScaffold(h5name,amShells,scaffoldJump,atAngles=[0,90,180,270],contApprox=['SIMPLE'],divider=1,iterations=1,randomShell=True,doTopBottom=True):
	#~ h5=tables.open_file(h5name,mode='r')
	
	#~ X=h5.root.X[:]
	#~ Y=h5.root.Y[:]
	#~ Z=h5.root.Z[:]
	#~ D=h5.root.data
	
	#~ allSliceList=[]
	#~ allSliceCodes=[]
	#~ allSliceZList=[]
	#~ allSliceDistList=[]
	#~ for slicInd in range(0,len(Z),1):
#~ #	for slicInd in range(0,65,1):
		#~ sList,sCodes,sZs,newDists=getSliceOfShellScaffold(D,X,Y,Z,slicInd,amShells,scaffoldJump,atAngles,contApprox,divider,iterations,randomShell,doTopBottom)
		#~ if (len(sCodes)==0):
			#~ continue
		#~ for j in range(len(sList)):
			#~ allSliceList.append(sList[j])
			#~ allSliceCodes.append(sCodes[j])
			#~ allSliceZList.append(sZs[j])
#~ #			print(slicInd,j,len(sList),len(newDists),sum(newDists))
			#~ allSliceDistList.append(newDists[j])
	#~ h5.close()
	#~ return allSliceList,allSliceCodes,allSliceZList,allSliceDistList
	
# amShells gives the amount of contours to include
# scaffoldStep gives the amount of lines in the stack to skip while making the solid infill => total dist = hatching of stack * scaffoldStep
# hatchStep is similar to the above, but now for the distance between contours and the distance for the top and bottom 'shells'
# codeSpeeds is typically a list of three speeds for the three line-types: contour, top/bottom shell, scaffold or inner infill
# codeIntens is a list giving the intensities to go with the above speeds
# atAngles: between each slicing layer, the view will be rotated by the next angle to get better uniformity
# sliceIndexList: when empty all slices in the stack are done, otherwise only the indices in this list are used. => Variable slicing distance!
# contApprox: contour approximation technique => see opencv documentation
# divider: reduce image stack size after first contour extraction to speed up execution; not really used anymore
# randomShell: whether to start at a random position along the contour
# doTopBottom: when writing contours or shells, should he also make upper/lower shells? Switch off when in solid mode
# writeHeader: make own gwl header or avoid when writing one of the blocks in a big stitch
# writeColourH5: can also create a colour-coded hdf5-file for visualisation of the instructions
#
# returns the lengths written per coded line-type: combined with the codeSpeeds, this allows for an easy time-estimation for the structure

# a simple Zscaler function giving 1 as an intensity multiplier independent of input height
def noBaseReflection(zheight):
	return 1.0
def stackToGwl(stackName,gwlName,amShells,scaffoldStep,hatchStep,atAngles=[0,90,180,270],codeSpeeds=[10000,10000,10000],codeIntens=[80,80,80],sliceIndexList=[],contApprox=['L1'],divider=1,randomShell=True,doTopBottom=True,writeHeader=True,writeColourH5=False,interfacePos=-1,ZscalerFunc=noBaseReflection):
	with tables.open_file(stackName,mode='r') as h5:
		X=h5.root.X[:]
		Y=h5.root.Y[:]
		Z=h5.root.Z[:]
		D=h5.root.data
		GC=h5.root.data.attrs.GalvoCentre
		
		allSliceList=[]
		allSliceCodes=[]
		allSliceZList=[]
		allSliceDistList=[]
		if (len(sliceIndexList)==0):
			sliceIndexList=range(0,len(Z))
		
		for slicInd in np.arange(len(sliceIndexList)):
			sList,sCodes,sZs,newDists=getSliceOfShellScaffold(D,X,Y,Z,sliceIndexList,slicInd,amShells,scaffoldStep,hatchStep,atAngles,contApprox,divider,randomShell,doTopBottom)
			if (len(sCodes)==0):
				continue
			for j in range(len(sList)):
				allSliceList.append(sList[j])
				allSliceCodes.append(sCodes[j])
				allSliceZList.append(sZs[j])
	#			print(slicInd,j,len(sList),len(newDists),sum(newDists))
				allSliceDistList.append(newDists[j])
				
	writeGwlFile(gwlName,allSliceList,allSliceCodes,allSliceZList,codeSpeeds,codeIntens,stagePosition=GC,writeHeader=writeHeader,interfacePos=interfacePos,debug=True,ZscalerFunc=ZscalerFunc)
	if (writeColourH5):
		Z2=Z[sliceIndexList]
		writeH5ColourGwl(gwlName[:-4]+'_col.h5',allSliceList,allSliceCodes,allSliceZList,X,Y,Z2)
	return getCodeLengths(allSliceDistList,allSliceCodes)

def getCodeLengths(allSliceDistList,allSliceCodes):
	if (len(allSliceCodes)!=0):
		amCodes=max(np.unique(allSliceCodes))
		codeLengths=np.zeros(amCodes)
		try: 
			for i in range(len(allSliceDistList)):
				codeLengths[allSliceCodes[i]-1]=codeLengths[allSliceCodes[i]-1]+allSliceDistList[i]
		except:
			print('Problem with codelength calculation!')
			print(allSliceDistList)
			print(allSliceCodes)
		return codeLengths
	else:
		return [0]



def writeGwlHeader(F,interfacePos,CodeSpeeds,CodeIntensities):
	F.write('GalvoScanMode\r')
	F.write('ContinuousMode\r\r')
#		F.write('ConnectPointsOff\r') # not used for GT in continuous mode
	F.write('InvertZAxis 1\r')  # DiLL configuration
	F.write('% DefocusFactor 0.6\r') # for air objective
	F.write('TiltCorrectionOff\r') # not available with GT
	F.write('StageVelocity 200\r')
	F.write('PowerScaling 1\r')
	F.write('GalvoAcceleration 2\r\r')

	F.write('var $folderNumber = 1\r\r')
	
	F.write('% 1 is for shell, 2 for interior and 3 for top/bottom\r')
	for j in range(len(CodeSpeeds)):
		F.write('var $Speed{} = {}\r'.format(j+1,CodeSpeeds[j]))
		F.write('var $Power{} = {}\r'.format(j+1,CodeIntensities[j]))
#            F.write('local $block = 0\r')
	
				
	F.write('\r')
	F.write('DebugModeOn\r')
	F.write('TimeStampOn\r')
		
	F.write ('var $perSliceDebug = 0\r')
	F.write ('var $perBlockDebug = 1\r')
	F.write ('var $writePics = 1\r\r')
	
	F.write('XOffset 0\r')  # Will use offsets later so initialise here
	F.write('YOffset 0\r')
	F.write('ZOffset 0\r\r') # Can add 3.5 here for a proper interface positioning
	
	F.write('var $doTiltCorrection = 0 % When >0, a correction for the interface location will be used instead of FindInterfaceAt\r')
	F.write('var $tiltdX=0  % The coefficient a in the fitted plane formula: Z=ax+by+c\r')
	F.write('var $tiltdY=0  % The coefficent b in the above formula\r')
	F.write('var $tiltInterface=3  % The equivalent FindInterfaceAt location\r')
	F.write('var $tiltMargin=20 % This is the amount of +- correction that will be allowed\r\r')
	F.write('if $doTiltCorrection>0 \r')
	F.write('	local $tiltInt=$tiltMargin+$tiltInterface\r')
	F.write('	FindInterfaceAt $tiltInt\r')
	F.write('	ZDrivePosition\r')
	F.write('	AddZOffset $tiltMargin\r')
	F.write('end\r\r')

	
	
	F.write('var $interfacePos = {} % if >=0 will be done for each block => clashes with TiltCorrection!!\r'.format(interfacePos))
	F.write('var $piezoInterf = 10  % piezo will be moved here to do the FindInterface\r')
	F.write('var $piezoWrite = 150 % afterwards it will come back here to write the block\r\r')
	
	F.write('var $currX = 0\r') # This will be used for relative positioning of different stitching blocks or arrays
	F.write('var $currY = 0\r')
	F.write('var $currZ = 0\r\r')
	
	F.write('% Write out the static parameters to the log file\r')
	F.write('ShowParameter\r')
	F.write('MessageOut "Speed1 %d Speed2 %d Speed3 %d" #($Speed1,$Speed2,$Speed3)\r')
	F.write('MessageOut "Power1 %d Power2 %d Power3 %d" #($Power1,$Power2,$Power3)\r')
	F.write('MessageOut "tiltMargin %d tiltInterface %d interfacePos %d" #($tiltMargin,$tiltInterface,$interfacePos)\r\r')
	F.write('% -----------   End of header    -----------\r\r')
	return
	


# Writes to GWL file
# The adaptation of intensity on Z-position will be handled here. 
# stagePosition is used for stitching multiple blocks. Will move stage to this position.
# In the header the current location is initialised to (0,0,0) and will then move stage to given location.
# X,Y and ZOffsets are next incremented to this stagePosition, so as to centre the coordinate system for Galvowriting.
# so make sure that all coordinates are reachable when the coordinates are centered on this stagePosition
# interfacePos <0 implies that FindInterfaceAt is not called for this block, if >=0 it is called with the given value
def writeGwlFile(filename,allSliceList,allSliceCodes,allSliceZList,CodeSpeeds=[10000,10000,10000],CodeIntensities=[80,80,80],stagePosition=[0,0,0],interfacePos=-1,ZscalerFunc=noBaseReflection,writeHeader=True,debug=True):
	with open(filename,mode="w") as F:
		if (writeHeader):
			writeGwlHeader(F,interfacePos,CodeSpeeds,CodeIntensities)
		
		
		# first handle relative movement of stage if necessary
		# also adapt X,Y and ZOffset parameters to have properly centered Galvo-coordinates for Nanoscribe
		# even if inside the gwl blocks, the original coordinate system is maintained.
		F.write('NewStructure\r')
		F.write('% The following commands are used for relative positioning.\r')
		F.write('% If you change them for manual positioning, do not forget to change (remove) the AddX(Y or Z)Offset values as well\r')
		F.write('local $locX = {:.0f}\r'.format(stagePosition[0])) # target values
		F.write('local $locY = {:.0f}\r'.format(stagePosition[1]))
		F.write('local $locZ = {:.0f}\r\r'.format(stagePosition[2]))
		
		F.write('local $xmov = $locX-$currX\r')
		F.write('MoveStageX $xmov\r')
		F.write('AddXOffset -$xmov\r')
		F.write('local $ymov = $locY-$currY\r')
		F.write('MoveStageY $ymov\r')
		F.write('AddYOffset -$ymov\r\r')
		
		F.write('local $zmov = $locZ-$currZ\r')
#		F.write('AddZDrivePosition $zmov\r')
		F.write('AddZOffset -$zmov\r\r')
		
		F.write('% Now find the interface if $interfacePos >=0\r')
		F.write('if $interfacePos >= 0 \r')
		F.write('	PiezoGotoY $piezoInterf \r')
		F.write('	FindInterfaceAt $interfacePos \r')
		F.write('	PiezoGotoY $piezoWrite \r')
		F.write('	ZDrivePosition\r')
		F.write('	set $currZ=0  % after a findinterface I am at zero again.\r')
		F.write('end\r\r')
		
		F.write('local $zmov = $locZ-$currZ\r')
		F.write('AddZDrivePosition $zmov\r')
#		F.write('AddZOffset -$zmov\r\r')
		
		F.write('set $currX=$locX\r')
		F.write('set $currY=$locY\r')
		F.write('set $currZ=$locZ\r\r')
		
		
		
		F.write('% Now handle the interface tilt correction => clashes with FindInterfaceAt\r')
		F.write('local $tiltCorr = 0\r')
		F.write('if $doTiltCorrection >0 \r')
		F.write('	set $tiltCorr = $tiltdX*$currX+$tiltdY*$currY\r')
		F.write('	AddZOffset $tiltCorr\r')
		F.write('end\r\r')
		
		F.write('% Now write out all the relevant numbers to the log file \r')
		F.write('MessageOut "tiltCorr %_2f " #($tiltCorr)\r')
		F.write('MessageOut "locX %_2f locY %_2f locZ %_2f" #($locX,$locY,$locZ)\r')
		
		
		Zscaler=1.0
		pointCounter=0
		lineCounter=0
		oldCode=-1
		oldZ=-1
		zCount=len(np.unique(allSliceZList))
		zIndex=0 # a counter to allow following the slice progress
		writeSliceInfo=False
		for i in np.arange(len(allSliceList)):
			sliceArray=allSliceList[i]
			zz=allSliceZList[i]
			cod=allSliceCodes[i] # can modify speed etc based on this code
			if ((cod!=oldCode) or (zz!=oldZ)):
				if (debug and (zz!=oldZ)):
					# first write old slice info
					F.write('if $perSliceDebug>0 \r')
					F.write('	MessageOut "Finished slice {0:05d} out of {1:d} for block %05d" #($folderNumber)\r'.format(zIndex,zCount))
					F.write('	if $writePics>0 \r')
					F.write('		Wait 0.1\r')
					filenam='PicPerSlice/Pic_Block_%05d_{0:05d}.tif'.format(zIndex) # I could also add the timestamp in the name or location (in case of stitching)
					F.write('		CapturePhoto "{}"  #($folderNumber)\r'.format(filenam))
					F.write('	end\r')
						
					
					zIndex=zIndex+1
					F.write('	MessageOut "Started slice {0:05d} out of {1:d} for block %05d" #($folderNumber)\r'.format(zIndex,zCount))
					F.write('end\r\r')
					
				F.write('ScanSpeed $Speed{}\r'.format(cod))
				F.write('LaserPower $Power{}\r'.format(cod))
				Zscaler=ZscalerFunc(zz)
				F.write('MultLaserPower {0:.2f}\r\r'.format(Zscaler))
				
					
				oldCode=cod
				oldZ=zz
				
			for j in np.arange(len(sliceArray)):
				for k in np.arange(len(sliceArray[j])):
					if (sliceArray[j].shape[1]==3):
						F.write('{0:.2f} {1:.2f} {2:.2f}\r'.format(sliceArray[j][k,0],sliceArray[j][k,1],sliceArray[j][k,2]))
					elif (sliceArray[j].shape[1]==4):
						F.write('{0:.2f} {1:.2f} {2:.2f} {3:.2f}\r'.format(sliceArray[j][k,0],sliceArray[j][k,1],sliceArray[j][k,2],sliceArray[j][k,3]))
					pointCounter=pointCounter+1
				F.write("write\r\r")
				lineCounter=lineCounter+1
		
		if (debug):
			F.write('if $perSliceDebug>0 \r')
			F.write('MessageOut "Finished slice {0:05d} out of {1:d} for block %05d" #($folderNumber)\r'.format(zIndex,zCount))
			F.write('	if $writePics>0 \r')
			F.write('		Wait 0.1\r')
			filenam='PicPerSlice/Pic_Block_%05d_{0:05d}.tif'.format(zIndex) # I could also add the timestamp in the name or location (in case of stitching)
			F.write('		CapturePhoto "{}"  #($folderNumber)\r'.format(filenam))
			F.write('	end\r')
			F.write('end\r\r')
		
		if (writeHeader): # if you don't want a header, than it is likely that you are not writing a single block => save the log only once!
			messfile='MyLogfile Block%d' # should make this name more unique in case of stitching but cannot make new folder
			if (debug):
				F.write('SaveMessages "{}" #($folderNumber)\r'.format(messfile))
				
		F.write('% Now remove the interface tilt correction again, if applied\r')
		F.write('if $doTiltCorrection >0 \r')
		F.write('	AddZOffset -$tiltCorr\r')
		F.write('end\r\r')

	return pointCounter,lineCounter
	
# Writes an empty GWL file, used for interface tilt determination
def writeEmptyGwlFile(filename,stagePosition=[0,0,0],interfacePos=0,writeHeader=True,debug=True):
	with open(filename,mode="w") as F:
		if (writeHeader):
			writeGwlHeader(F,interfacePos,[0],[0])
			
			
		# first handle relative movement of stage if necessary
		# also adapt X,Y and ZOffset parameters to have properly centered Galvo-coordinates for Nanoscribe
		# even if inside the gwl blocks, the original coordinate system is maintained.
		F.write('% The following commands are used for relative positioning.\r')
		F.write('% If you change them for manual positioning, do not forget to change (remove) the AddX(Y or Z)Offset values as well\r')
		F.write('local $locX = {:.0f}\r'.format(stagePosition[0])) # target values
		F.write('local $locY = {:.0f}\r'.format(stagePosition[1]))
		F.write('local $locZ = {:.0f}\r\r'.format(stagePosition[2]))
		
		F.write('local $xmov = $locX-$currX\r')
		F.write('MoveStageX $xmov\r')
		F.write('AddXOffset -$xmov\r')
		F.write('local $ymov = $locY-$currY\r')
		F.write('MoveStageY $ymov\r')
		F.write('AddYOffset -$ymov\r')
		F.write('local $zmov = $locZ-$currZ\r')
		F.write('AddZDrivePosition $zmov\r')
		F.write('AddZOffset -$zmov\r')
		
		F.write('set $currX=$locX\r')
		F.write('set $currY=$locY\r')
		F.write('set $currZ=$locZ\r\r')
		
		# handle using FindInterfaceAt
		F.write('if $interfacePos >= 0 \r')
		F.write('	PiezoGotoY $piezoInterf \r')
		F.write('	FindInterfaceAt $interfacePos \r')
		F.write('	PiezoGotoY $piezoWrite \r')
		F.write('	ZDrivePosition\r')
		F.write('end\r\r')
		
		
	return

def writeGwlForRegex(outGwlName,filePattern,codeSpeeds=[10000,10000,10000],codeIntensities=[80,80,80],interfacePos=-1,debug=True):
	fp=os.path.split(filePattern)
	if (fp[0]==''):
		fnames,fnumbs=findFiles(filePattern)
	else:
		fnames,fnumbs=findFiles(fp[1],fp[0])
	with open(outGwlName,mode="w") as F:
		writeGwlHeader(F,interfacePos,codeSpeeds,codeIntensities)
		
		blockCounter=0
		for ff in fnames:
			F.write('if $perBlockDebug>0 \r')
			F.write('	MessageOut "Starting block {0:d} out of {1:d}"\r'.format(blockCounter,len(fnames)))
			F.write('end \r')
			#F.write('include {}\r'.format(os.path.split(ff)[1])) # perhaps use relpath first and replace / with \
			tempPath=os.path.relpath(ff)
			startP=tempPath.find('/')
			F.write('include {}\r'.format(tempPath[startP+1:].replace('/','\\')))
			
			F.write('if $perBlockDebug>0 \r')
			F.write('	MessageOut "Finished block {0:d} out of {1:d}"\r'.format(blockCounter,len(fnames)))
			F.write('	if $writePics>0 \r')
			filenam='PicPerBlock/Block_%05d.tif'
			F.write('		Wait 1\r')
			F.write('		CapturePhoto "{}"  #($folderNumber)\r'.format(filenam))
			F.write('	end\r')
			
			F.write('end \r\r\r')
				
			F.write('set $folderNumber=$folderNumber+1\r')
			blockCounter=blockCounter+1
			
		if (debug):
			F.write('SaveMessages "MyLogfile"\r')
	
	return






# other good colourmaps 'viridis', 'magma','plasma','inferno', 'coolwarm', 'rainbow','spectral'
def writePointCloud(allSliceList,allSliceCodes,filename,divider=10,colorName='coolwarm',goodCodes=[0,1,2,3,4]):
	#just modified the input parameters
	pointCounter=0
	pointList=np.zeros((0,3))
	pointList=[]
	for i in np.arange(len(allSliceList)):
		cod=allSliceCodes[i] # can modify speed etc based on this code
		if (cod not in goodCodes):
			continue
		sliceArray=allSliceList[i]
		for contour in sliceArray:
			for pp in contour:
				if (np.random.randint(divider)==0):
#                pointList=np.append(pointList,contour,axis=0)
					pointList.append(pp)
	
	Pts=np.asarray(pointList)
			
	Hs=np.unique(np.asarray(Pts[:,2]))

	amPoints=len(Pts)
	F=open(filename,mode='w')
	F.write('# .PCD v.7 - Point Cloud Data file format\n')
	F.write('VERSION .7\n')
	F.write('FIELDS x y z rgb\n')
	F.write('SIZE 4 4 4 4\n')
	F.write('TYPE F F F I\n')
	F.write('COUNT 1 1 1 1\n')
	F.write('WIDTH {}\n'.format(amPoints))
	F.write('HEIGHT 1\n')
	F.write('VIEWPOINT 0 0 0 1 0 0 0\n')
	F.write('POINTS {}\n'.format(amPoints))
	F.write('DATA ascii\n')
	
	mappie = plt.get_cmap(colorName) 
	cNorm=colcol.Normalize(vmin=Hs.min(),vmax=Hs.max())
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=mappie)
	cc=scalarMap.to_rgba(Hs)
	rgb = (((cc[:,0]*255).astype(int) << 16) + ((cc[:,1]*255).astype(int) << 8) + (cc[:,2]*255).astype(int))
	
	
	for pp in Pts:
		gi=np.searchsorted(Hs,pp[2]) # only 1 height per slice
		col=rgb[gi]
		F.write('{:.4f} {:.4f} {:.4f} {:}\n'.format(pp[0],pp[1],pp[2],col))
			
	F.close()
	return amPoints

# other good colourmaps 'viridis','magma','plasma','inferno','rainbow','spectral'
def writePointCloud2(allSliceList,allSliceCodes,allSliceZList,filename,divider=10,colorName='coolwarm'):
	zmin=np.min(allSliceZList)
	zmax=np.max(allSliceZList)


	F=open(filename,mode='w')
	F.write('# .PCD v.7 - Point Cloud Data file format\n')
	F.write('VERSION .7\n')
	F.write('FIELDS x y z rgb\n')
	F.write('SIZE 4 4 4 4\n')
	F.write('TYPE F F F I\n')
	F.write('COUNT 1 1 1 1\n')
	F.write('WIDTH ')
	writeLoc1=F.tell()
	filler='                    ' # this should be larger than amount of points
	F.write(filler+'\n')
	F.write('HEIGHT 1\n')
	F.write('VIEWPOINT 0 0 0 1 0 0 0\n')
	F.write('POINTS ')
	writeLoc2=F.tell()
	F.write(filler+'\n')
	F.write('DATA ascii\n')

	Hs=np.arange(zmin,zmax,0.1)
	mappie = plt.get_cmap(colorName) 
	cNorm=colcol.Normalize(vmin=Hs.min(),vmax=Hs.max())
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=mappie)
	cc=scalarMap.to_rgba(Hs)
	rgb = (((cc[:,0]*255).astype(int) << 16) + ((cc[:,1]*255).astype(int) << 8) + (cc[:,2]*255).astype(int))    
	
	pointCounter=0
	for i in np.arange(len(allSliceList)):
		cod=allSliceCodes[i] # can modify speed etc based on this code
		if (cod!=2):
			continue
		sliceArray=allSliceList[i]
		gi=np.searchsorted(Hs,sliceArray[0][0,2]) # only 1 height per slice
		col=rgb[gi]
		for contour in sliceArray:
			for pp in contour:
				if (np.random.randint(divider)==0):
					F.write('{:.4f} {:.4f} {:.4f} {:}\n'.format(pp[0],pp[1],pp[2],col))
					pointCounter=pointCounter+1

	F.close()
	
	amPoints='{}'.format(pointCounter)
	neededSpaces=len(filler)-len(amPoints)
	F=open(filename,mode='r+')
	
	F.seek(writeLoc1)
	F.write(amPoints)
	for i in np.arange(neededSpaces):
		F.write(' ')
	F.write('\n')
	
	F.seek(writeLoc2)
	F.write(amPoints)
	for i in np.arange(neededSpaces):
		F.write(' ')
	F.write('\n')
	F.close()
	
	return pointCounter

# X, Y & Z are ranges with their values. Will make a 2d matrix of X & Y myself	
def writeEmptyStack(h5name,X,Y,Z,compression=True):
	zslices=len(Z)
	Xm,Ym=np.meshgrid(X,Y)
	height=Xm.shape[0]
	width=Xm.shape[1]

	compFilt=None
	if compression:
		compFilt=tables.Filters(complevel=3)

	with tables.open_file(h5name,filters=compFilt,mode="w") as h5file:
		#Array3D=h5file.createCArray("/",'data',tables.UInt8Atom(),(amount,height,width),filters=compFilt)
		Array3D=h5file.create_carray("/",'data',tables.BoolAtom(),(zslices,height,width),filters=compFilt)
		ArrayX=h5file.create_carray("/",'X',tables.Float32Atom(),(height,width))
		ArrayY=h5file.create_carray("/",'Y',tables.Float32Atom(),(height,width))
		ArrayZ=h5file.create_carray("/",'Z',tables.Float32Atom(),(zslices,))
		ArrayX[:]=Xm.copy()
		ArrayY[:]=Ym.copy()
		ArrayZ[:]=Z.copy()

	return
	
def formulaToStack(h5name,X,Y,Z,sliceFormula,extraParams=[],writingMask=np.ones((1,1)),compression=True):
	zslices=len(Z)
	height=X.shape[0]
	width=X.shape[1]

	if (writingMask.shape[0]==1):
		writingMask=np.ones(X.shape)

	compFilt=None
	if compression:
		compFilt=tables.Filters(complevel=3)

	with tables.open_file(h5name,filters=compFilt,mode="w") as h5file:
		#Array3D=h5file.createCArray("/",'data',tables.UInt8Atom(),(amount,height,width),filters=compFilt)
		Array3D=h5file.create_carray("/",'data',tables.BoolAtom(),(zslices,height,width),filters=compFilt)
		ArrayX=h5file.create_carray("/",'X',tables.Float32Atom(),(height,width))
		ArrayY=h5file.create_carray("/",'Y',tables.Float32Atom(),(height,width))
		ArrayZ=h5file.create_carray("/",'Z',tables.Float32Atom(),(zslices,))
		ArrayX[:]=X.copy()
		ArrayY[:]=Y.copy()
		ArrayZ[:]=Z.copy()
		for i in range(len(Z)):
			Array3D[i,:,:]=writingMask*sliceFormula(X,Y,Z[i],extraParams)
		# Still add GalvoCentre and PiezoCentre attributes
		h5file.root.data.attrs.GalvoCentre=[X[0].mean(),Y[:,0].mean(),Z.min()]
		h5file.root.data.attrs.PiezoCentre=[X[0].min(),Y[:,0].min(),Z.min()]
	return
	
def formulaToStack2(h5name,X,Y,Z,sliceFormula,maskFormula,extraParams,compression=True):
	zslices=len(Z)
	height=X.shape[0]
	width=X.shape[1]

	compFilt=None
	if compression:
		compFilt=tables.Filters(complevel=3)

	with tables.open_file(h5name,filters=compFilt,mode="w") as h5file:
		#Array3D=h5file.createCArray("/",'data',tables.UInt8Atom(),(amount,height,width),filters=compFilt)
		Array3D=h5file.create_carray("/",'data',tables.BoolAtom(),(zslices,height,width),filters=compFilt)
		ArrayX=h5file.create_carray("/",'X',tables.Float32Atom(),(height,width))
		ArrayY=h5file.create_carray("/",'Y',tables.Float32Atom(),(height,width))
		ArrayZ=h5file.create_carray("/",'Z',tables.Float32Atom(),(zslices,))
		ArrayX[:]=X.copy()
		ArrayY[:]=Y.copy()
		ArrayZ[:]=Z.copy()
		for i in range(len(Z)):
			mask=maskFormula(X,Y,Z[i],extraParams)
			Array3D[i,:,:]=mask*sliceFormula(X,Y,Z[i],extraParams)
		# Still add GalvoCentre and PiezoCentre attributes
		h5file.root.data.attrs.GalvoCentre=[X[0].mean(),Y[:,0].mean(),Z.min()]
		h5file.root.data.attrs.PiezoCentre=[X[0].min(),Y[:,0].min(),Z.min()]
	return
	
def writeH5ColourGwl(h5name,allSliceList,allSliceCodes,allSliceZList,X,Y,Z,compression=True,zprecision=0.05):
	zslices=len(Z)
	height,width=X.shape
	hashing=X[0,1]-X[0,0]
	compFilt=None
	if compression:
		compFilt=tables.Filters(complevel=3)
	
	with tables.open_file(h5name,filters=compFilt,mode="w") as h5file:
		
		Array3D=h5file.create_carray("/",'data',tables.UInt8Atom(),(zslices,height,width*3),filters=compFilt)
		ArrayX=h5file.create_carray("/",'X',tables.Float32Atom(),(height,width))
		ArrayY=h5file.create_carray("/",'Y',tables.Float32Atom(),(height,width))
		ArrayZ=h5file.create_carray("/",'Z',tables.Float32Atom(),(zslices,))
		ArrayX[:]=X.copy()
		ArrayY[:]=Y.copy()
		ArrayZ[:]=Z.copy()
		h5file.root.data.set_attr('PIXFORMAT',b'RGB8') # so it can be visualised by gigaviewer

		currZ=-1
		T=np.zeros((X.shape[0],X.shape[1],3),dtype=np.uint8)
		pos=0
		colours=[(255,0,0),(0,255,0),(0,0,255),(120,120,120)]
		for i in range(len(allSliceZList)):
			if (np.abs(allSliceZList[i]-currZ)>zprecision): # first time at this height => do initialisation
				pos=np.searchsorted(Z,allSliceZList[i]-zprecision)
				#print(i,allSliceZList[i],currZ,pos)
				Array3D[pos,:,:]=T.reshape((T.shape[0],T.shape[1]*3)) # now write old T (for that one is finished now), pos 0 will be overwritten but ok
				T=np.zeros((X.shape[0],X.shape[1],3),dtype=np.uint8) # make new T ready for next layers
				currZ=allSliceZList[i]
			
			# NEED to accumulate all images on a given height (add colours not replace)!
						
			profs=allSliceList[i]
			li=[(convertToPixels(X,Y,profs[ind][:,:2])).reshape(-1,1,2) for ind in range(len(profs))]
			#li=[profs[ind][:,:2].reshape(-1,1,2) for ind in range(len(profs))]
			
			#li2=(li-X.min())/hashing # Only correct when X.min() = Y.min() => should substract differently for X and Y
			li3=[ np.round(litem).astype(np.int32) for litem in li]
	#        li=[((profs[ind][:,:2]-X.min())/(X.max()-X.min())*X.shape[0]).astype(np.int32).reshape(-1,1,2) for ind in range(len(profs))]
			res=cv2.polylines(T,li3,False,colours[allSliceCodes[i]-1])
		Array3D[pos,:,:]=T.reshape((T.shape[0],T.shape[1]*3)) # don't forget last slice

def disk(radius, dtype=np.uint8):
	"""Generates a flat, disk-shaped structuring element.
	A pixel is within the neighborhood if the euclidean distance between
	it and the origin is no greater than radius.
	Parameters
	----------
	radius : int
		The radius of the disk-shaped structuring element.
	Other Parameters
	----------------
	dtype : data-type
		The data type of the structuring element.
	Returns
	-------
	selem : ndarray
		The structuring element where elements of the neighborhood
		are 1 and 0 otherwise.
	"""
	L = np.arange(-radius, radius + 1)
	X, Y = np.meshgrid(L, L)
	return np.array((X ** 2 + Y ** 2) <= radius ** 2, dtype=dtype)


### Typical regex would be "drops0ts(\d+).dat", placing the brackets will make the sort based on captured integer
def findFiles(regex,dirName=os.curdir,isFloat=False):
	filenames=glob.glob(os.path.join(os.path.abspath(dirName),"*"))
	if len(filenames) == 0:
		print("No pictures found")
		return (0,None)
	# two steps, first get a list of those filenames matching the given regex
	newList=[]
	numList=[]
	for gfile in filenames:
		mat=re.search(regex,gfile)
		if (mat):
			newList.append(gfile)
			if (len(mat.groups())>0):
				if (isFloat):
					nr=float(mat.groups()[-1])
				else:
					nr=int(mat.groups()[-1])
				numList.append(nr)
			
	# now sort the files according to their number. Always fitting regex \d+.dat to extract number for each file
	if (len(numList)>0):
		sI=np.argsort(numList)
		sortF= [newList[i] for i in sI]
	else:
		newList.sort()
		sortF=newList
	
	return sortF,np.sort(numList)

def getSquareStitchCentres(T,writeRad,X,Y):
	dX=X[0,1]-X[0,0]
	dist_pix=int(round(2*writeRad/dX))
	xx=np.arange(X.min(),X.max(),2*writeRad)
	yy=np.arange(Y.min(),Y.max(),2*writeRad)
	NX,NY=np.meshgrid(xx,yy)
	pp=np.column_stack((NX.ravel(),NY.ravel()))
	#NC_pix=FormulaSlicer.convertToPixels(X,Y,NC)
	
	#pp=calcHexCenters(amY,amX,writeRad)
	# now only keep those centres which have some work to do
	regions, vertices = makeVoronoiPolygons(pp)
	goodRegions=[]
	goodPoints=[]
	for i in np.arange(len(pp)):
		polygon = vertices[regions[i]]
		MM=makeVoronoiMask(polygon,X,Y)
		if ((MM*T).sum()!=0):
			goodPoints.append(pp[i])
			goodRegions.append(regions[i])

	hexPoints=np.asarray(goodPoints)
	return hexPoints
	
def getHexStitchCentres(T,writeRad,X,Y):
	amX=int(round((X.max()-X.min())*1.2/(2*writeRad)))
	amY=int(round((Y.max()-Y.min())*1.4/(2*writeRad)))
	pp=calcHexCenters(amY,amX,writeRad)
	# now only keep those centres which have some work to do
	regions, vertices = makeVoronoiPolygons(pp)
	goodRegions=[]
	goodPoints=[]
	for i in np.arange(len(pp)):
		polygon = vertices[regions[i]]
		MM=makeVoronoiMask(polygon,X,Y)
		if ((MM*T).sum()!=0):
			goodPoints.append(pp[i])
			goodRegions.append(regions[i])

	hexPoints=np.asarray(goodPoints)
	return hexPoints

# this will return the coords of the circle centers in a tesselation
def calcHexCenters(nrows,ncols,R):
	offx=R*np.sqrt(3)
	offy=1.5*R
	offExtra=np.sqrt(3)*R/2.0
	
	coordList=[]
	
	# draw circumscribing circles
	for i in np.arange(ncols):
		for j in np.arange(nrows):
			if np.logical_and(j%2==1,i==0):
				cx=i*offx-j%2*offExtra
				cy=j*offy
				coordList.append((cx,cy))
			cx=i*offx+j%2*offExtra
			cy=j*offy
			coordList.append((cx,cy))
	
	coords=np.asarray(coordList)
	centerx=(coords[:,0].max()+coords[:,0].min())/2.0
	centery=(coords[:,1].max()+coords[:,1].min())/2.0
	centCoord=coords-[centerx,centery]
	
	return centCoord

def getOptimalStitchCentres(T,VoteRadius,WriteRadius,X,Y):
	dPix=X[0,1]-X[0,0]
	VoteRad_pix=int(round(VoteRadius/dPix))
	WriteRad_pix=int(round(WriteRadius/dPix))
	
	In=T.copy()
	centList=[]
	while(In.sum()!=0):
		Vote=getVoteFullCircle(In,VoteRad_pix)#,oldT,oldVote)
		cp,Out=getVoteCentre(In,Vote,WriteRad_pix)
		diff=(In-Out).sum()
		if (diff==0): # something went wrong with differential voting => go back to absolute voting
			print('Something went wrong with voting')
		In=Out.copy()
		centList.append(cp)
	Cents=np.asarray(centList)
	Xc=X[Cents[:,1],Cents[:,0]]
	Yc=Y[Cents[:,1],Cents[:,0]]
	points=np.column_stack((Xc,Yc)).astype(np.float64)
	return points
	
def getOptimalStitchCentres_Bulk(T,VoteRadius,WriteRadius,X,Y):
	dPix=X[0,1]-X[0,0]
	VoteRad_pix=int(round(VoteRadius/dPix))
	WriteRad_pix=int(round(WriteRadius/dPix))
	
	In=T.copy()
	centList=[]
	while(In.sum()!=0):
		Vote=getVoteFullCircle_Bulk(In,VoteRad_pix)#,oldT,oldVote)
		cp,Out=getVoteCentre(In,Vote,WriteRad_pix)
		diff=(In-Out).sum()
		if (diff==0): # something went wrong with differential voting => go back to absolute voting
			print('Something went wrong with voting')
		In=Out.copy()
		centList.append(cp)
	Cents=np.asarray(centList)
	Xc=X[Cents[:,1],Cents[:,0]]
	Yc=Y[Cents[:,1],Cents[:,0]]
	points=np.column_stack((Xc,Yc)).astype(np.float64)
	return points

def getVoteFullCircle(T,VoteRad):
	Disc=disk(VoteRad)
#	E=cv2.Canny(T.astype(np.uint8),1,1)
	Sx=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,1,0)
	Sy=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,0,1)
	SS,Sang=cv2.cartToPolar(Sx,Sy)
	E=np.where(SS>1,1,0)
	RE=signal.fftconvolve(E.astype(np.uint8),Disc,'same')
	return RE/Disc.sum()
	
# apply gaussian filter on disk. But make sure that you are not voting outside of original disk
#  => should probably rescale the disk. => make smaller disk in large enough matrix and apply
# gaussian on that. 
def getVoteFullCircle_gaus(T,VoteRad,gaus):
	Disc=disk(VoteRad)
#	E=cv2.Canny(T.astype(np.uint8),1,1)
	Sx=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,1,0)
	Sy=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,0,1)
	SS,Sang=cv2.cartToPolar(Sx,Sy)
	E=np.where(SS>1,1,0)
	RE=signal.fftconvolve(E.astype(np.uint8),Disc,'same')
	return RE/Disc.sum()
	
def getVoteFullCircle_Bulk(T,VoteRad):
	Disc=disk(VoteRad)
#	E=cv2.Canny(T.astype(np.uint8),1,1)
	#~ Sx=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,1,0)
	#~ Sy=cv2.Sobel(T.astype(np.uint8),cv2.CV_64F,0,1)
	#~ SS,Sang=cv2.cartToPolar(Sx,Sy)
	#~ E=np.where(SS>1,1,0)
	RE=signal.fftconvolve(T.astype(np.uint8),Disc,'same')
	return RE/Disc.sum()


# This will cast votes in a circle of VoteRadius around all contourpoints
def getVoteFullCircle_old(T2,VoteRad,T=np.zeros((1,1)),Vote=np.zeros((1,1))):
	if (T.sum()==0):
		P=np.zeros(T2.shape,dtype=np.uint8)
		R=cv2.findContours(T2.astype(np.uint8),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
		Res=np.zeros(T2.shape,dtype=np.int32)
		for CC in R[1]:
			CC=CC.reshape((-1,2))

			for cc in CC:
				P[:]=0
			#    cc=CC[150]
				P=cv2.circle(P,(cc[0],cc[1]),VoteRad,1,-1)
				#Pb=cv2.GaussianBlur(P,(kernel,kernel),0)
				Res=Res+P
		return Res
	else:
		try:
			Tdif=T-T2
			Cdif=cv2.findContours((Tdif).astype(np.uint8),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
			Cold=cv2.findContours((T).astype(np.uint8),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)


			AllC_old=Cold[1][0].reshape((-1,2))
			for nC in Cold[1][1:]:
				AllC_old=np.row_stack((AllC_old,nC.reshape((-1,2))))
			#print(AllC_old.shape)

			#print(T2.sum(),len(Cdif[1]),Tdif.sum())
			AllC_dif=Cdif[1][0].reshape((-1,2))
			for nC in Cdif[1][1:]:
				AllC_dif=np.row_stack((AllC_dif,nC.reshape((-1,2))))
			

			P=np.zeros(T.shape,dtype=np.uint8)
			Sub=np.zeros(T.shape,dtype=np.int32)
			Add=np.zeros(T.shape,dtype=np.int32)
			for difE in AllC_dif:
				dis=np.abs(AllC_old-difE)
				if (len(np.argwhere(dis[:,0]+dis[:,1]<0.5))>0):
					P[:]=0
					P=cv2.circle(P,(difE[0],difE[1]),VoteRad,1,-1)
					Sub=Sub+P
				else:
					P[:]=0
					P=cv2.circle(P,(difE[0],difE[1]),VoteRad,1,-1)
					Add=Add+P

			NewVote=Vote.astype(np.int32)+Add-Sub
		except: # something went wrong with the differential voting system => use classic voting!
			NewVote=getVoteFullCircle(T2,VoteRad)
		return NewVote

# This will take the original matrix + the vote => find the elected 
# centre point and substract a zone with WriteRadius from the material matrix T => return T2
def getVoteCentre(T,Vote,WriteRad):
	mm=np.argmax(Vote*T)
	#mm=np.argmax(Vote) # While sometimes better, it can get stuck (happened for free-standing)
	xmax=mm%Vote.shape[1]
	ymax=np.int(mm/Vote.shape[1])
	R2=np.zeros(T.shape,dtype=np.float)
	R2=cv2.circle(R2,(xmax,ymax),WriteRad,1,-1)
	T2=np.where(T-R2>0,1,0)
	return np.array([xmax,ymax]),T2


def makeVoronoiPolygons(points):
	vor = Voronoi(points)
	regions, vertices = voronoi_finite_polygons_2d(vor)
	return regions,vertices

def voronoi_finite_polygons_2d(vor, radius=None):
	"""
	Reconstruct infinite voronoi regions in a 2D diagram to finite
	regions.

	Parameters
	----------
	vor : Voronoi
		Input diagram
	radius : float, optional
		Distance to 'points at infinity'.

	Returns
	-------
	regions : list of tuples
		Indices of vertices in each revised Voronoi regions.
	vertices : list of tuples
		Coordinates for revised Voronoi vertices. Same as coordinates
		of input vertices, with 'points at infinity' appended to the
		end.

	"""

	if vor.points.shape[1] != 2:
		raise ValueError("Requires 2D input")

	new_regions = []
	new_vertices = vor.vertices.tolist()

	center = vor.points.mean(axis=0)
	if radius is None:
		radius = vor.points.ptp().max()*2

	# Construct a map containing all ridges for a given point
	all_ridges = {}
	for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
		all_ridges.setdefault(p1, []).append((p2, v1, v2))
		all_ridges.setdefault(p2, []).append((p1, v1, v2))

	# Reconstruct infinite regions
	for p1, region in enumerate(vor.point_region):
		vertices = vor.regions[region]

		if all([v >= 0 for v in vertices]):
			# finite region
			new_regions.append(vertices)
			continue

		# reconstruct a non-finite region
		ridges = all_ridges[p1]
		new_region = [v for v in vertices if v >= 0]

		for p2, v1, v2 in ridges:
			if v2 < 0:
				v1, v2 = v2, v1
			if v1 >= 0:
				# finite ridge: already in the region
				continue

			# Compute the missing endpoint of an infinite ridge

			t = vor.points[p2] - vor.points[p1] # tangent
			t /= np.linalg.norm(t)
			n = np.array([-t[1], t[0]])  # normal

			midpoint = vor.points[[p1, p2]].mean(axis=0)
			direction = np.sign(np.dot(midpoint - center, n)) * n
			far_point = vor.vertices[v2] + direction * radius

			new_region.append(len(new_vertices))
			new_vertices.append(far_point.tolist())

		# sort region counterclockwise
		vs = np.asarray([new_vertices[v] for v in new_region])
		c = vs.mean(axis=0)
		angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
		new_region = np.array(new_region)[np.argsort(angles)]

		# finish
		new_regions.append(new_region.tolist())

	return new_regions, np.asarray(new_vertices)

def convertToPixels(X,Y,realList):
	X_start = X[0,0]
	X_delta = X[0,1]-X[0,0]
	Y_start = Y[0,0]
	Y_delta = Y[1,0]-Y[0,0]

	if (realList.shape[1]==2):
		CoordX_pix=(np.round((realList[:,0]-X_start)/X_delta)).astype(np.int)
		CoordY_pix=(np.round((realList[:,1]-Y_start)/Y_delta)).astype(np.int)
		coords_pix=np.column_stack((CoordX_pix,CoordY_pix))
		return coords_pix
	elif (realList.shape[1]==1):
		return (np.round((realList)/X_delta)).astype(np.int)  # can only assume X_delta and Y_delta are the same


def makeVoronoiMask(polygon,X,Y):
	M=np.zeros(X.shape,dtype=np.uint8)
#	sssx=[np.argmin(np.abs(xx-ppp[0])) for ppp in polygon]
#	sssy=[np.argmin(np.abs(yy-ppp[1])) for ppp in polygon]
#	kj=np.column_stack((sssx,sssy)) # or other order!!
	kj=convertToPixels(X,Y,polygon)
	M=cv2.fillPoly(M,[kj.astype(np.int32)],1,100)
	return M
	
def makeCircularMask(centre,radius,xx,yy):
	M=np.zeros((len(yy),len(xx)),dtype=np.uint8)
	kj=convertToPixels(X,Y,np.asarray(centre))
	dpix=xx[1]-xx[0]
	rad_pix=int(round(radius/dpix))
	M=cv2.circle(M,(kj[0],kj[1]),rad_pix,1,-1)
	return M
	
def getInfill(points,writeRadius,T,X,Y):
	regions, vertices = makeVoronoiPolygons(points)
	myL=[]
	for i in np.arange(len(points)):
		polygon = vertices[regions[i]]
		MM=makeVoronoiMask(polygon,X,Y)
		myMaterial=(MM*T).sum()
		myL.append(myMaterial)
		
	return np.asarray(myL)
	
def writeTSPfile(filename,pp):
	import scipy.spatial.distance as dist
	DD=dist.pdist(pp,'cityblock')
	
	
	F=open(filename,mode='w')
	F.write('NAME: {}\n'.format(filename[:-4]))
	F.write('TYPE: TSP\n')
	F.write('COMMENT: Nanoscribe Blocks in Ham Path\n')
	F.write('DIMENSION: {}\n'.format(len(pp)+1))
	F.write('EDGE_WEIGHT_TYPE: EXPLICIT\n')
	F.write('EDGE_WEIGHT_FORMAT: UPPER_ROW\n')
	F.write('DISPLAY_DATA_TYPE: NO_DISPLAY\n')
	F.write('EDGE_WEIGHT_SECTION\n')

	pcounter=0
	st=''
	for cc in np.arange(len(pp)): # start with a zero line to transform into ham path problem
		st=st+'0 '
	F.write(st+'\n')

	for c in np.arange(len(pp)-1,0,-1):
		st=''
		for i in np.arange(c):
			st=st+'{:d} '.format(int(DD[pcounter]))
			pcounter=pcounter+1
		F.write(st+'\n')


	F.write('EOF')

	F.close()
	return

def readTSPfile(filename):
	Sol=open(filename,mode='r').readlines()
	OrdList=[]
	for line in Sol[1:]:
		els=line.split()
		for el in els:
			OrdList.append(int(el))
	OO=np.asarray(OrdList)
	if (OO[0]==0):
		OO=OO[1:]-1 # cut away the fictitious first point
	else:
		print('First point is not the fictitious one: investigate')
	return OO
	
def TSPnearest(pp):
	import scipy.spatial.distance as dist
	DD=dist.pdist(pp,'cityblock')
	Dmat=dist.squareform(DD)
	order=(np.argsort(Dmat,axis=1))[:,1:] # cut out the first point as it is always the diagonal point

	pathList=[]
	distList=[]

	for stp in np.arange(len(pp)):
		visList=[stp]
		dist=0
		currP=stp
		while (len(visList)<len(pp)):
			foundNext=False
			runI=0
			while (not foundNext):
				targ=order[currP][runI]
				if (targ in visList):
					#print('already visited {}'.format(targ))
					runI=runI+1
				else:
					#print('found a candidate: {}'.format(targ))
					visList.append(targ)
					dist=dist+Dmat[currP,targ]
					currP=targ
					foundNext=True
	#    VV=np.asarray(visList)
		pathList.append(visList)
		distList.append(dist)

	PP=np.asarray(pathList)
	Dist=np.asarray(distList)

	bestInd=np.argmin(Dist)
	return PP[bestInd],Dist[bestInd]
	
def getVoteCentre_freestanding(T,L,Vote,WriteRad,interior=True):
	if (interior):
		mm=np.argmax(Vote*T)
	else:
		mm=np.argmax(Vote)
	xmax=mm%Vote.shape[1]
	ymax=np.int(mm/Vote.shape[1])
	LL=T*L
	R2=np.zeros(T.shape,dtype=np.float)
	R2=cv2.circle(R2,(xmax,ymax),WriteRad,1,-1)
	II,cnts,hier=cv2.findContours(R2.astype(np.uint8),cv2.RETR_LIST,cv2.CHAIN_APPROX_NONE)
	Conts=cnts[0].reshape((-1,2))
	crossing=np.unique(LL[Conts[:,1],Conts[:,0]])[1:]
	Cross=np.zeros(T.shape,dtype=np.uint8)
	for cc in crossing:
		CrossT=np.where(LL==cc,1,0)
		Cross=Cross+CrossT
	RemMask=np.where(R2-im.binary_dilation(Cross).astype(np.uint8)>0,1,0)
	goodLabs=np.unique(RemMask*LL)[1:]
	
	T2=np.where(T-RemMask>0,1,0)
	return np.array([xmax,ymax]),T2,goodLabs,crossing
	
def makeFreeMaterialMatrix(Tloc,seeds):
	TestIm=np.zeros(Tloc.shape,np.uint8)
	Mask=np.zeros((Tloc.shape[0]+2,Tloc.shape[1]+2),dtype=np.uint8)
	Mask[1:-1,1:-1]=Tloc[:]
	Mask=1-Mask
	for seed in seeds:
		rr,TestIm,Mask,rect=cv2.floodFill(TestIm,Mask,(seed[0].astype(np.int32),seed[1].astype(np.int32)),1)
	return TestIm
	
def getOptimalFreestandingStitchCentres(T,VoteRadius,WriteRadius,X,Y,closing=0):
	dPix=X[0,1]-X[0,0]
	VoteRad_pix=int(round(VoteRadius/dPix))
	WriteRad_pix=int(round(WriteRadius/dPix))
	
	L,n=im.label(T,np.ones((3,3)))
	posLabels=[]
	for i in np.arange(1,n+1):
		Lcoords=np.argwhere(L==i)
		meanLabel=Lcoords.mean(axis=0)
		posLabels.append([meanLabel[1],meanLabel[0],i])
	Pos=np.round(np.asarray(posLabels)).astype(np.int)
	#print(Pos.shape,X.shape)
	realPos=np.column_stack((X[0][Pos[:,0]],Y[:,0][Pos[:,1]]))
	if (len(Pos)==0):
		print('Big problem with algorithm, check debugoutput')
		return L,posLabels,T
	

	In=T.copy()
	if (closing==0):
		Voters=In.copy()
	else:
		#Voters=im.binary_closing(In,iterations=closing).astype(np.uint8)
		Voters=im.binary_dilation(In,iterations=closing).astype(np.uint8)
	
	centList=[]
	labPosList=[]
	labsList=[]
	#VoteList=[]
	while(In.sum()!=0):
		Vote=getVoteFullCircle(Voters,VoteRad_pix)
	 #   VoteList.append(Vote)
		cp,Out,goodLabs,crossing=getVoteCentre_freestanding(In,L,Vote,WriteRad_pix,True)
		diff=(In-Out).sum()
		if (diff==0): # something went wrong with differential voting => go back to absolute voting
			print('Something wrong with voting, will do plan B')
			TempIn=np.where(L==crossing[0],1,0)
			Vote=getVoteFullCircle(TempIn,VoteRad_pix)
			cp,Out,goodLabs,crossing=getVoteCentre_freestanding(In,L,Vote,WriteRad_pix,True)
			diff=(In-Out).sum()
			if (diff==0):
				print('Do not know what to do, will abort')
				return In,Voters,Vote
				break
		In=Out.copy()
		if (closing==0):
			Voters=In.copy()
		else:
			Voters=im.binary_closing(In,iterations=closing).astype(np.uint8)
		centList.append(cp)
		labPosList.append(realPos[goodLabs-1])
		labsList.append(goodLabs)
	Cents=np.asarray(centList)
	Xc=X[Cents[:,1],Cents[:,0]]
	Yc=Y[Cents[:,1],Cents[:,0]]
	points=np.column_stack((Xc,Yc)).astype(np.float64)
	return points,labPosList,labsList#,VoteList,Voters
