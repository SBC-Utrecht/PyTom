def refineMarkerPositions( tiltSeriesName, markerFileName, firstProj, 
        lastProj, finealigfile, 
        dimBox=32):
    """
    refine coordinates of markers
    """
    from pytom.image2D.imageStructures import ImageStack
    from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.basic.functions import initSphere, taper_edges

    # prepare mask
    #mask = initSphere(sizeX=dimBox, sizeY=dimBox, sizeZ=1, radius=dimBox/5.,
    #          smooth=dimBox/5., maxradius=0, cent=None)
    mask = None

    #read data
    MyTiltAlignmentParas=TiltAlignmentParameters(
        dmag=False, drot=False, dbeam=False,
	finealig=True, finealigfile=finealigfile,
	grad=False,
	irefmark=1, ireftilt=21, r=None, cent= [1025,1025],
	handflip=False,
	optimizer='leastsq', maxIter=0)
    MyTiltSeries= TiltSeries(tiltSeriesName=tiltSeriesName,
        TiltAlignmentParas=MyTiltAlignmentParas,
        alignedTiltSeriesName='dummy',
        markerFileName=markerFileName,
        firstProj=firstProj, lastProj=lastProj)
    MyTiltAlignment = TiltAlignment(MyTiltSeries)
    MyTiltAlignment.computeCoarseAlignment( MyTiltSeries)
    #MyTiltAlignment.alignFromFiducials()

    for (imark, marker) in enumerate(MyTiltSeries._Markers):
        markerImageStack = ImageStack(verbose=True)
        for itilt in range(0,len(marker._projIndices)):
            x = round(marker.get_xProj(itilt))
            y = round(marker.get_yProj(itilt))
	    if (x > -1) and ( y > -1):
	        proj = MyTiltSeries._ProjectionList[itilt]
	        filename = proj.getFilename()
                #copy data to ImageStack
                markerImageStack.addImageFromFile(filename=filename, 
	            boxCoords=[x-dimBox/2-2,y-dimBox/2-2,0], 
	            dims=[dimBox,dimBox,1], shiftX=0, shiftY=0, rotation=0, 
		    appliedShiftX=0, appliedShiftY=0, index=itilt)
	    else:
	        proj = MyTiltSeries._ProjectionList[itilt]
	        print("marker not clicked in "+proj.getFilename())
        #markerImageStack.normalize( normtype="StdMeanInMask", mask=mask)
        markerImageStack.bandpass( lowfreq=dimBox/32., hifreq=6.*dimBox/32., 
	    smooth=2.*dimBox/32., bpf=None)
        markerImageStack.normalize( normtype="StdMean")
	# smoothen edges
        markerImageStack.taper_edges( width=dimBox/8.)
        markerImageStack.exMaxAlign( niter=10, mask=mask)
	markerImageStack.writeWorkingCopies( filename="Marker"+str(imark)+'.em')
	markerImageStack.writeAverage( filename='av_'+str(imark)+".em")

	# set new coordinates 
	irun = 0
        for itilt in range(0,len(marker._projIndices)):
            x = round(marker.get_xProj(itilt))
            y = round(marker.get_yProj(itilt))
	    if ( x > -1 ) and ( y > -1 ):
	        marker.set_xProj(itilt, x+ markerImageStack.images[irun].shiftX)
	        marker.set_yProj(itilt, y+ markerImageStack.images[irun].shiftY)
	        irun = irun + 1
    MyTiltAlignment.computeCoarseAlignment( MyTiltSeries)
	


