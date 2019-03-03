from math import sin, cos, pi
from sys import exit
import numpy

def markerResidual(cent, Markers_, cTilt, sTilt, 
        transX, transY, rotInPlane, 
        isoMag=None, 
        dBeam=None, 
        dMagnFocus=None, dRotFocus=None,
        equationSet=False):
    """
    calculate residual of a marker model given the marker coords

    @param cent: x,y coordinates of tiltaxis rotation stage.
    @type cent: numpy.array
    @param Markers_: Markers
    @type Markers_: Markers
    @param transX: translations in X 
    @type transX: vector (ntilt)
    @param transY: translations in Y 
    @type transY: vector (ntilt)
    @param rotInPlane: rotation of tilt axis in each projection (against X-axis in deg)
    @type rotInPlane: array
    @param dmag: delta magnification (vector; dim: ntilt)
    @type dmag: vector (ntilt)
    @param dbeam: beam inclination
    @type dbeam: float
    @param dmagnfoc: linear magnification change as function of defocus (experimental)
    @type dmagnfoc: float
    @param drotfoc: linear image rotation change as function of defocus
    @type drotfoc: float
    @param equationSet: Compute differences instead of scalar function (sum of squares)
    @type quationSet: logical
    @return: squared error of fit in 3D (or array of deviations for equationSet==True)
    @rtype: float (or array for equationSet==True)
    @author: foerster

    tilt geom:

    [cos(the_itilt)  0  -sin(the_itilt) ]  [ x_imark ]   
    [ 0              1             0    ]  [ y_imark ]
                                           [ z_imark ]


    = [ cos(-psi+dpsi_itilt-pi/2) sin(-psi+dpsi_itilt-pi/2)] isoMag_itilt [markx_itilt,imark] - [transX_itilt]
      [-sin(-psi+dpsi_itilt-pi/2) cos(-psi+dpsi_itilt-pi/2)]              [marky_itilt,imark]   [transY_itilt]


    """
    

    from math import sin, cos, pi
    from sys import exit
    from numpy import mean
    from pytom.tools.maths import rotate_vector2d


    rotInPlane *=0.

    ntilt = len(transX)
    nmark = len(Markers_)
    if equationSet:
        Dev = []

    # approximate tilt axis = mean
    meanpsi = mean(rotInPlane)
    cmeanpsi = cos(- meanpsi/180.*pi - pi/2.)
    smeanpsi = sin(- meanpsi/180.*pi - pi/2.)

    # pre-compute sin and cos
    cpsi = numpy.array(ntilt*[0.])
    spsi = numpy.array(ntilt*[0.])
    for (ii, psi) in enumerate(rotInPlane):
        cpsi[ii] = cos(- psi/180.*pi - pi/2.)
        spsi[ii] = sin(- psi/180.*pi - pi/2.)
    if dBeam:
        cdbeam = cos(dBeam)
        sdbeam = sin(dBeam)

    ndif =0 #counter for datapoints

    # marker loop
    residual = 0.
    errors = []
    for (imark, Marker) in enumerate(Markers_):
        errors.append(0.)
        ind_dif = 0
        markCoord = Marker.get_r()
        # marker positions in 3d model rotated on approximate tilt axis
        zmod = markCoord[2]
        [xmod, ymod] = rotate_vector2d([markCoord[0],markCoord[1]], 
	        cmeanpsi, smeanpsi)

        # tilt loop
        for iproj in range(0, ntilt):
            # check for not-clicked Markers (have coordinates -1,-1)
            if ( (Markers_[imark].xProj[iproj] > -1.) and 
                    (Markers_[imark].yProj[iproj] > -1.) ):
                ndif= ndif +1

                # rotate clicked markers back and inverse shift
                xmark = Markers_[imark].xProj[iproj] - cent[0]- transX[iproj]
                ymark = Markers_[imark].yProj[iproj] - cent[1]- transY[iproj]
                
                [x_meas, y_meas] = rotate_vector2d([xmark,ymark], cpsi[iproj], spsi[iproj])

                # additional variable magnification
                try:
                    if isoMag:
                        x_meas = x_meas * isoMag[iproj]
                        y_meas = y_meas * isoMag[iproj]
                
                # model projection

                # beam inclination
                    if dBeam:
                        x_proj = cTilt[iproj]*xmod - sTilt[iproj]*sdbeam*ymod - sTilt[iproj]*cdbeam*zmod
                        y_proj = sTilt[iproj]*sdbeam*xmod + ( cdbeam**2+sdbeam**2*cTilt[iproj]*ymod + 
                                                              cdbeam*sdbeam*(1-cTilt[iproj])*zmod )
                except:
                    x_proj = cTilt[iproj] * xmod - sTilt[iproj]*zmod
                    y_proj = ymod
        
                # x-dependent magnification of model
                try:
                    if dMagnFocus:
                        y_proj_dmag = y_proj
                        tmp         = dMagnFocus*x_proj
                        x_proj      = (1+.5*tmp)*x_proj
                        y_proj_dmag = (1+tmp)*y_proj
                        y_proj      = y_proj_dmag
                except:
                    continue
                # score

                deltaX = x_meas-x_proj
                deltaY = y_meas-y_proj
                ## Warning for large differences
                #if (abs(deltaX) >100.) or (abs(deltaY) > 100.):
                #    print "Huge Difference for marker oberved: "
                #    print "  iproj="+str(iproj)+", imark="+str(imark+1)
                #    print "  deltaX = "+str(deltaX)
                #    print "  deltaY = "+str(deltaY)
                
                residual = deltaX**2 + deltaY**2 + residual
                if equationSet:
                    Dev.append(deltaX)
                    Dev.append(deltaY)
                ind_dif += 1
                
                
                errors[-1] +=  numpy.sqrt(deltaX**2 + deltaY**2) 
        errors[-1] = errors[-1]/float(ind_dif)**2
    # normalize and prepare output
    normf = 1. /(ndif - nmark)
    # residual = sqrt(residual/(ndif - size(Matrixmark,3)));
    residual = residual*normf
    #print error, residual
    if equationSet:
        return Dev
    else:
        return errors # residual


def alignmentFixMagRot( Markers_, cTilt, sTilt,
        ireftilt, irefmark=1, r=None, imdim=2048, handflip=False,
	mute=True):
    """
    compute alignment analytically (constant mag. and tilt axis)

    @param cTilt: cosine of tilt angles
    @type cTilt: numpy array
    @param sTilt: sine of tilt angles
    @type sTilt: numpy array
    @param r: coordinates of reference marker
    @type r: numpy.array
    @param handflip: flip handedness 
    @type handflip: logical
    @param mute: turn output silent
    @type mute: L{bool}

    @author: FF
    """
    from math import atan, tan, pi, sqrt

    ntilt = len(cTilt)
    nmark = len(Markers_)
    
    if ( r==None):
        if not mute:
            print "Assign reference marker to default value"
	r = numpy.array(3*[0.])
        r[0] = Markers_[irefmark-1].xProj[ireftilt-1]
        r[1] = Markers_[irefmark-1].yProj[ireftilt-1]
        r[2] = float(imdim/2 +1)
    else:
	r = Markers_[irefmark-1].get_r()
        
    #   calculate means of difference vectors
    meanx=numpy.array(nmark*[0.])
    meany=numpy.array(nmark*[0.])
    norm =numpy.array(nmark*[0.])
    for (imark,Marker) in enumerate(Markers_):
        for itilt in range(0,ntilt):
            if ( (Marker.xProj[itilt] > -1.) and (Markers_[irefmark-1].xProj[itilt] > -1.)): #allow overlapping MPs
                meanx[imark] = meanx[imark] + Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt]
                meany[imark] = meany[imark] + Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt]
                norm[imark]  = norm[imark] +1;
        meanx[imark] = meanx[imark] / norm[imark]
        meany[imark] = meany[imark] / norm[imark]
    #   calculate some sums for determination of tilt axis azimuth
    #   e.g. sumxx = sum(imark) sum(itilt) (delta(imark,itilt)-deltaaverage(imark))^2 
    sumxx=0.
    sumyy=0.
    sumxy=0.
    for (imark,Marker) in enumerate(Markers_):
        for itilt in range(0,ntilt):
            
            if ( (Marker.xProj[itilt] > -1.) and (Markers_[irefmark-1].xProj[itilt] > -1.)): #allow overlapping MPs
                sumxx = (Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt] -meanx[imark])**2 + sumxx
                sumyy = (Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt] -meany[imark])**2 + sumyy
                sumxy = ((Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt] -meanx[imark]) * 
		         (Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt] -meany[imark]) + sumxy)
    
    #  calculate azimuth :
    #  minimize sum( (x(i)*sin(psi) + y(i)*cos(psi))^2) =: Min(F)  --- linear regression  
    #  d/d(psi) (F) leads to:
    psi = 0.
    if sumxx-sumyy > 0.0000001:
        psi = 0.5*atan(2*sumxy/(sumxx-sumyy))
    if (sumxx > sumyy):
        psi = psi - 0.5*pi*cmp(psi, 0)
    if handflip:
        psi = psi + pi
    psiindeg = psi*180/pi
    
    #  psi 2?!
    #fact = 1 /(1 + (sumxy/sumxx)**2)
    
    # calculate deviations
    cpsi = cos(psi)
    spsi = sin(psi)
    tpsi = tan(psi)
    sumt=0.
    sumxx=0.
    ndif =0
    distLine = numpy.array(nmark*[ntilt*[0.]])
    for (imark,Marker) in enumerate(Markers_):
        for itilt in range(0,ntilt):
            if ( (Marker.xProj[itilt] > -1.) and (Markers_[irefmark-1].xProj[itilt] > -1.)): #allow overlapping MPs
                if (imark != irefmark-1):
                    ndif= ndif +1 # count all markers except for refmark
		distLine[imark,itilt] = ((Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt] -meanx[imark])*cpsi + 
		        (Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt] -meany[imark])*spsi)
                sumxx = distLine[imark,itilt]**2 +sumxx
    sigma = sqrt(sumxx/(ndif - nmark ));
    #   deviation as angle in deg
    if not mute:
        print('Number of tilts:.............. = %3d' %ntilt)
        print('reference tilt index:......... = %3d' %ireftilt)
        print('Total number of marker points  = %3d' %nmark)
        print('Index of reference point:..... = %3d' %irefmark)
        print('Tilt axis azimuth:............ = %6.2f deg' %psiindeg)
        print('RMS error fit:................ = %4.2f pix' %sigma)
    
    #   ---- 2nd part: determination of shifts ----
    
    #determine 3D coordinates of markers -> solve eqn system
    x = numpy.array(nmark*[0.])
    y = numpy.array(nmark*[0.])
    z = numpy.array(nmark*[0.])
    if not mute:
        print('Coordinates of reference marker: %7.1f, %7.1f, %7.1f'%(r[0], r[1], r[2]))
        print('Difference vectors of marker points:')
    for (imark,Marker) in enumerate(Markers_):
        sumxx=0.; sumyy=0.; sumxy=0.; sumyx=0.; salpsq = 0.; scalph = 0.;
        P = numpy.array(3*[3*[0.]])
        P_t = numpy.array(3*[3*[0.]])
        temp = numpy.array(3*[0.])
	norm[imark]=0.


        for itilt in range(0,ntilt):

            if ( (Marker.xProj[itilt] > -1.) and (Markers_[irefmark-1].xProj[itilt] > -1.)): #allow overlapping MPs
                norm[imark]= norm[imark]+1;
                salpsq = salpsq + sTilt[itilt]**2; #sum sin^2
                scalph = scalph + cTilt[itilt]*sTilt[itilt] #sum cos*sin
                #sum delta x * cos
                sumxx = sumxx+ (Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt])* cTilt[itilt]
                #sum delta(y)*cos
                sumyy = sumyy + (Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt])* cTilt[itilt]
                #sum delta(x)*sin
                sumxy = sumxy+ (Marker.xProj[itilt] - Markers_[irefmark-1].xProj[itilt])* sTilt[itilt]
                #sum delta(y)*sin
                sumyx = sumyx + (Marker.yProj[itilt] - Markers_[irefmark-1].yProj[itilt])* sTilt[itilt]
        P[0,0] = norm[imark] - salpsq*spsi**2
        P[0,1] = salpsq*cpsi*spsi
        P[0,2] = scalph*spsi
        P[1,0] = P[0,1]
        P[1,1] = norm[imark] - salpsq*cpsi**2
        P[1,2] = -scalph*cpsi
        P[2,0] = P[0,2]
        P[2,1] = P[1,2]
        P[2,2] = salpsq

        dt = numpy.linalg.det(P)
        temp[0] = ( (sumxx*spsi-sumyy*cpsi)*spsi + (cpsi*meanx[imark]+
	                  spsi*meany[imark])*cpsi*norm[imark] )
        temp[1] = ( -(sumxx*spsi-sumyy*cpsi)*cpsi + (cpsi*meanx[imark]+
	                  spsi*meany[imark])* spsi*norm[imark] )
        temp[2] = sumxy*spsi - sumyx*cpsi
        if (dt != 0):
            P_t=P.__copy__()
            P_t[0,0] = temp[0]
            P_t[1,0] = temp[1]
            P_t[2,0] = temp[2]
            x[imark] = numpy.linalg.det(P_t)/dt
            P_t=P.__copy__()
            P_t[0,1] = temp[0]
            P_t[1,1] = temp[1]
            P_t[2,1] = temp[2]
            y[imark] = numpy.linalg.det(P_t)/dt
            P_t=P.__copy__()
            P_t[0,2]=temp[0]
            P_t[1,2]=temp[1]
            P_t[2,2]=temp[2]
            z[imark] = numpy.linalg.det(P_t)/dt;
            if not mute:
                print('     %3d - %3d :.............. = %7.1f, %7.1f, %7.1f'
                    %(imark+1, irefmark, x[imark], y[imark], z[imark]))
            x[imark] = x[imark] + r[0] - imdim/2. -1. # move to center
            y[imark] = y[imark] + r[1] - imdim/2. -1.
            z[imark] = z[imark] + r[2] - imdim/2. -1.
        else:
            if not mute:
                print('Marker '+str(imark)+' : undefined! det = 0! Click more!')
            x[imark] = 1000000
    
    if (len(y) < len(x)):
        for ii in range(len(y-1),len(x)):
	    y[ii] = 1000000
    if (len(z) < len(x)):
        for ii in range(len(z-1),len(x)):
	    z[ii] = 1000000
    
    # determination of shifts
    projX=numpy.array(ntilt*[0.])
    projY=numpy.array(ntilt*[0.])
    shiftX=numpy.array(ntilt*[0.])
    shiftY=numpy.array(ntilt*[0.])
    diffX = numpy.array(ntilt*[nmark*[0.]])
    diffY = numpy.array(ntilt*[nmark*[0.]])
    shiftVarX  = numpy.array(ntilt*[0.])
    shiftVarY  = numpy.array(ntilt*[0.])
    for itilt in range(0,ntilt):
        sumxx = 0.
        sumyy = 0.
        ndif  = 0
        for (imark,Marker) in enumerate(Markers_):
            projX[itilt] = (x[imark]*(spsi**2*cTilt[itilt]+cpsi**2) + 
                             y[imark]*spsi*cpsi*(1-cTilt[itilt]) + 
                             z[imark]*spsi*sTilt[itilt] + imdim/2. + 1.)
            projY[itilt] = (x[imark]*spsi*cpsi*(1-cTilt[itilt]) +
                             y[imark]*(cpsi**2*cTilt[itilt]+spsi**2) -
                             z[imark]*cpsi*sTilt[itilt] + imdim/2. + 1.)
            
            
            if ( (Marker.xProj[itilt] > -1.) and (x[imark]!=1000000) ):
                ndif = ndif + 1
                diffX[itilt,imark] = Marker.xProj[itilt] - projX[itilt]
                diffY[itilt,imark] = Marker.yProj[itilt] - projY[itilt]
            else:
                diffX[itilt,imark] = 0.
                diffY[itilt,imark] = 0.
            sumxx = sumxx + diffX[itilt,imark]
            sumyy = sumyy + diffY[itilt,imark]
        # mean values of shifts
        if (ndif > 0):
            shiftX[itilt] = sumxx / ndif
            shiftY[itilt] = sumyy / ndif
        else:
            shiftX[itilt] = 1000000
            shiftY[itilt] = 1000000
        
        # deviations of individual shift from mean
        sumxx = 0.
        sumyy = 0.
        for (imark,Marker) in enumerate(Markers_):
            if ( (Marker.xProj[itilt] > -1.) and (x[imark]!=1000000) ):
                sumxx = sumxx + (diffX[itilt,imark]-shiftX[itilt])**2
                sumyy = sumyy + (diffY[itilt,imark]-shiftY[itilt])**2
        if (ndif > 1):
            shiftVarX[itilt]=sqrt(sumxx/(ndif-1))
            shiftVarY[itilt]=sqrt(sumyy/(ndif-1))
        else:
            shiftVarX[itilt]=None
            shiftVarY[itilt]=None

        

    return(psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY, 
           shiftVarX, shiftVarY)


def simulate_markers(markCoords, tiltAngles, tiltAxis=-76.71, ireftilt=None,
        ampTrans=100., ampRot=.5, ampMag=.01, dBeam=None, dMagnFocus=None):
    """
    simulate marker coordinates for testing alignment routines

    (Markers, TiltSeries, TiltAlignmentParas) = simulate_markers(markCoords, \
    tiltAngles, tiltAxis=-76.71, ireftilt=None,\
    ampTrans=100., ampRot=.5, ampMag=.01, dBeam=None)

    @param markCoords: marker coordinates
    @type markCoords: array
    @param tiltAngles: tilt angles
    @type tiltAngles: array
    @param tiltAxis: approximate tilt axis
    @type tiltAxis: float
    @param ampTrans: amplitude of translation distortions (gaussian)
    @type ampTrans: float
    @param ampRot: amplitude of rotation distortions (gaussian)
    @type ampRot: float
    @param ampMag: amplitude of maginification distortions (gaussian)
    @type ampMag: float
    @param dBeam: beam inclination
    @type dBeam: float
    @param dMagnFocus: magnification-dependent on distance from tilt axis
    @type dMagnFocus: float
    @return: Markers, TiltSeries, TiltAlignmentParas

    """
    from random import gauss
    from math import cos, sin, pi
    from pytom.reconstruction import TiltAlignmentStructures
    from pytom.tools.maths import rotate_vector2d

    if dBeam:
        print('dBeam=%5.3f' %dBeam)
    cent= [1025,1025]
    if dBeam:
        cdbeam = cos(dBeam*pi/180.)
        sdbeam = sin(dBeam*pi/180.)

    ntilt = len(tiltAngles)
    rot = ntilt*[tiltAxis]
    srot = ntilt*[0.]
    crot = ntilt*[1.]
    transX = ntilt*[0.]
    transY = ntilt*[0.]
    mag = ntilt*[1.]
    cTilt = ntilt*[0.]
    sTilt = ntilt*[0.]

    if not ireftilt:
        ireftilt = tiltAngles.argmin()+1
    #generate Markers
    nmarks = len(markCoords)
    Markers = []
    for (imark, markCoord) in enumerate(markCoords):
        Marker = TiltAlignmentStructures.Marker(range(0,ntilt))
	Marker.set_r(markCoord)
	Markers.append(Marker)

    TiltAlignmentParas = TiltAlignmentStructures.TiltAlignmentParameters( 
            dmag=True, drot=True, 
            dbeam=bool(dBeam),
            finealig=False, finealigfile=None,
            grad=False,
            irefmark=1, ireftilt=ireftilt, r= [0., 0., 0.],
            cent= cent,
            handflip=False,
            optimizer='fmin', maxIter=2000, leastsq=False)
    TiltSeries = TiltAlignmentStructures.TiltSeries( 
            tiltSeriesName=None,
            TiltAlignmentParas=TiltAlignmentParas, 
	    alignedTiltSeriesName=None,
            markerFileName=None, firstProj=1, lastProj=ntilt, projIndices=None, 
	    tiltSeriesFormat='em')
    TiltSeries.createEmptyProjections()

    #compute 'alignment'
    for (itilt, tiltAngle) in enumerate(tiltAngles):
        rot[itilt] = tiltAxis + gauss(0.,ampRot)
	crot[itilt] = cos(rot[itilt]/180.*pi + pi/2.)
	srot[itilt] = sin(rot[itilt]/180.*pi + pi/2.)
	if itilt+1 == ireftilt:
	    mag[itilt] = 1.
	    transX[itilt] = 0.
	    transY[itilt] = 0.
	else:
            mag[itilt] = gauss(1.,ampMag)
	    transX[itilt] = gauss(0., ampTrans)
	    transY[itilt] = gauss(0., ampTrans)
	cTilt[itilt] = cos(tiltAngle/180.*pi)
	sTilt[itilt] = sin(tiltAngle/180.*pi)
	TiltSeries._ProjectionList[itilt].setTiltAngle( tiltAngle)
	TiltSeries._ProjectionList[itilt].setAlignmentTransX(transX[itilt])
	TiltSeries._ProjectionList[itilt].setAlignmentTransY(transY[itilt])
	TiltSeries._ProjectionList[itilt].setAlignmentRotation(rot[itilt])
	TiltSeries._ProjectionList[itilt].setAlignmentMagnification(mag[itilt])

    # convert marker 3d coordinates to approximate tilt axis
    cpsi = cos(-(tiltAxis/180.*pi+pi/2.))
    spsi = sin(-(tiltAxis/180.*pi+pi/2.))

    for (imark, Marker) in enumerate(Markers):
        markCoord = Marker.get_r()
        # marker positions in 3d model rotated on approximate tilt axis
        zmod = markCoord[2]
        [xmod, ymod] = rotate_vector2d(markCoord, cpsi, spsi)

        # simulate corresponding marker positions in projs
        for (itilt, tiltAngle) in enumerate(tiltAngles):
            # model projection
            ## beam inclination
            if dBeam:
                x_proj = cTilt[itilt]*xmod - sTilt[itilt]*sdbeam*ymod - sTilt[itilt]*cdbeam*zmod
                y_proj = sTilt[itilt]*sdbeam*xmod + ( cdbeam**2+sdbeam**2*cTilt[itilt]*ymod + 
                         cdbeam*sdbeam*(1-cTilt[itilt])*zmod )
            else:
                x_proj = cTilt[itilt] * xmod - sTilt[itilt]*zmod
                y_proj = ymod
       
            # x-dependent magnification of model
            if dMagnFocus:
                y_proj_dmag = y_proj
                tmp         = dMagnFocus*x_proj
                x_proj      = (1+.5*tmp)*x_proj
                y_proj_dmag = (1+tmp)*y_proj
                y_proj      = y_proj_dmag
	    
            # isotropic magnification
	    if mag[itilt] != 1.:
                x_proj = x_proj * 1./mag[itilt]
                y_proj = y_proj * 1./mag[itilt]

            # image rotation
            [x_proj, y_proj] = rotate_vector2d([x_proj, y_proj], 
	                         crot[itilt], srot[itilt])

	    # image shift and back from center
	    x_proj = x_proj + transX[itilt] + cent[0]
	    y_proj = y_proj + transY[itilt] + cent[1]

	    #set marker coordinates in projection
            Markers[imark].set_xProj( itilt, x_proj)
            Markers[imark].set_yProj( itilt, y_proj)
    TiltSeries._Markers = Markers

    return Markers, TiltSeries, TiltAlignmentParas


def refineMarkerPositions(tiltSeriesName, markerFileName, firstProj, lastProj, finealigfile,
                          ireftilt=0, irefmark=1, cent=[1025,1025],dimBox=32, verbose=False):
    """
    refine coordinates of markers

    @param tiltSeriesName: name of tilt series (<Name>_index.em or <Name>_index.mrc)
    @type tiltSeriesName: L{str}
    @param markerFileName: name of marker file
    @type markerFileName: L{str}
    @param firstProj: index of first projection
    @type firstProj: L{int}
    @param lastProj: index of last projection
    @type lastProj: L{int}
    @param finealigfile: output file with refined marker coordinates
    @type finealigfile: L{str}
    @param dimBox: dimension of box for markers
    @type dimBox: L{int}
    @param verbose: verbose for debugging
    @type verbose: L{bool}

    @author: FF
    """
    from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment

    #read data
    MyTiltAlignmentParas=TiltAlignmentParameters(
        dmag=False, drot=False, dbeam=False,
	finealig=True, finealigfile=finealigfile,
	grad=False,
	irefmark=irefmark, ireftilt=ireftilt, r=None, cent=cent,
	handflip=False, optimizer='leastsq', maxIter=0)
    MyTiltSeries= TiltSeries(tiltSeriesName=tiltSeriesName,
        TiltAlignmentParas=MyTiltAlignmentParas,
        alignedTiltSeriesName='dummy',
        markerFileName=markerFileName,
        firstProj=firstProj, lastProj=lastProj)
    MyTiltAlignment = TiltAlignment(MyTiltSeries)
    #MyTiltAlignment.computeCarseAlignment( MyTiltSeries, mute=True)
    MyTiltAlignment.refineMarkerPositions( TiltSeries_=MyTiltSeries, dimBox=dimBox, 
        finealigfile=finealigfile, verbose=False)


def readIMODmarkerfile(markerfile, binning=1):
    """
    read Markers from WIMP file generated by IMOD, taking binning into account

    @param markerfile: wimp file
    @type markerfile: L{str}
    @param binning: IMOD pre-binning (default: 1 = no binning). binning=2: 2x2 pixels -> 1 pixel, \
    binning=3: 3x3 pixels -> 1 pixel, etc. Pre-binning is applied on projections prior to marker tracking in IMOD.
    @type binning: int or float

    @return: Markers - list of markers
    @rtype: list
    """
    from pytom.reconstruction.TiltAlignmentStructures import Marker
    if binning < 1:
        raise ValueError('readIMODmarkerfile: binning must be >= 1!')
    scalefac = float(binning)
    # read markerfile
    fh=open(markerfile)
    # skip first 6 lines
    for ii in range(0,6):
        lines = fh.readline()
    tmp = fh.readline().split()
    ntilt = int(tmp[3])

    markers = []
    marker = Marker(projIndices=range(1,ntilt+1))
    markers.append(marker)

    lines = fh.readlines()
    fh.close()
    ii = -1
    for line in lines:
        tmp = line.split()
	if len(tmp) == 0:
	    continue
	elif tmp[0] == 'Object':
	    ii = ii + 1
            marker = Marker(projIndices=range(1,ntilt+1))
	    markers.append(marker)
        elif tmp[0] == '#':
	    continue
	elif tmp[0] == 'Display':
	    continue
	elif tmp[0] == 'END':
	    continue
	else:
	    marker.set_xProj( iproj=int(float(tmp[3])), xProj=float(tmp[1]) * scalefac)
	    marker.set_yProj( iproj=int(float(tmp[3])), yProj=float(tmp[2]) * scalefac)
    return markers

def readIMODtiltAngles(tltfile):
    """
    read tilt angles from IMOD file

    @param tltfile: file containing tilt angles
    @type tltfile: str

    @return: list
    @author: FF
    """
    # read tiltangles
    th = open(tltfile)
    lines = th.readlines()
    tiltangles = []
    for line in lines:
        tiltangles.append(float(line))
    th.close()
    return tiltangles

def getIMODpreshifts(prexgfile):
    """
    read preshifts from IMOD prexgfile

    @param prexgfile: file containing shifts
    @type prexgfile: str
    @return: shiftX, shiftY
    @rtype: lists

    @author: FF
    """
    fh = open(prexgfile)
    lines = fh.readlines()
    shiftX = []
    shiftY = []
    for line in lines:
        tmp = line.split()
        shiftX.append(float(tmp[4]))
        shiftY.append(float(tmp[5]))
    return shiftX, shiftY

def applyPreshiftsToMarkers( markers, shiftX, shiftY):
    """
    apply inverted IMOD preshifts to markers - in imod these are always unbinned
    @param markers: markers (list or Marker objects)
    @type markers: list
    @param shiftX: shifts in X
    @type shiftX: list
    @param shiftY: shifts in Y
    @type shiftY: list

    @author: FF
    """
    if markers[0].get_numberOfProjections() != len(shiftX):
        raise ValueError('Dimensions (=number of projections) of shiftX and markers do not match')
    for marker in markers:
        for ii in range(0,len(shiftX)):
	    x = marker.get_xProj( ii)
	    y = marker.get_yProj( ii)
	    if (x > -1) and (y > -1):
	        x = x - shiftX[ii]
	        y = y - shiftY[ii]
		marker.set_xProj(ii, x)
		marker.set_yProj(ii, y)
    return markers

