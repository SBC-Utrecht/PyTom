'''
Created July/Aug 2014

@author: FF
'''

from pytom.basic.correlation import nxcc
import scipy.optimize

class Alignment:
    def __init__(self, vol1, vol2, score, mask=None, iniRot=None, iniTrans=None,
                 opti='fmin_powell', interpolation='linear', verbose=False):
        """
        alignment of a particle against a reference

        @param vol1: (constant) volume
        @type vol1: L{pytom.lib.pytom_volume.vol}
        @param vol2: volume that is matched to reference
        @type vol2: L{pytom.lib.pytom_volume.vol}
        @param score: score for alignment - e.g., pytom.basic.correlation.nxcc
        @type score: L{pytom.basic.correlation}
        @param mask: mask correlation is constrained on
        @type mask: L{pytom.lib.pytom_volume.vol}
        @param iniRot: initial rotation of vol2
        @type iniRot: L{pytom.basic.Rotation}
        @param iniTrans: initial translation of vol2
        @type iniTrans: L{pytom.basic.Shift}
        @param opti: optimizer ('fmin_powell', 'fmin', 'fmin_cg', 'fmin_slsqp', 'fmin_bfgs')
        @param interpolation: interpolation type - 'linear' (default) or 'spline'
        @type interpolation: str
        @type opti: L{str}

        @author: FF
        """
        from pytom.basic.normalise import normaliseUnderMask, mean0std1
        from pytom.tools.macros import volumesSameSize
        from pytom.lib.pytom_volume import vol
        from pytom.basic.structures import Rotation, Shift
        assert isinstance(interpolation, str), "interpolation must be of type str"

        self.verbose = verbose
        if not volumesSameSize(vol1,vol2):
            raise RuntimeError('Vol1 and vol2 must have same size!')

        # normalize constant volume
        if mask: 
            (v,p) = normaliseUnderMask( vol1,mask)
        else:
            v = mean0std1( vol1, True)

        self.vol1 = v
        self.vol2 = vol2
        self.rotvol2 = vol(self.vol1.sizeX(), self.vol2.sizeY(), self.vol2.sizeZ())
        self.mask = mask
        
        if iniRot is None:
            iniRot=Rotation()
        if iniTrans is None:
            iniTrans=Shift()
        self.rot_trans = self.transRot2vector( rot=iniRot, trans=iniTrans)

        self.score = score
        self.val = -100000.
        self.centX = int(self.vol1.sizeX()//2)
        self.centY = int(self.vol1.sizeY()//2)
        self.centZ = int(self.vol1.sizeZ()//2)
        self.binning = 1
        self.interpolation = interpolation

        # set optimizer
        self.opti = opti
        if opti=='fmin':
            self.optimizer = scipy.optimize.fmin
        elif opti=='fmin_slsqp':
            self.optimizer = scipy.optimize.fmin_slsqp
        elif opti=='fmin_cg':
            self.optimizer = scipy.optimize.fmin_cg
        elif opti=='fmin_bfgs':
            self.optimizer = scipy.optimize.fmin_bfgs
        elif opti=='fmin_powell':
            self.optimizer = scipy.optimize.fmin_powell
        else:
            raise TypeError('opti must be of type str')

    def transRot2vector(self, rot, trans):
        """
        convert rotation and translation to 6-dimensional vector
        @param rot: rotation of vol2
        @type rot: L{pytom.basic.Rotation}
        @param trans: translation of vol2
        @type trans: L{pytom.basic.Shift}
        @return: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @rtype: L{list}
        @author: FF

        """
        self.rot_trans = 6*[0.]
        for ii in range(0,3):
            self.rot_trans[ii] = rot[ii]
        for ii in range(3,6):
            self.rot_trans[ii] = trans[ii-3]
        return self.rot_trans

    def vector2transRot(self, rot_trans):
        """
        convert 6-dimensional vector to rotation and translation

        @param rot_trans: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @type rot_trans: L{list}
        @return: rotation of vol2, translation of vol2
        @rtype: L{pytom.basic.Rotation}, L{pytom.basic.Shift}
        @author: FF
        """
        from pytom.basic.structures import Rotation, Shift

        rot=Rotation()
        trans=Shift()
        for ii in range(0,3):
            rot[ii] = self.rot_trans[ii]
        for ii in range(3,6):
            trans[ii-3] = self.rot_trans[ii]
        return rot, trans

    def set_rot(self, rot):
        """
        set initial rotation

        @param rot: rotation of vol2
        @type rot: L{pytom.basic.Rotation}
        """
        self.rot_trans[0] = rot[0]
        self.rot_trans[1] = rot[1]
        self.rot_trans[2] = rot[2]

    def set_trans(self, trans):
        """
        set initial translation

        @param trans: translation of vol2
        @type trans: L{pytom.basic.Shift}
        """
        self.rot_trans[3] = trans[0]
        self.rot_trans[4] = trans[1]
        self.rot_trans[5] = trans[2]

    def set_searchVol(self, vol1):
        """
        set search volume (vol1 internally)

        @param vol1: search volume
        @type vol1: L{pytom.lib.pytom_volume.vol}

        """
        from pytom.basic.normalise import normaliseUnderMask, mean0std1

        if self.mask:
            (self.vol1,p) = normaliseUnderMask( vol1, self.mask)
        else:
            self.vol1 = mean0std1( vol1, True)


    def evalScore(self, rot_trans):
        """
        evaluate score for given rotation and translation - NEGATIVE cc because\
        optimizer MINIMIZES

        @param rot_trans: rotation and translation vector (6-dim: z1,z2,x \
            angles, x,y,z, translations)
        @type rot_trans: L{list}
        @author: FF
        """
        from pytom.lib.pytom_volume import transformSpline, transform
        if self.interpolation.lower() == 'spline':
            transformSpline(self.vol2, self.rotvol2, rot_trans[0], rot_trans[1], rot_trans[2],
                            self.centX,self.centY,self.centZ,0,0,0,
                            rot_trans[3]/self.binning, rot_trans[4]/self.binning, rot_trans[5]/self.binning)
        else:
            transform(self.vol2, self.rotvol2, rot_trans[0], rot_trans[1], rot_trans[2],
                      self.centX,self.centY,self.centZ,0,0,0,
                      rot_trans[3]/self.binning,rot_trans[4]/self.binning,rot_trans[5]/self.binning)
        self.val = -1.*(self.score(volume=self.vol1, template=self.rotvol2, mask=self.mask, volume_is_normalized=True))
        return self.val

    def cc(self, rot, trans):
        """
        evaluate the CC for given rotation and translation - here directly CC not multiplied by -1 as in evalScore

        @param rot: rotation of vol2
        @type rot: L{pytom.basic.Rotation}
        @param trans: translation of vol2
        @type trans: L{pytom.basic.Shift}
        @return: cc
        @rtype: L{float}
        @author: FF
        """
        from pytom.lib.pytom_volume import transformSpline

        #transform vol2
        transformSpline( self.vol2, self.rotvol2, rot[0], rot[1], rot[2],
            self.centX,self.centY,self.centZ,0,0,0,
            trans[0]/self.binning, trans[1]/self.binning, trans[2]/self.binning)
        # compute CCC
        cc = self.score(volume=self.vol1, template=self.rotvol2,
                        mask=self.mask, volume_is_normalized=True)
        return cc

    def localOpti( self, iniRot=None, iniTrans=None):
        """
        @param iniRot: initial rotation of vol2
        @type iniRot: L{pytom.basic.Rotation}
        @param iniTrans: initial translation of vol2
        @type iniTrans: L{pytom.basic.Shift}
        @return: opti_score, opti_rot, opti_trans
        @rtype: L{float}, L{pytom.basic.structures.Rotation}, L{pytom.basic.structures.Shift}
        @author: FF
        """
        from pytom.basic.structures import Rotation, Shift
    
        if not (type(iniRot) == type(None)):
            self.set_rot( rot=iniRot)
        if not (type(iniTrans) == type(None)):
            self.set_trans( trans=iniTrans)

        if self.verbose:
            # alignment score before optimization
            print("CC before optimization %1.3f" % (-1.*self.evalScore(self.rot_trans)))
    
        # optimize scoring function
        maxiter=20
        if ((self.opti == 'fmin') or (self.opti == 'fmin_powell')):
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                xtol=0.05, ftol=0.001, maxiter=maxiter, maxfun=maxiter*20)
                #xtol=0.0001, ftol=0.0001, maxiter=maxiter, maxfun=None)
        elif self.opti == 'fmin_cg':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans, 
                gtol=0.0000001,
                maxiter=maxiter)
        elif self.opti == 'fmin_slsqp':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans, 
                iter=maxiter, acc=1e-03)
        elif self.opti == 'fmin_bfgs':
            rot_trans = self.optimizer(self.evalScore, self.rot_trans,
                maxiter=maxiter, epsilon=1e-06)
        elif self.opti == 'leastsq':
            rot_trans, success = self.optimizer(self.evalScore, self.rot_trans, 
                maxfev=maxiter, epsfcn=0.0, factor=10)
        self.rot_trans = rot_trans
    
        # alignment score before optimization
        finscore = self.evalScore(self.rot_trans)
        rot, trans = self.vector2transRot(rot_trans)
        if self.verbose:
            print("CC after optimization %1.3f" % (-1.*finscore))
            print("rot_trans = ", rot_trans)
            print(rot, trans)

        return -1.*finscore, rot, trans
    

def alignVolumesAndFilterByFSC(vol1, vol2, mask=None, nband=None, iniRot=None, iniTrans=None, interpolation='linear',
                               fsc_criterion=0.143, verbose=0):
    """
    align two volumes, compute their FSC, and filter by FSC
    @param vol1: volume 1
    @param vol2: volume 2
    @mask: mask volume
    @type mask: L{pytom.lib.pytom_volume.vol}
    @param nband: Number of bands
    @type nband: L{int}
    @param iniRot: initial guess for rotation
    @param iniTrans: initial guess for translation
    @param interpolation: interpolation type - 'linear' (default) or 'spline'
    @param fsc_criterion: filter -> 0 according to resolution criterion
    @type fsc_criterion: float
    @param verbose: verbose level (0=mute, 1 some output, 2=talkative)
    @type verbose: int
    @type interpolation: str
    @return: (filvol1, filvol2, fsc, fsc_fil, optiRot, optiTrans) i.e., filtered volumes, their FSC, the corresponding\
        filter that was applied to the volumes, and the optimal rotation and translation of vol2 with respect to vol1\
        note: filvol2 is NOT rotated and translated!
    @author: FF
    """
    from pytom.lib.pytom_volume import transformSpline, vol
    from pytom.basic.correlation import fsc
    from pytom.basic.filter import filter_volume_by_profile
    from pytom.alignment.localOptimization import Alignment
    from pytom.basic.correlation import nxcc

    assert isinstance(vol1, vol), "alignVolumesAndFilterByFSC: vol1 must be of type vol"
    assert isinstance(vol2, vol), "alignVolumesAndFilterByFSC: vol2 must be of type vol"
    # filter volumes prior to alignment according to SNR
    fsc_prior = fsc(volume1=vol1, volume2=vol2, number_bands=nband)
    fil = design_fsc_filter(fsc=fsc_prior, fildim=int(vol2.sizeX()//2))
    #filter only one volume so that resulting CCC is weighted by SNR only once
    filvol2 = filter_volume_by_profile(volume=vol2, profile=fil)
    # align vol2 to vol1
    if verbose == 2:
        alignment = Alignment(vol1=vol1, vol2=filvol2, score=nxcc, mask=mask,
                              iniRot=iniRot, iniTrans=iniTrans, opti='fmin_powell', interpolation=interpolation,
                              verbose=verbose)
    else:
        alignment = Alignment(vol1=vol1, vol2=filvol2, score=nxcc, mask=mask,
                              iniRot=iniRot, iniTrans=iniTrans, opti='fmin_powell', interpolation=interpolation,
                              verbose=False)
    optiScore, optiRot, optiTrans = alignment.localOpti( iniRot=iniRot, iniTrans=iniTrans)
    if verbose:
        from pytom.angles.angleFnc import differenceAngleOfTwoRotations
        from pytom.basic.structures import Rotation
        diffAng = differenceAngleOfTwoRotations(rotation1=Rotation(0,0,0), rotation2=optiRot)
        print("Alignment densities: Rotations: %2.3f, %2.3f, %2.3f; Translations: %2.3f, %2.3f, %2.3f " % (optiRot[0],
                                    optiRot[1], optiRot[2], optiTrans[0], optiTrans[1], optiTrans[2]))
        print("Orientation difference: %2.3f deg" % diffAng)
    vol2_alig = vol(vol2.sizeX(), vol2.sizeY(), vol2.sizeZ())
    transformSpline(vol2, vol2_alig, optiRot[0], optiRot[1], optiRot[2],
                    int(vol2.sizeX()/2),int(vol2.sizeY()/2),int(vol2.sizeY()/2),
                    0, 0, 0,
                    optiTrans[0], optiTrans[1], optiTrans[2])
    # finally compute FSC and filter of both volumes
    if not nband:
        nband = int(vol2.sizeX()/2)
    fsc_out = fsc(volume1=vol1, volume2=vol2_alig, number_bands=nband)
    fil = design_fsc_filter(fsc=fsc_out, fildim=int(vol2.sizeX()/2), fsc_criterion=fsc_criterion)
    filvol1 = filter_volume_by_profile(volume=vol1, profile=fil)
    #filvol2 = filter_volume_by_profile( volume=vol2_alig, profile=fil)
    filvol2 = filter_volume_by_profile(volume=vol2, profile=fil)

    return (filvol1, filvol2, fsc_out, fil, optiRot, optiTrans)


def design_fsc_filter(fsc, fildim=None, fsc_criterion=0.143):
    """
    design spectral filter to weight by SNR of frequency
    @param fsc: input fsc
    @type fsc: 1-d list
    @param fildim: filter dimension
    @type fildim: int
    @return: filter
    @rtype: list
    @author: FF
    """
    from math import sqrt
    from pytom.basic.resolution import getResolutionBandFromFSC
    if not fildim:
        fildim = len(fsc)
    nband = len(fsc)
    if fsc_criterion != 0.0:
        resolutionBand = getResolutionBandFromFSC(fsc, criterion=fsc_criterion)
        smooth = max(resolutionBand/5,2)
    else:
        resolutionBand = len(fsc)
        smooth = 1
    #print "filter: fsc_criterion %2.3f, resolutionBand %d" % (fsc_criterion, resolutionBand)
    # filter by sqrt(FSC)
    fsc_fil = len(fsc)*[0.]
    for ii in range(0,len(fsc)):
        fscval = fsc[ii]
        if fscval > 0.:
            #fsc_fil[ii] = sqrt(fsc[ii])
            if ii <= resolutionBand:
                fsc_fil[ii] = sqrt(fsc[ii])
            elif (ii > resolutionBand) and (ii <= resolutionBand + smooth):
                fsc_fil[ii] = sqrt(fsc[ii]) * (((resolutionBand + smooth) - ii)/smooth)**2  # squared filter
            else:
                fsc_fil[ii] = 0.
        else:
            fsc_fil[ii] = 0.
    #design filter
    fil = fildim*[0.]
    if nband != len(fil):
        shrinkfac = 1./float(len(fil)/nband)
    for ii in range(len(fil)-1):
        # linear resample fsc if nband ~= size(fil)
        if nband != len(fil):
            ilow = int(shrinkfac*ii)
            dlow = shrinkfac*ii - ilow
            ihi  = ilow+1
            dhi  = ihi - shrinkfac*ii
            fil[ii] = fsc_fil[ilow]*dhi + fsc_fil[ihi]*dlow
        else:
            fil[ii] = fsc_fil[ii]
    return fil
