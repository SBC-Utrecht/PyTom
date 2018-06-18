'''
Created on Aug 6, 2010

@author: chen
'''

class ExPeakResult:
    """
    ExPeakResult: ExPeakResult of scoring function, used to retrieve score and orientation values.
    """
    
    def __init__(self, volFilename='', resultFilename='', orientFilename='', angleListFilename='', score=None):
        """
        @param volFilename: filename of original volume
        @type volFilename: string
        @param resultFilename: filename of result volume
        @type resultFilename: string
        @param orientFilename: filename of orientation volume
        @type orientFilename: string
        @param angleListFilename: filename of used angle list
        @type angleListFilename: string
        @param score: score type of getting this result
        @type score: L{pytom.alignment.score}
        
        @author: chen
        """
        self.volFilename = volFilename
        self.resultFilename = resultFilename
        self.orientFilename = orientFilename
        self.angleListFilename = angleListFilename
        self.score = score
        
        self.result = None
        self.orient = None
        self.angleList = None
        
        
    def readAll(self):
        """
        readAll: Read all the specified files from the disk (except the original volume)
        """
        from pytom.localization.structures import Volume, Orientation
        self.result = Volume(self.resultFilename).getVolume()
        self.orient = Orientation(self.orientFilename).getVolume()
        
        if self.angleListFilename!='':
            from pytom.angles.globalSampling import GlobalSampling
            self.angleList = GlobalSampling(self.angleListFilename).getRotations()
        
    def get(self, x, y, z):
        """
        get : Get both the score and orientation values according to index
        
        @param x: x dimension
        @type x: integer
        @param y: y dimension
        @type y: integer
        @param z: z dimension
        @type z: integer
        @return: both score and orientation
        @rtype: [float, [float, float, float] ]
        
        @author: chen
        """
        if self.result==None or self.orient==None or self.angleList==None:
            self.readAll()
        
        score = self.result.getV(x, y, z)
        orientIdx = self.orient.getV(x, y, z)
        if int(orientIdx)<len(self.angleList):
            orientation = self.angleList[int(orientIdx)]
        else:
            orientation = self.angleList[-1]
        
        return [score, orientation]
    
    def maskOut(self, mask, center, size):
        """
        maskOut: Set part of mask volume to all zero. The region is specified by center and size.
        @param mask: volume that you handle with
        @type mask: L{pytom_volume.vol}
        @param center: center of the region
        @type center: [x,y,z]
        @param size: size of the region
        @type size: [sizeX, sizeY, sizeZ] or radius
        """
        
        from pytom_volume import vol, putSubVolume
            
        if size.__class__ == list:
            p_sizeX = size[0]
            p_sizeY = size[1]
            p_sizeZ = size[2]
        elif size.__class__ == vol:
            mm = size
            p_sizeX = mm.sizeX()
            p_sizeY = mm.sizeY()
            p_sizeZ = mm.sizeZ()
        else:
            radius = size
            p_sizeX = radius*2
            p_sizeY = radius*2
            p_sizeZ = radius*2
            
        
        maskSize = [mask.sizeX(), mask.sizeY(), mask.sizeZ()]
        
        if maskSize < center:
            raise RuntimeError('Center out of range!')
        
        # [)
        # mask out double size. CHANGED!!!
        startX = int(center[0]-p_sizeX/2)
        endX = int(center[0]+p_sizeX/2)
        startY = int(center[1]-p_sizeY/2)
        endY = int(center[1]+p_sizeY/2)
        startZ = int(center[2]-p_sizeZ/2)
        endZ = int(center[2]+p_sizeZ/2)
        
        # only used for radius
        sub_startX = 0; sub_startY = 0; sub_startZ = 0
        
        if startX < 0:
            sub_startX = -startX
            startX = 0
        if endX > maskSize[0]:
            endX = maskSize[0]
        if startY < 0:
            sub_startY = -startY
            startY = 0
        if endY > maskSize[1]:
            endY = maskSize[1]
        if startZ < 0:
            sub_startZ = -startZ
            startZ = 0
        if endZ > maskSize[2]:
            endZ = maskSize[2]
            
        sizeX = endX - startX
        sizeY = endY - startY
        sizeZ = endZ - startZ
        
        if size.__class__ == list:
            subV = vol(sizeX, sizeY, sizeZ)
            subV.setAll(0)
        elif size.__class__ == vol:
            from pytom_volume import limit, subvolume
            subV = (mm-1)/-1
            limit(subV, 0.999, 0, 0,0, True, False)
            subV = subvolume(subV, sub_startX, sub_startY, sub_startZ, sizeX, sizeY, sizeZ)
            tempV = subvolume(mask, startX, startY, startZ, sizeX, sizeY, sizeZ)
            subV = subV*tempV # AND operation
        else:
            from pytom_volume import initSphere, subvolume
            subV = vol(radius*2, radius*2, radius*2)
            initSphere(subV, radius, 0, 0, radius, radius, radius)
            tempV = vol(radius*2, radius*2, radius*2)
            tempV.setAll(1)
            subV = tempV - subV
            subV = subvolume(subV, sub_startX, sub_startY, sub_startZ, sizeX, sizeY, sizeZ)
            tempV = subvolume(mask, startX, startY, startZ, sizeX, sizeY, sizeZ)
            subV = subV*tempV # AND operation
            
        putSubVolume(subV, mask, startX, startY, startZ)
        
    
    def findParticles(self, sizeParticle, maxNumParticle=0, minScore=-1, write2disk=0, margin=None):
        """
        findParticles: Find particles in target volume according to the result volume.
        @param sizeParticle: size or radius of searched particle
        @type sizeParticle: [x,y,z] or integer
        @param maxNumParticle: maximal number of particles you want to pick
        @type maxNumParticle: integer
        @param minScore: minimal score as threshold
        @type minScore: float 
        @param write2disk: write the found particles to the disk or not (0: do not write, otherwise the length of each dimension)
        @type write2disk: integer
        @param margin: set the margin of the score volume
        @param margin: [x,y,z] or integer
        
        @return: list of found particles
        @rtype: L{pytom.localization.structures.FoundParticle}
        """
        from pytom_volume import vol, peak, putSubVolume, read
        from pytom.localization.structures import FoundParticle
        
        # prepare the mask
        x = self.result.sizeX(); y = self.result.sizeY(); z = self.result.sizeZ()
        
        if sizeParticle.__class__ == list:
            xP = sizeParticle[0]; yP = sizeParticle[1]; zP = sizeParticle[2]
        elif sizeParticle.__class__ == vol:
            xP = sizeParticle.sizeX(); yP = sizeParticle.sizeY(); zP = sizeParticle.sizeZ()
        else:
            radius = sizeParticle
            xP = 2*sizeParticle; yP = 2*sizeParticle; zP = 2*sizeParticle
        
        if margin:
            if margin.__class__ == list:
                marginX, marginY, marginZ = margin
            else:
                marginX = marginY = marginZ = margin
        else: # no margin given, set automatically
            marginX = int(xP/2); marginY = int(yP/2); marginZ = int(zP/2)
        
        mask = vol(x,y,z)
        mask.setAll(0)
        
        maskIn = vol(x-2*marginX, y-2*marginY, z-2*marginZ)
        maskIn.setAll(1)
        putSubVolume(maskIn, mask, marginX, marginY, marginZ)
        
        # progress bar
        from pytom.tools.ProgressBar import FixedProgBar
        prog = FixedProgBar(0, maxNumParticle-1, '')
        
        # find the particles
        resList = []
        for i in xrange(maxNumParticle):
            prog.update(i)
            
            try:
                posV = peak(self.result, mask)
            except:
                break # the mask is all zero
            
            [scoreV, orientV]= self.get(posV[0],posV[1],posV[2])
            # test if the peak score is bigger than minimal threshold
            if scoreV > minScore:
                particleFilename = 'particle_'+str(i)+'.em'
                if write2disk:
                    # write the found particle to the disk
                    l = write2disk
                    v = read(self.volFilename, posV[0]-l/2, posV[1]-l/2, posV[2]-l/2, l, l, l, 0,0,0,0,0,0)
                    v.write(particleFilename)
                
                score = self.score()
                score.setValue(scoreV)
                from pytom.basic.structures import PickPosition, Rotation
                pos = PickPosition(posV, originFilename=self.volFilename)
                orientation = Rotation(orientV)
                p = FoundParticle(pos, orientation, score, particleFilename)
                
                resList.append(p)
                
                if sizeParticle.__class__ == list:
                    self.maskOut(mask, posV, [xP, yP, zP])
                elif sizeParticle.__class__ == vol:
                    self.maskOut(mask, posV, sizeParticle)
                else:
                    self.maskOut(mask, posV, radius)
            else:
                break
            
        return resList
