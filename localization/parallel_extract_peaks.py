'''
Created on May 20, 2010

@author: chen
'''
import numpy

def getMsgStr():
    '''
    getMsgStr: Use MPI to receive the message
    @return: message string
    @rtype: str
    '''
    import pytom_mpi
    mpi_msgString = pytom_mpi.receive()
    return mpi_msgString

class PeakWorker(object):
    '''
    PeakWorker: Worker to calculate the score function
    @author: Chen
    '''
    def __init__(self):
        import pytom_mpi
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
            
        self.mpi_id = pytom_mpi.rank()
        self.name = 'node_' + str(self.mpi_id)
        self.runtimes = 0
        
    def getJobMsg(self, mpi_msgString):
        '''
        getJobMsg: Get the job message from received string
        @param mpi_msgString: message string
        @type mpi_msgString: string
        @return: message
        @rtype: L{pytom.localization.peak_job_msg.PeakJobMsg}
        '''
        from pytom.localization.peak_job_msg import PeakJobMsg
        
        msg = PeakJobMsg()
        msg.fromStr(mpi_msgString)
        
        return msg
    
    def setJob(self, job):
        self.volume = job.volume
        self.reference = job.reference
        self.mask = job.mask
        self.rotations = job.rotations
        self.wedge = job.wedge
        self.score = job.score
        self.jobID = job.jobID
        self.bandpass = job.bandpass
        
    def jobFromMsg(self, msg):
        """
        jobFromMsg: Set the job information from received message
        @param msg: received message
        @type msg: L{pytom.localization.peak_job_msg.PeakJobMsg}
        """
        from pytom.localization.peak_job_msg import PeakJobMsg
        from pytom.localization.peak_job import PeakJob
        job = msg.getJob()
        self.setJob(job)
        if msg.getSender() != '':
            self.backTo = int(msg.getSender())

    
    def run(self, verbose=True, moreInfo=False, gpuID=-1):
        """
        run: Run the worker and return the result
        @param verbose: verbose mode
        @type verbose: boolean
        @return: result of calculation
        @rtype: L{pytom.localization.peak_job.PeakResult}
        """
#        from pytom.tools.timing import timing
#        t = timing(); t.start()
        
        # read the necessary files on the disk

        v = self.volume.getVolume(self.volume.subregion)
        ref = self.reference.getVolume()
        maskFilename = self.mask.getFilename()
        if maskFilename == '': # No mask is given, use the default
            m = None
            mIsSphere = True
        else:
            m = self.mask.getVolume()
            mIsSphere = self.mask.isSphere()
        rot = self.rotations
        wedg = self.wedge
        scoreFnc = self.score.getScoreFunc()
        
        # apply the bandpass
        if hasattr(self, 'bandpass') and self.bandpass:
            if verbose: print('Bandpass added!')
            ref = self.bandpass.filter(ref)
        
        # calculate the result volume
        if gpuID is None:
            from pytom.localization.extractPeaks import extractPeaks
        else:
            from pytom.localization.extractPeaks import extractPeaksGPU as extractPeaks

        if verbose==True:
            print(self.name + ': starting to calculate %d rotations' % rot.numberRotations())
        [resV, orientV, sumV, sqrV] = extractPeaks(v, ref, rot, scoreFnc, m, mIsSphere, wedg, nodeName=self.name,
                                                   verboseMode=verbose, moreInfo=moreInfo, gpuID=gpuID)
        
        self.runtimes = self.runtimes + 1
        
#        time = t.end(); print '%s calculated in: %f' % (self.name, time)
        return [resV, orientV, sumV, sqrV]
    
    def parallelRun(self, verbose=True, **kwargs):
        '''
        parallelRun: Run the worker in parallel status and send the result message back
        @param verbose: verbose mode
        @type verbose: boolean
        '''
        
        from pytom.parallel.messages import StatusMessage,MessageError
        from pytom.basic.exceptions import ParameterError
        from pytom.basic.structures import PyTomClassError
        
        
        end = False
        
        while not end:
            # get the message string
            mpi_msgString = getMsgStr()
            
            try:
                msg = self.getJobMsg(mpi_msgString)
                self.jobFromMsg(msg)
                
                if verbose==True:
                    print(self.name + ': running...')
                [resV, orientV, sumV, sqrV] = self.run(verbose, moreInfo=True)
                
                # write the result volume back to the disk
                resFilename = self.name + '_job' + str(self.jobID) + '_res.em'
                orientFilename = self.name + '_job' + str(self.jobID) + '_orient.em'
                resV.write(resFilename)
                orientV.write(orientFilename)
                
                if sumV and sqrV:
                    sumFilename = self.name + '_job' + str(self.jobID) + '_sum.em'
                    sqrFilename = self.name + '_job' + str(self.jobID) + '_sqr.em'
                    sumV.write(sumFilename)
                    sqrV.write(sqrFilename)
                
                from pytom.localization.structures import Volume, Orientation
                res = Volume(resFilename, self.volume.subregion)
                orient = Orientation(orientFilename)
                
                # construct the result
                from pytom.localization.peak_job import PeakResult
                result = PeakResult(res, orient, self.jobID)
                
                # sent the result back
                if verbose==True:
                    print(self.name + ': sending back result')
                result.send(self.mpi_id, self.backTo)
                
            except (MessageError,PyTomClassError,ParameterError):
                try:
                    # get the message as StatusMessage and finish
                    msg = StatusMessage('','')
                    msg.fromStr(mpi_msgString)
                    if msg.getStatus() == 'End':
                        end = True
                    if verbose==True:
                        print(self.name + ': ending...')
                except (MessageError,PyTomClassError,ParameterError):
                    print('Error parsing message. Message either unknown or invalid.')
                    assert False
            
        
class PeakManager():
    '''
    PeakManager: Manager to manage the workers
    '''
    def __init__(self, suffix=''):
        self.name = 'node_0'
        self.numWorkers = 0
        self.jobInfo = {} # store the information of distributed jobs
        self.suffix = suffix

    def sequentialStart(self, job):
        '''
        sequentialStart: Start the worker to sequentially finish the job
        @param job: job to be done
        @type job: L{pytom.localization.peak_job.PeakJob}
        '''
        # sequentially start the job without using MPI
        from pytom.localization.peak_job_msg import PeakJobMsg
        jobMsg = PeakJobMsg()
        jobMsg.setJob(job)
        
        worker = PeakWorker()
        worker.jobFromMsg(jobMsg)
        result = worker.run()
        
        result[0].write("res.em")
        result[1].write("orient.em")
            
    def getResMsg(self, mpi_msgString):
        '''
        getResMsg: Get the result message from received string
        @param mpi_msgString: message string
        @type mpi_msgString: string
        @return: message
        @rtype: L{pytom.localization.peak_job_msg.PeakResultMsg}
        '''
        from pytom.localization.peak_job_msg import PeakResultMsg
        from pytom.tools.files import dumpMsg2Log
        
        logFileName = self.name + '.log'
        dumpMsg2Log(logFileName,mpi_msgString)
        
        msg = PeakResultMsg()
        msg.fromStr(mpi_msgString)
        
        return msg
    
    def resFromMsg(self, msg):
        return msg.getResult()
    
#    def procResFromWorker(self, v, orient, resV, resO):
#        """
#        procResFromWorker: process the results from the workers
#        @param v: to be updated result volume
#        @type v: L{pytom_volume.vol}
#        @param orient: to be updated orientation volume
#        @type orient: L{pytom_volume.vol}
#        @param resV: result volume from the worker
#        @type resV: L{pytom_volume.vol}
#        @param resO: result orientation volume from the worker
#        @type resO: L{pytom_volume.vol}
#        """
#        for i in xrange(v.sizeX()):
#            for j in xrange(v.sizeY()):
#                for k in xrange(v.sizeZ()):
#                    vT = v.getV(i,j,k); vN = resV.getV(i,j,k)
#                    if vT < vN:
#                        v.setV(vN, i,j,k)
#                        orient.setV(resO.getV(i,j,k), i,j,k)
        
    def parallelInit(self):
        '''
        parallelInit: Initialization for the parallelization
        '''
        import pytom_mpi
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
            
        if pytom_mpi.size() < 2:
            raise RuntimeError('Number of available cluster nodes is less than 2.')
        
        self.numWorkers = pytom_mpi.size()-1

    def splitAngles(self, job, verbose=True):
        '''
        splitAngles: Distribute the job to the workers by splitting angles
        @param job: job to be done
        @type job: L{pytom.localization.peak_job.PeakJob}
        @param verbose: verbose mode
        @type verbose: boolean
        '''
        import pytom_mpi
        mpi_myid = pytom_mpi.rank()
        
        if not mpi_myid == 0:
            raise RuntimeError('This function can only be processed by mpi_id = 0! ID == ' + str(mpi_myid) +' Aborting!')
        
        # split the job into smaller ones 
        rotationsPerWorker = job.rotations.numberRotations()//self.numWorkers
        if rotationsPerWorker == 0:
            raise RuntimeError("Not enough angles to split!")
        
        if verbose==True:
            print('\n\nManager: distribute number of %d rotations to %d workers' % (job.rotations.numberRotations(), self.numWorkers))
        
        for i in range(1, self.numWorkers+1):
            # split the rotations
            if i != self.numWorkers:
                subRot = job.rotations[(i-1)*rotationsPerWorker : i*rotationsPerWorker]
            else: # the last node will take all the rest of rotations
                subRot = job.rotations[(i-1)*rotationsPerWorker : ]
                
            self.jobInfo[i] = (i-1)*rotationsPerWorker
            
            from pytom.angles.angleList import AngleList
            rot = AngleList(subRot)
            
            from pytom.localization.peak_job import PeakJob
            subJob = PeakJob(volume = job.volume, reference = job.reference, 
                            mask = job.mask, wedge = job.wedge, rotations = rot, 
                            score = job.score, jobID = i, 
                            dstDir = job.dstDir, bandpass = job.bandpass)
            subJob.send(0, i)
    
    def parallelStart_splitAng(self, job, verbose=True):
        '''
        parallelStart_splitAng: Start the parallel working for the job in the way of splitting angles
        @param job: job to be done
        @type job: L{pytom.localization.peak_job.PeakJob}
        @param verbose: verbose mode
        @type verbose: boolean
        '''
        
        self.parallelInit()
        
        import pytom_mpi
        import os
        mpi_myid = pytom_mpi.rank()
        if mpi_myid == 0: # manager
            # distribute the job to the workers
            self.splitAngles(job, verbose)
            
            # gather results and post processing
            i = 0
            volume = None
            orient = None
            while i < self.numWorkers:
                mpi_msgString = getMsgStr()
                msg = self.getResMsg(mpi_msgString)
                resFromWorker = self.resFromMsg(msg)
                
                if verbose == True:
                    print("Manager: processing result from worker " + msg.getSender())
                    
                resV = resFromWorker.result.getVolume()
                resO = resFromWorker.orient.getVolume()
                jobID = resFromWorker.jobID
                
                # correct the orientation information
#                n = int(msg.getSender())
#                rotationsPerWorker = job.rotations.numberRotations()/self.numWorkers
#                offset = (n-1)*rotationsPerWorker
#                print "Manager: Offset is %d and %d" % (offset, self.jobInfo[jobID])
                offset = self.jobInfo[jobID]
                resO = resO + offset
                
                if volume==None or orient==None:
                    from pytom_volume import vol
                    volume = vol(resV.sizeX(), resV.sizeY(), resV.sizeZ())
                    orient = vol(resO.sizeX(), resO.sizeY(), resO.sizeZ())
                    volume.copyVolume(resV)
                    orient.copyVolume(resO)
                else:
#                    self.procResFromWorker(volume, orient, resV, resO)
                    from pytom_volume import updateResFromVol
                    updateResFromVol(volume, resV, orient, resO)
                    
                i = i+1
            
            # write the final result back to the disk
            volume.write(self.name + '_res.em')
            orient.write(self.name + '_orient.em')
            
            # post processing the sum and sqr volume
            print("Start post processing the sum and sqr volume ...")
            sumList = []
            sqrList = []
            filenames = os.listdir('.')
            for name in filenames:
                if name[-6:] == 'sum.em':
                    sumList.append(name)
                elif name[-6:] == 'sqr.em':
                    sqrList.append(name)
                else:
                    pass
            
            sumV = None
            sqrV = None
            from pytom_volume import read
            for name in sumList:
                v = read(name)
                if sumV:
                    sumV = sumV + v
                else:
                    sumV = v
            
            for name in sqrList:
                v = read(name)
                if sqrV:
                    sqrV = sqrV + v
                else:
                    sqrV = v
            
            sumV = sumV / job.rotations.numberRotations()
            sqrV = sqrV / job.rotations.numberRotations()
            
            sumV.write('node_0_sum.em')
            sqrV.write('node_0_sqr.em')
            
            # delete the temporary files on the disk
            import os
            files = os.listdir('.')
            for name in files:
                if 'job' in name and '.em' in name:
                    os.remove(name)
                    
            # rename the result files name
            os.rename('node_0_res.em', 'scores.em')
            os.rename('node_0_orient.em', 'angles.em')
            
            # finishing, stop all workers
            self.parallelEnd(verbose)
            
            if verbose == True:
                print("Manager: end")
            
        else: # worker
            worker = PeakWorker()
            worker.parallelRun(verbose)


    def splitVolumes(self, job, splitX, splitY, splitZ, verbose=True):
        """
        splitVolumes: Split the target volume in job into smaller ones.
        """
        # check if the split is feasible
        v = job.volume.getVolume()
        r = job.reference.getVolume()

        print(r.sizeX(), r.sizeY())
        
        vsizeX = v.sizeX(); vsizeY = v.sizeY(); vsizeZ = v.sizeZ()
        sizeX = vsizeX//splitX; sizeY = vsizeY//splitY; sizeZ = vsizeZ//splitZ
        rsizeX = r.sizeX(); rsizeY = r.sizeY(); rsizeZ = r.sizeZ()
        if rsizeX>sizeX or rsizeY>sizeY or rsizeZ>sizeZ:
            raise RuntimeError("Not big enough volume to split!")
        
        self.jobInfo["originalSize"] = [vsizeX, vsizeY, vsizeZ]
        self.jobInfo["splitSize"] = [sizeX, sizeY, sizeZ]
        
        # read the target volume, calculate the respective subregion
        from pytom.localization.peak_job import PeakJob
        from pytom.localization.structures import Volume
        _start = [-rsizeX//2,-rsizeY//2,-rsizeZ//2]
        _size = [sizeX+rsizeX, sizeY+rsizeY, sizeZ+rsizeZ]
        
        for i in range(splitX*splitY*splitZ):
            strideZ = splitX*splitY; strideY = splitX
            incZ = i//strideZ; incY = (i%strideZ)//strideY; incX = i%strideY
            _start = [-rsizeX//2+incX*sizeX,-rsizeY//2+incY*sizeY,-rsizeZ//2+incZ*sizeZ]
            
            start = _start[:]
            end = [start[j]+_size[j] for j in range(len(start))]
            
            if start[0] < 0:
                start[0]=0
            if start[1] < 0:
                start[1]=0
            if start[2] < 0:
                start[2]=0
            if end[0] > vsizeX:
                end[0] = vsizeX
            if end[1] > vsizeY:
                end[1] = vsizeY
            if end[2] > vsizeZ:
                end[2] = vsizeZ
            
            size = [end[j]-start[j] for j in range(len(start))]
#            print start[0], start[1], start[2], size[0], size[1], size[2]
            print(size)
            # for reassembling the result
            whole_start = start[:]
            sub_start = [0, 0, 0]
            if start[0] != 0:
                whole_start[0] = start[0]+rsizeX//2
                sub_start[0] = rsizeX//2
            if start[1] != 0:
                whole_start[1] = start[1]+rsizeY//2
                sub_start[1] = rsizeY//2
            if start[2] != 0:
                whole_start[2] = start[2]+rsizeZ//2
                sub_start[2] = rsizeZ//2
            self.jobInfo[i+1] = [sub_start, whole_start]
            
            subVol = Volume(job.volume.getFilename(),
                            [start[0], start[1], start[2], size[0], size[1], size[2]])
            subJob = PeakJob(subVol, job.reference, job.mask, job.wedge, job.rotations, job.score, i+1, job.dstDir, job.bandpass)
            subJob.send(0, i%self.numWorkers+1)
        
    def parallelStart_splitVol(self, job, splitX, splitY, splitZ, verbose=True):
        """
        """
        
        self.parallelInit()
        
        import pytom_mpi
        mpi_myid = pytom_mpi.rank()
        if mpi_myid == 0: # manager
            # distribute the job to the workers and get the original size of the volume
            self.splitVolumes(job, splitX, splitY, splitZ, verbose)

            [vsizeX,vsizeY,vsizeZ] = self.jobInfo["originalSize"]
            
            # gather results and post processing
            from pytom_volume import vol
            volume = vol(vsizeX, vsizeY, vsizeZ)
            volume.setAll(0)
            orient = vol(vsizeX, vsizeY, vsizeZ)
            orient.setAll(0)
            
            i = 0
            while i < splitX*splitY*splitZ:
                mpi_msgString = getMsgStr()
                msg = self.getResMsg(mpi_msgString)
                resFromWorker = self.resFromMsg(msg)
                
                if verbose == True:
                    print("Manager: processing result from worker " + msg.getSender())
                    
                resV = resFromWorker.result.getVolume()
                resO = resFromWorker.orient.getVolume()
                jobID = resFromWorker.jobID
                
                [sizeX ,sizeY, sizeZ] = self.jobInfo["splitSize"]
                
                [sub_start, start] = self.jobInfo[jobID]
                
                from pytom_volume import subvolume, putSubVolume
                sub_resV = subvolume(resV, sub_start[0],sub_start[1],sub_start[2], sizeX,sizeY,sizeZ)
                sub_resO = subvolume(resO, sub_start[0],sub_start[1],sub_start[2], sizeX,sizeY,sizeZ)
                
                putSubVolume(sub_resV, volume, start[0],start[1],start[2])
                putSubVolume(sub_resO, orient, start[0],start[1],start[2])
                              
                i = i+1
            
            # write the final result back to the disk
            volume.write(self.name + '_res.em')
            orient.write(self.name + '_orient.em')
            
            # delete the temporary files on the disk
            import os
            files = os.listdir('.')
            for name in files:
                if 'job' in name and '.em' in name and not 'sum' in name and not 'sqr' in name:
                    os.remove(name)
            
            # rename the result files name
            os.rename('node_0_res.em', 'scores{}.em'.format(suffix))
            os.rename('node_0_orient.em', 'angles{}.em'.format(suffix))
            
            # finishing, stop all workers
            self.parallelEnd(verbose)
            
            if verbose == True:
                print("Manager: end")
            
        else: # worker
            worker = PeakWorker()
            worker.parallelRun(verbose)
        
        pytom_mpi.finalise()
        
#    def distributeJobs(self, job, splitX, splitY, splitZ, verbose=True):
#        pass
#    
#    def parallelStart(self, job, splitX, splitY, splitZ, verbose=True):
#        pass

    def parallelEnd(self, verbose=True):
        """
        parallelEnd : Sends status message = end to all workers.
        @param verbose: verbose mode
        @type verbose: boolean
        """
        
        if verbose == True:
            print('Manager: sending end messages to workers')
        
        import pytom_mpi
        from pytom.parallel.messages import StatusMessage
        
        mpi_numberNodes = pytom_mpi.size()
        mpi_myid = pytom_mpi.rank()
        
        for i in range(1,mpi_numberNodes):
            msg = StatusMessage(str(mpi_myid),str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg),i)


class PeakLeader(PeakWorker):
    """
    PeakLeader: Class for parallel running of jobs (new architecture)
    """
    def __init__(self,suffix=''):
        import pytom_mpi
        
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
        self.suffix=suffix
        self.mpi_id = pytom_mpi.rank()
        self.name = 'node_' + str(self.mpi_id)
        self.clean()
        
    def clean(self):
        """
        clean: Clean itself for new jobs
        """
        self.members = 1 # including itself
        self.backTo = None
        self.jobInfoPool = {}
        self.runtimes = 0
        self.dstDir = './'
        
        # result stored in itself
        self.resVolA = None
        self.resOrientA = None
        self.resVol = None
        self.resOrient = None
        
    def getMsgType(self, mpi_msgString):
        """
        getMsgType: Determine the message type
        
        @rtype: 2-PeakJobMsg, 1-PeakResultMsg, 0-StatusMessage
        """
        from lxml import etree
        
        xmlObj = etree.fromstring(mpi_msgString)
        
        if len(xmlObj.xpath('/PeakJobMsg')) > 0:
            return 2
        elif len(xmlObj.xpath('/PeakResultMsg')) > 0:
            return 1
        elif len(xmlObj.xpath('/StatusMessage')) > 0:
            return 0
        else:
            return -1
        
    def setJob(self, job):
        """
        setJob: Set the relevant fields before calculation
        """
        self.volume = job.volume
        self.reference = job.reference
        self.mask = job.mask
        self.rotations = job.rotations
        self.wedge = job.wedge
        self.score = job.score
        
        self.jobID = job.jobID
        self.members = job.members
        self.dstDir = job.dstDir
        
    def jobFromMsg(self, msg):
        """
        jobFromMsg: Set the job information from received message
        @param msg: received message
        @type msg: L{pytom.localization.peak_job_msg.PeakJobMsg}
        """
        from pytom.localization.peak_job_msg import PeakJobMsg
        from pytom.localization.peak_job import PeakJob
        job = msg.getJob()
        sendID = int(msg.getSender())
        if sendID != self.mpi_id: # make sure that the result will not send back to itself
            self.backTo = sendID
        
        return job

    def getResMsg(self, mpi_msgString):
        '''
        getResMsg: Get the result message from received string
        @param mpi_msgString: message string
        @type mpi_msgString: string
        @return: message
        @rtype: L{pytom.localization.peak_job_msg.PeakResultMsg}
        '''
        from pytom.localization.peak_job_msg import PeakResultMsg
        
        msg = PeakResultMsg()
        msg.fromStr(mpi_msgString)
        
        return msg
    
    def resFromMsg(self, msg):
        """
        resFromMsg: Get the result from result message
        
        @rtype: L{pytom.localization.peak_job.PeakResult}
        """
        return msg.getResult()

    def splitAngles(self, job, verbose=True):
        """
        splitAngles: Split the job in the "angle" way (binary split)
        @param job: job
        @type job: L{pytom.localization.peak_job.PeakJob}
        """
        from pytom.localization.peak_job import PeakJob
        
        totalNum = job.rotations.numberRotations()
        
#        # initialize the jobInfo structure
#        self.jobInfo["numDoneJobs"] = 0
#        self.jobInfo["numJobs"] = 0
#        self.jobInfo["originalJobID"] = job.jobID
#        self.jobInfo["splitType"] = "Ang"
        originalJobID = job.jobID
        
        while self.members > 1 and totalNum > 1:
            numEach = totalNum//self.members
            subMem1 = self.members//2
            subMem2 = self.members-self.members//2
            
            from pytom.angles.angleList import AngleList
            subRot1 = AngleList(job.rotations[:numEach*subMem1])
            subRot2 = AngleList(job.rotations[numEach*subMem1:])
            

            # avoid the collision of the job id
            subJob1 = PeakJob(job.volume, job.reference, job.mask, job.wedge, subRot1, job.score, job.jobID*10+1, subMem1, self.dstDir, job.bandpass)
            subJob2 = PeakJob(job.volume, job.reference, job.mask, job.wedge, subRot2, job.score, job.jobID*10+2, subMem2, self.dstDir, job.bandpass)
            
            if verbose==True:
                print(self.name+': send number of %d rotations to node %d' % (subJob2.rotations.numberRotations(), self.mpi_id+subMem1))
            subJob2.send(self.mpi_id, self.mpi_id+subMem1)
            
#            self.jobInfo["numJobs"] = self.jobInfo["numJobs"] + 1
#            # set offset
#            self.jobInfo[subJob1.jobID] = 0
#            self.jobInfo[subJob2.jobID] = numEach*subMem1
            self.jobInfoPool["numJobsA"] = self.jobInfoPool["numJobsA"] + 1
            from pytom.localization.peak_job import JobInfo
            info = JobInfo(subJob1.jobID, originalJobID, "Ang")
            info.angleOffset = 0
            self.jobInfoPool[subJob1.jobID] = info
            
            info = JobInfo(subJob2.jobID, originalJobID, "Ang")
            info.angleOffset = numEach*subMem1
            self.jobInfoPool[subJob2.jobID] = info
            
            job = subJob1
            totalNum = job.rotations.numberRotations()
            self.setJob(job)
            
#        self.jobInfo["numJobs"] = self.jobInfo["numJobs"] + 1
        self.jobInfoPool["numJobsA"] = self.jobInfoPool["numJobsA"] + 1
            

    def splitVolumes(self, job, splitX, splitY, splitZ, verbose=True):
        """
        splitVolumes: Split the job in the "volume" way (sequentially)
        @param job: job
        @type job: L{pytom.localization.peak_job.PeakJob}
        @param splitX: split part along the x dimension
        @type splitX: integer
        @param splitY: split part along the y dimension
        @type splitY: integer
        @param splitZ: split part along the z dimension
        @type splitZ: integer
        """
        # check if the split is feasible
        if job.volume.subregion == [0,0,0,0,0,0]:
            v = job.volume.getVolume()
            origin = [0,0,0]
            vsizeX = v.sizeX(); vsizeY = v.sizeY(); vsizeZ = v.sizeZ()
        else:
            origin = job.volume.subregion[0:3]
            vsizeX = job.volume.subregion[3]
            vsizeY = job.volume.subregion[4]
            vsizeZ = job.volume.subregion[5]
            
        sizeX = vsizeX//splitX; sizeY = vsizeY//splitY; sizeZ = vsizeZ//splitZ
        r = job.reference.getVolume()
        rsizeX = r.sizeX(); rsizeY = r.sizeY(); rsizeZ = r.sizeZ()
        if rsizeX>sizeX or rsizeY>sizeY or rsizeZ>sizeZ:
            raise RuntimeError("Not big enough volume to split!")
        
#        # initialize the jobInfo structure
#        self.jobInfo["numDoneJobs"] = 0
#        self.jobInfo["numJobs"] = 0
#        self.jobInfo["originalJobID"] = job.jobID
#        self.jobInfo["splitType"] = "Vol"
#
#        self.jobInfo["originalSize"] = [vsizeX, vsizeY, vsizeZ]
#        self.jobInfo["splitSize"] = [sizeX, sizeY, sizeZ]
        originalJobID = job.jobID
            
        # read the target volume, calculate the respective subregion
        from pytom.localization.peak_job import PeakJob
        from pytom.localization.structures import Volume
        _start = [-rsizeX//2+origin[0],-rsizeY//2+origin[1],-rsizeZ//2+origin[2]]
        _size = [sizeX+rsizeX, sizeY+rsizeY, sizeZ+rsizeZ]
        
        numPieces = splitX*splitY*splitZ
        totalMem = self.members
        numMemEach = totalMem//numPieces
        targetID = self.mpi_id
        
        for i in range(numPieces):
            strideZ = splitX*splitY; strideY = splitX
            incZ = i//strideZ; incY = (i%strideZ)//strideY; incX = i%strideY
            _start = [-rsizeX//2+origin[0]+incX*sizeX,-rsizeY//2+origin[1]+incY*sizeY,-rsizeZ//2+origin[2]+incZ*sizeZ]
            
            start = _start[:]
            end = [start[j]+_size[j] for j in range(len(start))]
            
            if start[0] < origin[0]:
                start[0]=origin[0]
            if start[1] < origin[1]:
                start[1]=origin[1]
            if start[2] < origin[2]:
                start[2]=origin[2]
            if end[0] > vsizeX+origin[0]:
                end[0] = vsizeX+origin[0]
            if end[1] > vsizeY+origin[1]:
                end[1] = vsizeY+origin[1]
            if end[2] > vsizeZ+origin[2]:
                end[2] = vsizeZ+origin[2]
            
            size = [end[j]-start[j] for j in range(len(start))]
            
#            # make sure that the last dimension is not odd
#            if size[2]%2 == 1:
#                size[2] = size[2] - 1
#                end[2] = end[2] - 1
            
            # for reassembling the result
            whole_start = start[:]
            sub_start = [0, 0, 0]
            if start[0] != origin[0]:
                whole_start[0] = start[0]+rsizeX//2
                sub_start[0] = rsizeX//2
            if start[1] != origin[1]:
                whole_start[1] = start[1]+rsizeY//2
                sub_start[1] = rsizeY//2
            if start[2] != origin[2]:
                whole_start[2] = start[2]+rsizeZ//2
                sub_start[2] = rsizeZ//2
            
#            self.jobInfo[subJobID] = [sub_start, whole_start]
            subJobID = job.jobID+i+1
            subVol = Volume(job.volume.getFilename(),
                            [start[0], start[1], start[2], size[0], size[1], size[2]])
            if i == 0:
                numMem = totalMem - (numPieces-1)*numMemEach
            else:
                numMem = numMemEach
            subJob = PeakJob(subVol, job.reference, job.mask, job.wedge, job.rotations, job.score, subJobID, numMem, self.dstDir, job.bandpass)
            
            from pytom.localization.peak_job import JobInfo
            info = JobInfo(subJob.jobID, originalJobID, "Vol")
            info.sub_start = sub_start
            info.whole_start = whole_start
            info.originalSize = [vsizeX, vsizeY, vsizeZ]
            info.splitSize = [sizeX, sizeY, sizeZ]
            info.origin = origin
            self.jobInfoPool[subJobID] = info
            
            if targetID == self.mpi_id:
                self.setJob(subJob)
                if self.members > 1:
                    self.splitAngles(subJob, verbose)
            else:
                if verbose==True:
                    print(self.name + ' : send part of the volume to ' + str(targetID))
                print(targetID, type(targetID))
                subJob.send(self.mpi_id, targetID)
            
            targetID = targetID + numMem
            self.jobInfoPool["numJobsV"] = self.jobInfoPool["numJobsV"] + 1
            
            
    def distributeJobs(self, job, splitX=0, splitY=0, splitZ=0):
        """
        distributeJobs: Distribute the jobs to the members (if any), as well as itself. \
        This function should be called only once for each node.\
        Depending on the members number, the node will decide how to split the job \
	automatically.

        @param job: job
        @type job: L{pytom.localization.peak_job.PeakJob}
        @param splitX: split part along the x dimension
        @type splitX: integer
        @param splitY: split part along the y dimension
        @type splitY: integer
        @param splitZ: split part along the z dimension
        @type splitZ: integer
        """
        
        numPieces = splitX*splitY*splitZ
        
        # see how many members do I have
        self.setJob(job)
        
        if self.members == 1:
            pass
        else:
            self.jobInfoPool["numDoneJobsA"] = 0
            self.jobInfoPool["numJobsA"] = 0
            self.jobInfoPool["numDoneJobsV"] = 0
            self.jobInfoPool["numJobsV"] = 0
            if numPieces==0 or numPieces==1 or numPieces > self.members: # node num not enough for split vol
                self.splitAngles(job)
            else:
                print(type(splitX),type(splitY),type(splitZ))
                self.splitVolumes(job, splitX, splitY, splitZ)
        
    
    def writeRes(self, resV, orientV, jobID=None):
        """
        writeRes: Write the result back to the disk, and return the PeakJobResult.
        @param resV: result volume
        @type resV: L{pytom_volume.vol}
        @param orientV: orientation volume
        @type orientV: L{pytom_volume.vol}
        @param jobID: ID of job
        @type jobID: integer
        
        @rtype: L{pytom.localization.peak_job.PeakResult}
        """
        if jobID != None:
            resFilename = self.dstDir + self.name + '_job' + str(jobID) + '_res.em'
            orientFilename = self.dstDir + self.name + '_job' + str(jobID) + '_orient.em'
        else:
            resFilename = self.dstDir + self.name + '_res.em'
            orientFilename = self.dstDir + self.name + '_orient.em'
        resV.write(resFilename)
        orientV.write(orientFilename)
        
        from pytom.localization.structures import Volume, Orientation
        res = Volume(resFilename)
        orient = Orientation(orientFilename)
        
        # construct the result
        from pytom.localization.peak_job import PeakResult
        result = PeakResult(res, orient, jobID)
        
        return result
    
    def summarize(self, result, jobID, verbose=True):
        """
        summarize: Get the result and do the summarization accordingly.\
        (Send back, or update result, or write to the disk).

        @param result: result
        @type result: L{pytom.localization.peak_job.PeakResult}
        @param jobID: ID of job
        @type jobID: integer
        """
        resV = result[0]; orientV = result[1]
        
        if self.jobInfoPool == {}: # leaf node
            assert self.backTo != None
            
            result = self.writeRes(resV, orientV, jobID)
            
            result.send(self.mpi_id, self.backTo)
            return
        
        # non leaf node, update the result
        if self.jobInfoPool[jobID].splitType == "Ang":
            self.jobInfoPool["numDoneJobsA"] = self.jobInfoPool["numDoneJobsA"] + 1
            offset = self.jobInfoPool[jobID].angleOffset
#            print self.name + ': JobID ' + str(jobID) + ' Offset ' + str(offset)
            orientV = orientV + offset
            
            if self.resVolA==None or self.resOrientA==None:
                from pytom_volume import vol
                self.resVolA = vol(resV.sizeX(), resV.sizeY(), resV.sizeZ())
                self.resOrientA = vol(orientV.sizeX(), orientV.sizeY(), orientV.sizeZ())
                self.resVolA.copyVolume(resV)
                self.resOrientA.copyVolume(orientV)
            else:
                from pytom_volume import updateResFromVol
                updateResFromVol(self.resVolA, resV, self.resOrientA, orientV)
                
            if self.jobInfoPool["numDoneJobsA"] == self.jobInfoPool["numJobsA"] and self.jobInfoPool["numJobsV"]>0:
                self.summarize([self.resVolA, self.resOrientA], self.jobInfoPool[jobID].originalJobID, verbose)
                return
        elif self.jobInfoPool[jobID].splitType == "Vol":
            self.jobInfoPool["numDoneJobsV"] = self.jobInfoPool["numDoneJobsV"] + 1
            [originX, originY, originZ] = self.jobInfoPool[jobID].origin
            if self.resVol==None or self.resOrient==None:
                from pytom_volume import vol
                [vsizeX, vsizeY, vsizeZ] = self.jobInfoPool[jobID].originalSize

                self.resVol = vol(vsizeX, vsizeY, vsizeZ)
                self.resVol.setAll(0)
                self.resOrient = vol(vsizeX, vsizeY, vsizeZ)
                self.resOrient.setAll(0)
            
            [vsizeX, vsizeY, vsizeZ] = self.jobInfoPool[jobID].originalSize
            [sizeX ,sizeY, sizeZ] = self.jobInfoPool[jobID].splitSize
            sub_start = self.jobInfoPool[jobID].sub_start
            start = self.jobInfoPool[jobID].whole_start

            stepSizeX = min(vsizeX-sub_start[0], sizeX)
            stepSizeY = min(vsizeY-sub_start[1], sizeY)
            stepSizeZ = min(vsizeZ-sub_start[2], sizeZ)

            from pytom_volume import subvolume, putSubVolume
            sub_resV = subvolume(resV, sub_start[0],sub_start[1],sub_start[2], stepSizeX,stepSizeY,stepSizeZ)
            sub_resO = subvolume(orientV, sub_start[0],sub_start[1],sub_start[2], stepSizeX,stepSizeY,stepSizeZ)

            putSubVolume(sub_resV, self.resVol, start[0]-originX,start[1]-originY,start[2]-originZ)
            putSubVolume(sub_resO, self.resOrient, start[0]-originX,start[1]-originY,start[2]-originZ)
        else:
            raise RuntimeError("Unclear split type!")
        
        # if all results are there, write back to disk and return to high level
        if self.jobInfoPool["numDoneJobsA"]+self.jobInfoPool["numDoneJobsV"] == self.jobInfoPool["numJobsA"]+self.jobInfoPool["numJobsV"]:
            if self.jobInfoPool["numJobsV"]==0:
                self.resVol = self.resVolA
                self.resOrient = self.resOrientA
                
            if self.backTo != None:
                result = self.writeRes(self.resVol, self.resOrient, self.jobInfoPool[jobID].originalJobID)
                if verbose==True:
                    print(self.name + ': sending back result to ' + str(self.backTo))
                result.send(self.mpi_id, self.backTo)
            else:
                # write the final result to the disk
                self.writeRes(self.resVol, self.resOrient)
                self.parallelEnd(verbose)
                
    
    def parallelRun(self, job, splitX=0, splitY=0, splitZ=0, verbose=True, gpuID=-1):
        """
        parallelRun: Parallel run the job on the computer cluster.
        @param job: job
        @type job: L{pytom.localization.peak_job.PeakJob}
        @param splitX: split part along the x dimension
        @type splitX: integer
        @param splitY: split part along the y dimension
        @type splitY: integer
        @param splitZ: split part along the z dimension
        @type splitZ: integer
        """
        import pytom_mpi
        if self.mpi_id == 0: # send the first message
#            if not pytom_mpi.isInitialised():
#                pytom_mpi.init()
            job.members = pytom_mpi.size()
            print('job members', job.members)
            job.send(0, 0)
            print("\n")
        
        end = False
        while not end:
            # get the message string
            mpi_msgString = getMsgStr()
            
            msgType = self.getMsgType(mpi_msgString)
            
            if msgType == 2: # Job msg
                msg = self.getJobMsg(mpi_msgString)
                job = self.jobFromMsg(msg) # set members
                
                if self.mpi_id == 0:
                    self.distributeJobs(job, splitX, splitY, splitZ)
                else:
                    self.distributeJobs(job)
                
                result = self.run(verbose, gpuID=gpuID)
                self.summarize(result, self.jobID)
                
            elif msgType == 1: # Result msg
                msg = self.getResMsg(mpi_msgString)
                res = self.resFromMsg(msg)
                
                if verbose == True:
                    print(self.name + ": processing result from worker " + msg.getSender())
                    
                resV = res.result.getVolume()
                resO = res.orient.getVolume()
                jobID = res.jobID

                self.summarize([resV, resO], jobID)
                
            elif msgType == 0: # Status msg
                # get the message as StatusMessage and finish
                from pytom.parallel.messages import StatusMessage
                msg = StatusMessage('','')
                msg.fromStr(mpi_msgString)
                if msg.getStatus() == 'End':
                    end = True
                    if verbose==True:
                        print(self.name + ': end')
            else: # Error
                raise RuntimeError("False message type!")
        
        if self.mpi_id == 0:
            # delete the temporary files on the disk
            import os
            files = os.listdir(self.dstDir)
            for name in files:
                if 'job' in name and '.em' in name and not 'sum' in name and not 'sqr' in name:
                    os.remove(self.dstDir+'/'+name)
            
            # rename the result files name
            os.rename(self.dstDir+'/'+'node_0_res.em', self.dstDir+'/'+'scores_{}.em'.format(self.suffix))
            os.rename(self.dstDir+'/'+'node_0_orient.em', self.dstDir+'/'+'angles_{}.em'.format(self.suffix))
        
        self.clean() # clean itself
        pytom_mpi.finalise()

        
    def parallelEnd(self, verbose=True):
        """
        parallelEnd: End the parallel running of the program.
        """
        import pytom_mpi
        from pytom.parallel.messages import StatusMessage
        
        if verbose == True:
            print(self.name + ': sending end messages to all')
        
        for i in range(pytom_mpi.size()):
            msg = StatusMessage(str(self.mpi_id), str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg), i)
            
            
    def parallelRunMultiJobs(self, jobs, splitX=0, splitY=0, splitZ=0, verbose=True):
        for i in range(len(jobs)):
            job = jobs[i]
            
            if self.mpi_id == 0:
                from pytom.tools.files import checkDirExists
                new_dir = job.dstDir + 'Job_' + str(i)
                if not checkDirExists(new_dir):
                    import os
                    os.mkdir(new_dir)
                job.dstDir = new_dir
            
            self.parallelRun(job, splitX, splitY, splitZ, verbose)
