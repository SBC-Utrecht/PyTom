'''
Created on May 20, 2010

@author: chen
'''
import numpy as np
import os
from pytom.lib.pytom_volume import vol, updateResFromVol
from pytom.lib.pytom_numpy import vol2npy, npy2vol
from pytom.agnostic.tools import subvolume, putSubVolume


def getMsgStr():
    '''
    getMsgStr: Use MPI to receive the message
    @return: message string
    @rtype: str
    '''
    import pytom.lib.pytom_mpi as pytom_mpi
    mpi_msgString = pytom_mpi.receive()
    return mpi_msgString


class PeakWorker(object):
    '''
    PeakWorker: Worker to calculate the score function
    @author: Chen
    '''
    def __init__(self):
        import pytom.lib.pytom_mpi as pytom_mpi
        
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
        
        g = None
        if gpuID is None:  # gpu dependent import of functions
            from pytom.localization.extractPeaks import extractPeaks
        else:
            from pytom.localization.extractPeaks import templateMatchingGPU as extractPeaks
            v, ref, m = [vol2npy(x).copy() for x in (v, ref, m)]
            g = gpuID[self.mpi_id]

        if verbose==True:
            print(self.name + ': starting to calculate %d rotations' % rot.numberRotations() )
        [resV, orientV, sumV, sqrV] = extractPeaks(v, ref, rot, scoreFnc, m, mIsSphere, wedg, nodeName=self.name,
                                                   verboseMode=verbose, moreInfo=moreInfo, gpuID=g,
                                                   jobid=self.mpi_id)
        self.runtimes = self.runtimes + 1

        if g is not None:  # convert the results back to pytomvol for merging between procs
            # vol2npy returns a pytom volume view of a numpy array
            # need to ensure the returned variable does not refer to a locally existing variable
            copy1, copy2 = np.asfortranarray(resV), np.asfortranarray(orientV)
            tmp_res, tmp_orient = npy2vol(copy1, 3), npy2vol(copy2, 3)
            resV, orientV = vol(*resV.shape), vol(*orientV.shape)
            [resV.copyVolume(tmp_res), orientV.copyVolume(tmp_orient)]
        return [resV, orientV, sumV, sqrV]


class PeakLeader(PeakWorker):
    """
    PeakLeader: Class for parallel running of jobs (new architecture)
    """
    def __init__(self,suffix=''):
        import pytom.lib.pytom_mpi as pytom_mpi

        if not pytom_mpi.isInitialised():
            pytom_mpi.init()

        self.suffix=suffix
        self.mpi_id = pytom_mpi.rank()
        self.name = 'node_' + str(self.mpi_id)

        self.size = pytom_mpi.size()

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
        
        # initialize the jobInfo structure
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
            vsize_x = v.size_x(); vsize_y = v.size_y(); vsize_z = v.size_z()
        else:
            origin = job.volume.subregion[0:3]
            vsize_x = job.volume.subregion[3]
            vsize_y = job.volume.subregion[4]
            vsize_z = job.volume.subregion[5]
            
        size_x = vsize_x//splitX; size_y = vsize_y//splitY; size_z = vsize_z//splitZ
        r = job.reference.getVolume()
        rsize_x = r.size_x(); rsize_y = r.size_y(); rsize_z = r.size_z()
        if rsize_x>size_x or rsize_y>size_y or rsize_z>size_z:
            raise RuntimeError("Not big enough volume to split!")
        
        # initialize the jobInfo structure
        originalJobID = job.jobID
            
        # read the target volume, calculate the respective subregion
        from pytom.localization.peak_job import PeakJob
        from pytom.localization.structures import Volume
        _start = [-rsize_x//2+origin[0],-rsize_y//2+origin[1],-rsize_z//2+origin[2]]
        _size = [size_x+rsize_x, size_y+rsize_y, size_z+rsize_z]
        
        numPieces = splitX*splitY*splitZ
        totalMem = self.members
        numMemEach = totalMem//numPieces
        targetID = self.mpi_id
        
        for i in range(numPieces):
            strideZ = splitX*splitY; strideY = splitX
            incZ = i//strideZ; incY = (i%strideZ)//strideY; incX = i%strideY
            _start = [-rsize_x//2+origin[0]+incX*size_x,-rsize_y//2+origin[1]+incY*size_y,-rsize_z//2+origin[2]+incZ*size_z]
            
            start = _start[:]
            end = [start[j]+_size[j] for j in range(len(start))]
            
            if start[0] < origin[0]:
                start[0]=origin[0]
            if start[1] < origin[1]:
                start[1]=origin[1]
            if start[2] < origin[2]:
                start[2]=origin[2]
            if end[0] > vsize_x+origin[0]:
                end[0] = vsize_x+origin[0]
            if end[1] > vsize_y+origin[1]:
                end[1] = vsize_y+origin[1]
            if end[2] > vsize_z+origin[2]:
                end[2] = vsize_z+origin[2]
            
            size = [end[j]-start[j] for j in range(len(start))]
            
            # make sure that the last dimension is not odd
            # if size[2]%2 == 1:
            #     size[2] = size[2] - 1
            #     end[2] = end[2] - 1
            
            # for reassembling the result
            whole_start = start[:]
            sub_start = [0, 0, 0]
            if start[0] != origin[0]:
                whole_start[0] = start[0]+rsize_x//2
                sub_start[0] = rsize_x//2
            if start[1] != origin[1]:
                whole_start[1] = start[1]+rsize_y//2
                sub_start[1] = rsize_y//2
            if start[2] != origin[2]:
                whole_start[2] = start[2]+rsize_z//2
                sub_start[2] = rsize_z//2
            
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
            info.originalSize = [vsize_x, vsize_y, vsize_z]
            info.splitSize = [size_x, size_y, size_z]
            info.origin = origin
            self.jobInfoPool[subJobID] = info
            
            if targetID == self.mpi_id:
                self.setJob(subJob)
                if self.members > 1:
                    self.splitAngles(subJob, verbose)
            else:
                if verbose==True:
                    print(self.name + ' : send part of the volume to ' + str(targetID))
                subJob.send(self.mpi_id, targetID)
            
            targetID = targetID + numMem
            self.jobInfoPool["numJobsV"] = self.jobInfoPool["numJobsV"] + 1
            
            
    def distributeJobs(self, job, splitX=0, splitY=0, splitZ=0):
        """
        distributeJobs: Distribute the jobs to the members (if any), as well as itself. \
        This function should be called only once for each node.\
        Depending on the members number, the node will decide how to split the job automatically.

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
                self.splitVolumes(job, splitX, splitY, splitZ)
        
    
    def writeRes(self, resV, orientV, jobID=None):
        """
        writeRes: Write the result back to the disk, and return the PeakJobResult.
        @param resV: result volume
        @type resV: L{pytom.lib.pytom_volume.vol}
        @param orientV: orientation volume
        @type orientV: L{pytom.lib.pytom_volume.vol}
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
            if not self.backTo is None:
                result =  self.writeRes(resV, orientV, jobID)
                result.send(self.mpi_id, self.backTo)
            else:
                self.writeRes(resV,orientV,None)
                self.parallelEnd(verbose)
            return
        
        # non leaf node, update the result
        if self.jobInfoPool[jobID].splitType == "Ang":
            self.jobInfoPool["numDoneJobsA"] = self.jobInfoPool["numDoneJobsA"] + 1
            offset = self.jobInfoPool[jobID].angleOffset
#            print self.name + ': JobID ' + str(jobID) + ' Offset ' + str(offset)
            orientV = orientV + offset
            
            if self.resVolA is None or self.resOrientA is None:
                self.resVolA = vol(resV.size_x(), resV.size_y(), resV.size_z())
                self.resOrientA = vol(orientV.size_x(), orientV.size_y(), orientV.size_z())
                self.resVolA.copyVolume(resV)
                self.resOrientA.copyVolume(orientV)
            else:
                updateResFromVol(self.resVolA, resV, self.resOrientA, orientV)
                
            if self.jobInfoPool["numDoneJobsA"] == self.jobInfoPool["numJobsA"] and self.jobInfoPool["numJobsV"]>0:
                self.summarize([self.resVolA, self.resOrientA], self.jobInfoPool[jobID].originalJobID, verbose)
                return
        elif self.jobInfoPool[jobID].splitType == "Vol":
            self.jobInfoPool["numDoneJobsV"] = self.jobInfoPool["numDoneJobsV"] + 1
            [originX, originY, originZ] = self.jobInfoPool[jobID].origin
            if self.resVol is None or self.resOrient is None:
                [vsize_x, vsize_y, vsize_z] = self.jobInfoPool[jobID].originalSize
                self.resVol = vol(vsize_x, vsize_y, vsize_z)
                self.resVol.setAll(0)
                self.resOrient = vol(vsize_x, vsize_y, vsize_z)
                self.resOrient.setAll(0)
            
            [vsize_x, vsize_y, vsize_z] = self.jobInfoPool[jobID].originalSize
            [size_x ,size_y, size_z] = self.jobInfoPool[jobID].splitSize
            sub_start = self.jobInfoPool[jobID].sub_start
            start = self.jobInfoPool[jobID].whole_start

            stepSizeX = min(vsize_x-sub_start[0], size_x)
            stepSizeY = min(vsize_y-sub_start[1], size_y)
            stepSizeZ = min(vsize_z-sub_start[2], size_z)

            sub_resV = subvolume(resV, sub_start[0],sub_start[1],sub_start[2], stepSizeX,stepSizeY,stepSizeZ)
            sub_resO = subvolume(orientV, sub_start[0],sub_start[1],sub_start[2], stepSizeX,stepSizeY,stepSizeZ)

            # print('aaa:  ', self.resVol.__class__, self.resVol.sum(), sub_resV.shape, sub_resV.sum(),  start[0]-originX, start[1]-originY,start[2]-originZ)

            putSubVolume(sub_resV, self.resVol, start[0]-originX, start[1]-originY,start[2]-originZ)
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
        import pytom.lib.pytom_mpi as pytom_mpi
        if self.mpi_id == 0: # send the first message
            if not pytom_mpi.isInitialised():
                pytom_mpi.init()
            job.members = pytom_mpi.size()
            print('job members', job.members)
            self.distributeJobs(job, splitX, splitY, splitZ)
            result = self.run(verbose, gpuID=gpuID)
            self.summarize(result, self.jobID)
            #job.send(0, 0)

        self.gpuID = gpuID
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
            files = os.listdir(self.dstDir)
            for name in files:
                if 'job' in name and '.em' in name and not 'sum' in name and not 'sqr' in name:
                    os.remove(self.dstDir+'/'+name)

            if self.gpuID:
                gpuflag = ''
            else:
                gpuflag = ''
            # rename the result files name
            os.rename(self.dstDir+'/'+'node_0_res.em', self.dstDir+'/'+'scores_{}{}.em'.format(self.suffix, gpuflag))
            os.rename(self.dstDir+'/'+'node_0_orient.em', self.dstDir+'/'+'angles_{}{}.em'.format(self.suffix, gpuflag))
        
        self.clean() # clean itself
        pytom_mpi.finalise()

        
    def parallelEnd(self, verbose=True):
        """
        parallelEnd: End the parallel running of the program.
        """
        import pytom.lib.pytom_mpi as pytom_mpi
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
                    os.mkdir(new_dir)
                job.dstDir = new_dir
            
            self.parallelRun(job, splitX, splitY, splitZ, verbose)
