'''
FRM alignment classes
Created on Mar 5, 2012

@author: yuxiangchen
'''

from pytom.basic.structures import PyTomClass
import pytom.lib.pytom_mpi as pytom_mpi
import os

class FRMJob(PyTomClass): # i need to rename the class, but for now it works
    def __init__(self, pl=None, ref=None, mask=None, peak_offset=0, sample_info=None, 
            bw_range=None, freq=None, dest='.', max_iter=10, r_score=False, 
            weighting=False, bfactor=None, symmetries=None, adaptive_res=0.1, 
            fsc_criterion=0.5, constraint=None, binning=1):
        """
        initiate FRM job
        @param pl: particle list
        @type ps: L{pytom.basic.structures.ParticleList}
        @param ref: reference density
        @type ref: L{pytom.basic.structures.Reference}
        @param mask: mask 
        @type ref: L{pytom.basic.structures.Mask}
        @param peak_offset: peak offset in voxel
        @type peak_offset: C{int}
        @param sample_info: ?? (Default: None)
        @type sample_info: ??
        @param bw_range: bandwidth range in pixel (2-dim vector)
        @type bw_range: C{list}
        @param freq: frequency (default: None)
        @type: C{int}
        @param dest: distination directory (default: '.')
        @type: C{str}
        @param max_iter: maximum number of iterations
        @type max_iter: C{int}
        @param r_score: use r_score (??) (default: False)
        @type r_score: C{bool}
        @param weighting: weighting (default: False)
        @type weighting: C{bool}
        @param bfactor: B-factor (default: None)
        @type bfactor: C{float}?
        @param symmetries: symmetry (default: None)
        @type L{pytom.basic.structures.Symmetries}
        @param adaptive_res: adaptive resolution - add to resolution for filtering
        @type adaptive_res: C{float}
        @param fsc_criterion: FSC criterion (default: 0.5)
        @type fsc_criterion: C{float}
        @param constraint: Constraint on orientations (deafult: None)
        @type constraint: ??
        @param binning: Perform binning (downscale) of subvolumes by factor. Default=1.
        @type binning C{float}
        """
        self.particleList = pl
        self.reference = ref
        self.mask = mask
        self.peak_offset = peak_offset
        self.sampleInformation = sample_info
        self.bw_range = bw_range
        self.freq = freq
        self.destination = dest
        self.max_iter = max_iter
        self.r_score = r_score
        self.weighting = weighting
        self.bfactor = bfactor
        self.symmetries = symmetries
        self.adaptive_res = adaptive_res
        self.fsc_criterion = fsc_criterion
        self.constraint = constraint
        self.binning = binning
    
    def fromXML(self, xmlObj):
        """
        read from xml file
        @param xmlObj: xml object
        @type xmlObj: L{lxml.etree.Element}
        """
        from lxml.etree import _Element


        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "FRMJob":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('FRMJob')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a FRMJob.")
            jobDescription = jobDescription[0]
        
        from pytom.basic.structures import ParticleList, Reference, Mask, SampleInformation, MultiSymmetries
        
        pl = ParticleList('.')
        particleList_element = jobDescription.xpath('ParticleList')
        if len(particleList_element) > 0:
            pl.fromXML(particleList_element[0])
        else:
            list_elements = jobDescription.xpath('ParticleListLocation')
            for e in list_elements:
                sub_pl = ParticleList()
                sub_pl.fromXMLFile(e.get('Path'))
                pl += sub_pl
        self.particleList = pl
        
        r = jobDescription.xpath('Reference')[0]
        self.reference = Reference('')
        self.reference.fromXML(r)
        
        m = jobDescription.xpath('Mask')[0]
        self.mask = Mask('')
        self.mask.fromXML(m)
        
        try:
            si = jobDescription.xpath('SampleInformation')[0]
            self.sampleInformation = SampleInformation()
            self.sampleInformation.fromXML(si)
        except:
            self.sampleInformation = SampleInformation()
        
        try:
            syms = jobDescription.xpath('MultiSymmetries')[0]
            self.symmetries = MultiSymmetries()
            self.symmetries.fromXML(syms)
        except:
            self.symmetries = MultiSymmetries()

        self.peak_offset = int(jobDescription.get('PeakOffset'))
        self.bw_range = [int(i) for i in jobDescription.get('BandwidthRange')[1:-1].split(',')]
        self.freq = int(jobDescription.get('Frequency'))
        self.destination = jobDescription.get('Destination')
        self.max_iter = int(jobDescription.get('MaxIterations'))
        self.r_score = jobDescription.get('RScore')=='True'
        self.weighting = jobDescription.get('WeightedAverage')=='True'
        self.bfactor = jobDescription.get('BFactor')
        self.binning = int(jobDescription.get('Binning'))
        if jobDescription.get('AdaptiveResolution'):
            adaptive_resolution = jobDescription.get('AdaptiveResolution')
            if adaptive_resolution == '+1':
                self.adaptive_res = False # always increase by 1
            else:
                self.adaptive_res = float(adaptive_resolution)
        else:
            self.adaptive_res = 0.0 # default, if not specified
        if jobDescription.get('FSC'):
            self.fsc_criterion = float(jobDescription.get('FSC'))
        else:
            self.fsc_criterion = 0.5 # default value
        
        # for the constraint
        try:
            from sh_alignment.constrained_frm import AngularConstraint
            con = jobDescription.xpath('AngularConstraint')
            if len(con) != 0:
                ac = AngularConstraint()
                c = ac.fromXML(con[0])
                self.constraint = c
            else:
                self.constraint = None
        except:
            self.constraint = None

    
    def toXML(self):
        """
        copy to xml structure
        @return: xml object for job
        @rtype L{lxml.etree.Element}
        """
        from lxml import etree

        jobElement = etree.Element("FRMJob")
        jobElement.append(self.particleList.toXML())
        jobElement.append(self.reference.toXML())
        jobElement.append(self.mask.toXML())
        jobElement.append(self.sampleInformation.toXML())
        if self.symmetries is not None:
            jobElement.append(self.symmetries.toXML())
        jobElement.set("PeakOffset", str(self.peak_offset))
        jobElement.set("BandwidthRange", str(self.bw_range))
        jobElement.set("Frequency", str(self.freq))
        jobElement.set("Destination", self.destination)
        jobElement.set("MaxIterations", str(self.max_iter))
        jobElement.set("RScore", str(self.r_score))
        jobElement.set("WeightedAverage", str(self.weighting))
        jobElement.set("BFactor", str(self.bfactor))
        jobElement.set("Binning", str(self.binning))
        if self.adaptive_res is False:
            jobElement.set("AdaptiveResolution", '+1')
        else:
            jobElement.set("AdaptiveResolution", str(self.adaptive_res))
        jobElement.set("FSC", str(self.fsc_criterion))
        if self.constraint:
            jobElement.append(self.constraint.toXML())
        
        return jobElement
    
    def check(self):
        from pytom.tools.files import checkDirExists
        self.particleList.check()
        self.reference.check()
        self.mask.check()
        if not checkDirExists(self.destination):
            raise RuntimeError('Destination path not found! ' + self.destination)

from pytom.basic.score import Score
class FRMScore(Score):
    def __init__(self,value=0):
        self.ctor()
        self._type = 'FRMScore'
        self.setValue(value)
    
    def getWorstValue(self):
        return -100000000000.0
    
class FRMResult(PyTomClass):
    def __init__(self, name='', pl=None, worker_id=0):
        self.name = name
        self.pl = pl
        self.worker_id = worker_id
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "FRMResult":
            result = xmlObj
        else:  
            result = xmlObj.xpath('FRMResult')
            if len(result) == 0:
                raise Exception("This XML is not a FRMResult.")
            result = result[0]
        
        from pytom.basic.structures import ParticleList
        
        particleList_element = result.xpath('ParticleList')[0]
        pl = ParticleList('.')
        pl.fromXML(particleList_element)
        self.pl = pl
        
        self.name = result.get('Name')
        self.worker_id = int(result.get('WorkerID'))
    
    def toXML(self):
        from lxml import etree

        jobElement = etree.Element("FRMResult")
        jobElement.append(self.pl.toXML())
        jobElement.set("Name", str(self.name))
        jobElement.set("WorkerID", str(self.worker_id))
        
        return jobElement

class FRMWorker():
    def __init__(self):
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()
            
        self.mpi_id = pytom_mpi.rank()
        self.num_workers = pytom_mpi.size()-1
        self.node_name = 'node_' + str(self.mpi_id)
        
        if self.num_workers < 1:
            raise RuntimeError("Not enough nodes to parallelize the job!")
        
    def start(self, job, verbose=False):
        """
        start FRM job
        @param job: FRM job
        @type job: L{FRMJob}
        @param verbose: print stuff (default: False)
        @type verbose: C{bool}
        """
        if self.mpi_id == 0:
            from pytom.basic.structures import ParticleList, Reference
            from pytom.basic.resolution import bandToAngstrom
            from pytom.basic.filter import lowpassFilter
            from math import ceil
            self.destination = job.destination
            new_reference = job.reference
            old_freq = job.freq
            new_freq = job.freq
            #print(f"reference  = {job.reference}")
            #print(f"particlelist = {job.particleList}")
            print(f"iterations = {job.max_iter:d}")
            print(f"binning    = {job.binning:d}")
            #print(f"mask       = {job.mask}")
            #print(f"peak_offset= {job.peak_offset:f2.1}")
            print(f"destination= {job.destination:s}")
            print(f"freq cut   = {job.freq:d}")
            # main node
            for i in range(job.max_iter):
                if verbose:
                    print(self.node_name + ': starting iteration %d ...' % i)
                
                # construct a new job by updating the reference and the frequency
                new_job = FRMJob(job.particleList, new_reference, job.mask, 
                                 job.peak_offset, job.sampleInformation, job.bw_range, 
                                 new_freq, job.destination, job.max_iter-i, job.r_score, 
                                 job.weighting, constraint=job.constraint, 
                                 binning=job.binning)
                
                # distribute it
                self.distribute_job(new_job, verbose)
                
                # get the result back
                all_even_pre = None
                all_even_wedge = None
                all_odd_pre = None
                all_odd_wedge = None
                pl = ParticleList()
                for j in range(self.num_workers):
                    result = self.get_result()
                    pl += result.pl
                    even_pre, even_wedge, odd_pre, odd_wedge = self.retrieve_res_vols(result.name)
                    
                    if all_even_pre:
                        all_even_pre += even_pre
                        all_even_wedge += even_wedge
                        all_odd_pre += odd_pre
                        all_odd_wedge += odd_wedge
                    else:
                        all_even_pre = even_pre
                        all_even_wedge = even_wedge
                        all_odd_pre = odd_pre
                        all_odd_wedge = odd_wedge
                
                # write the new particle list to the disk
                pl.toXMLFile(os.path.join(job.destination, 'aligned_pl_iter'+str(i)+'.xml'))
                
                # create half sets
                even = self.create_average(all_even_pre, all_even_wedge)
                odd = self.create_average(all_odd_pre, all_odd_wedge)
                
                # apply symmetries before determine resolution
                even = job.symmetries.applyToParticle(even)
                odd = job.symmetries.applyToParticle(odd)
                resNyquist, resolutionBand, numberBands = self.determine_resolution(even, odd, job.fsc_criterion, None, job.mask, verbose)
                
                # write the half set to the disk
                even.write(os.path.join(self.destination, 'fsc_'+str(i)+'_even.em'))
                odd.write(os.path.join(self.destination, 'fsc_'+str(i)+'_odd.em'))
                all_even_wedge.write(os.path.join(self.destination, 'fsc_' + str(i) + '_even_wedge.em'))
                all_odd_wedge.write(os.path.join(self.destination, 'fsc_' + str(i) + '_odd_wedge.em'))
                
                # determine the resolution
                if verbose:
                    print(self.node_name + ': determining the resolution ...')
                current_resolution = bandToAngstrom(resolutionBand, job.sampleInformation.getPixelSize(), numberBands, 1)
                if verbose:
                    print(self.node_name + ': current resolution ' + str(current_resolution), resNyquist)
                
                # create new average
                all_even_pre += all_odd_pre
                all_even_wedge += all_odd_wedge
                average = self.create_average(all_even_pre, all_even_wedge)
                
                # apply symmetries
                average = job.symmetries.applyToParticle(average)
                
                # filter average to resolution and update the new reference
                average_name = os.path.join(self.destination, 'average_iter'+str(i)+'.em')
#                pl.average(average_name, True)
                average.write(average_name)
                new_reference = Reference(average_name)
                
                # low pass filter the reference and write it to the disk
                filtered = lowpassFilter(average, ceil(resolutionBand), ceil(resolutionBand)/10)
                filtered_ref_name = os.path.join(self.destination, 'average_iter'+str(i)+'_res'+
                                                 f'{current_resolution:.2f}'+'.em')
                filtered[0].write(filtered_ref_name)
                
                # if the position/orientation is not improved, break it
                
                # change the frequency to a higher value
                new_freq = int(ceil(resolutionBand))+1
                if new_freq <= old_freq:
                    if job.adaptive_res is not False: # two different strategies
                        print(self.node_name + ': Determined resolution gets worse. Include additional %f percent frequency to be aligned!' % job.adaptive_res)
                        new_freq = int((1+job.adaptive_res)*new_freq)
                        old_freq = new_freq
                    else: # always increase by 1
                        print(self.node_name + ': Determined resolution gets worse. Increase the frequency to be aligned by 1!')
                        new_freq = old_freq+1
                        old_freq = new_freq
                else:
                    old_freq = new_freq
                if new_freq >= numberBands:
                    print(self.node_name + ': New frequency too high. Terminate!')
                    break
                
                if verbose:
                    print(self.node_name + ': change the frequency to ' + str(new_freq))
            
            # send end signal to other nodes and terminate itself
            self.end(verbose)
        else:
            # other nodes
            self.run(verbose)
    
    def end(self, verbose=False):
        if verbose == True:
            print(self.node_name + ': sending end messages to others')
        
        from pytom.parallel.messages import StatusMessage
        
        mpi_numberNodes = pytom_mpi.size()
        mpi_myid = pytom_mpi.rank()
        
        for i in range(1, mpi_numberNodes):
            msg = StatusMessage(str(mpi_myid),str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg),i)
        
        pytom_mpi.finalise()
    
    def run(self, verbose=False):
        from sh_alignment.frm import frm_align
        from sh_alignment.constrained_frm import frm_constrained_align, AngularConstraint
        from pytom.basic.structures import Shift, Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.transformations import resize, resizeFourier
        binningType = 'Fourier'
        
        while True:
            # get the job
            try:
                job = self.get_job()
                print('get job particle')
            except Exception as e:
                print("Error: ", e)
                if verbose:
                    print(self.node_name + ': end')
                break # get some non-job message, break it
            
            if verbose:
                prog = FixedProgBar(0, len(job.particleList)-1, self.node_name+':')
                i = 0
            ref = job.reference.getVolume()
            if job.binning > 1:
                ref = resize(volume=ref, factor=1./job.binning, interpolation=binningType)
                if type(ref) == tuple:
                    ref = ref[0]
            # re-set max frequency in case it exceeds Nyquist - a bit brute force
            job.freq = min(job.freq, ref.sizeX()//2-1)
            # run the job
            for p in job.particleList:
                if verbose:
                    prog.update(i)
                    i += 1
                v = p.getVolume()
                if job.binning > 1:
                    v = resize(volume=v, factor=1./job.binning, interpolation=binningType)
                    if type(v) == tuple:
                        v = v[0]
                mask = job.mask.getVolume()
                if job.binning > 1:
                    mask = resize(volume=mask, factor=1./job.binning, interpolation='Spline')
                    if type(mask) == tuple:
                        mask = mask[0]
                if job.constraint:
                    constraint = job.constraint
                    if job.constraint.type == AngularConstraint.ADP_ANGLE: # set the constraint around certain angle
                        rot = p.getRotation()
                        constraint.setAngle(rot.getPhi(), rot.getPsi(), rot.getTheta())
                    #pos, angle, score = frm_constrained_align(v, p.getWedge(), ref, None, job.bw_range, job.freq, job.peak_offset, job.mask.getVolume(), constraint)
                    if job.binning >1:
                        pos, angle, score = frm_constrained_align(v, p.getWedge(), ref, None, 
                                                job.bw_range, job.freq, job.peak_offset/job.binning,
                                                mask, constraint)
                    else:
                        pos, angle, score = frm_constrained_align(v, p.getWedge(), ref, None, 
                                                job.bw_range, job.freq, job.peak_offset, mask, constraint)
                else:
                    #pos, angle, score = frm_align(v, p.getWedge(), ref, None, job.bw_range, job.freq, job.peak_offset, job.mask.getVolume())
                    #if job.binning >1:
                    #    print(job.peak_offset)
                    #    print(type(job.peak_offset))
                    #    print(job.peak_offset/job.binning)
                    #    print(type(job.binning))
                    #    pos, angle, score = frm_align(v, p.getWedge(), ref, None, job.bw_range, job.freq, 
                    #                            job.peak_offset/job.binning, mask)
                    #else:
                    pos, angle, score = frm_align(v, p.getWedge(), ref, None, job.bw_range, job.freq, 
                                                job.peak_offset, mask)

                if job.binning > 1:
                    pos[0] = job.binning*(pos[0]-v.sizeX()/2)
                    pos[1] = job.binning*(pos[1]-v.sizeY()/2) 
                    pos[2] = job.binning*(pos[2]-v.sizeZ()/2)
                    p.setShift(Shift([pos[0], pos[1], pos[2]]))
                else:     
                    p.setShift(Shift([pos[0]-v.sizeX()/2, pos[1]-v.sizeY()/2, pos[2]-v.sizeZ()/2]))
                p.setRotation(Rotation(angle))
                p.setScore(FRMScore(score))
                
            # average the particle list
            name_prefix = os.path.join(job.destination, self.node_name+'_'+str(job.max_iter))
            self.average_sub_pl(job.particleList, name_prefix, job.weighting)
            
            # send back the result
            self.send_result(FRMResult(name_prefix, job.particleList, self.mpi_id))
        
        pytom_mpi.finalise()
    
    def average_sub_pl(self, pl, name_prefix, weight_average):
        """For worker node, do two things, averaging & obtaining the even/odd partitions to save some time.
           @param pl: particle list
           @type ps: L{pytom.basic.structures.ParticleList}
           @param name_prefix: name prefix output densities
           @type name_prefix: C{str}
           @param weight_average: weighted average
           @type weight_average: C{str}
        """
        from pytom.basic.structures import ParticleList
        even = ParticleList('.')
        odd = ParticleList('.')
        
        for i in range(len(pl)):
            if i%2 == 0:
                even.append(pl[i])
            else:
                odd.append(pl[i])
        
        even.average(name_prefix+'even.em', progressBar=False, createInfoVolumes=False, _mpiParallel=False, weighting=weight_average)
        odd.average(name_prefix+'odd.em', progressBar=False, createInfoVolumes=False, _mpiParallel=False, weighting=weight_average)
    
    def retrieve_res_vols(self, name_prefix):
        """For master node, retrieve the even/odd sub-averages and do the cleaning.
           @param name_prefix: name prefix output densities
           @type name_prefix: C{str}
        """
        from pytom.lib.pytom_volume import read
        even_pre = read(name_prefix+'even'+'-PreWedge.em')
        even_wedge = read(name_prefix+'even'+'-WedgeSumUnscaled.em')
        odd_pre = read(name_prefix+'odd'+'-PreWedge.em')
        odd_wedge = read(name_prefix+'odd'+'-WedgeSumUnscaled.em')
        
        # delete the volumes from disk
        import os
        try:
            os.remove(name_prefix+'even'+'-PreWedge.em')
            os.remove(name_prefix+'even'+'-WedgeSumUnscaled.em')
            
            os.remove(name_prefix+'odd'+'-PreWedge.em')
            os.remove(name_prefix+'odd'+'-WedgeSumUnscaled.em')
            
            os.remove(name_prefix+'even'+'.em')
            os.remove(name_prefix+'odd'+'.em')
        except:
            pass
        
        return (even_pre, even_wedge, odd_pre, odd_wedge)
    
    def create_average(self, pre, wedge):
        """For the master node, create the average according to the pre-wedge and wedge volumes.
           @param pre: density prior to weighting
           @type pre: L{pytom.lib.pytom_volume.vol}
           @param wedge: wedge 
           @type wedge: L{pytom.basic.Wedge}
           @return: wedge-weighted density
           @rtype: L{pytom.lib.pytom_volume.vol}
        """
        from pytom.lib.pytom_volume import complexDiv, limit
        from pytom.basic.fourier import fft,ifft
        
        limit(wedge, 0.1, 0, 0,0,True,False) # set all the values below the specified value to 0
        
        f_pre = fft(pre)
        r = complexDiv(f_pre, wedge)
        average = ifft(r)
        average.shiftscale(0.0,1/float(average.sizeX()*average.sizeY()*average.sizeZ()))
        
        return average
    
    def determine_resolution(self, even, odd, criterion, numberBands, mask, verbose=False):
        """For the master node, determine the resolution.
           @param even: particle list even
           @type even: L{pytom.basic.structures.ParticleList}
           @param odd: particle list odd
           @type odd: L{pytom.basic.structures.ParticleList}
        """
        from pytom.basic.correlation import FSC, determineResolution
        
        if not numberBands:
            numberBands = even.sizeX()/2
        
        fsc = FSC(even, odd, numberBands, mask, verbose=False)
        if verbose:
            print(self.node_name + ': FSC: ' + str(fsc))
        
        return determineResolution(fsc, criterion, verbose=False)
    
    def send_job(self, job, dest):
        """
        @param job: FRM job
        @type job: L{FRMJob}
        """
        pytom_mpi.send(str(job), dest)
    
    def get_job(self):
        """
        @return: FRM job
        @rtype: L{FRMJob}
        """
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        job = FRMJob()
        job.fromStr(mpi_msgString)
        
        return job
    
    def send_result(self, result):
        pytom_mpi.send(str(result), 0)
    
    def get_result(self):
        """
        @return: FRM result
        @rtype: L{FRMResult}
        """
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        result = FRMResult()
        result.fromStr(mpi_msgString)
        
        return result
    
    def distribute_job(self, job, verbose=False):
        """
        @param job: FRM job
        @type job: L{FRMJob}
        """
        n = len(job.particleList)
        particlesPerNode = int(n/self.num_workers)
        residual = n-particlesPerNode*self.num_workers
        start_idx = 0
        for i in range(1, self.num_workers+1):
            if i < residual+1: # since i starts from 1
                l = particlesPerNode+1
            else:
                l = particlesPerNode
            
            if l < 2:
                raise RuntimeError("Particles per node is less than 2. Please use less nodes!")
            
            subPL = job.particleList[start_idx : start_idx+l]
            start_idx += l
            
            subJob = FRMJob(subPL, job.reference, job.mask, job.peak_offset, 
                            job.sampleInformation, job.bw_range, job.freq, 
                            job.destination, job.max_iter, job.r_score, 
                            job.weighting, constraint=job.constraint, 
                            binning=job.binning)
            self.send_job(subJob, i)
            
            if verbose:
                print(self.node_name + ': distributed %d particles to node %d' % (len(subPL), i))

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Subtomogram alignment by Fast Rotational Matching.',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-j'], 'Job xml file.', True, False),
                                    ScriptOption(['-v'], 'Verbose mode.', False, True),
                                    ScriptOption(['--help'], 'Help info.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    
    try:
        job_filename, verbose, bHelp = parse_script_options(sys.argv[1:], helper)     
    except Exception as e:
        print(e)
        sys.exit()
        
    if bHelp is True:
        print(helper)
        sys.exit()

    print(job_filename)
    # check the job
    job = FRMJob()
    job.fromXMLFile(job_filename)
    job.check()
    
    worker = FRMWorker()
    
    if verbose:
        from pytom.tools.timing import Timing
        t = Timing()
        t.start()
    
    # start it
    worker.start(job, verbose)
    
    if verbose:
        print('Overall execution time: %f s.' % t.end())
        
        
