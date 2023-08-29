'''
Created on Sep 3, 2012

@author: yuxiangchen
'''

from pytom.basic.structures import PyTomClass
import pytom.lib.pytom_mpi as pytom_mpi
import os

class FRMJob(PyTomClass): # i need to rename the class, but for now it works
    def __init__(self, pl=None, ref=None, mask=None, peak_offset=0, sample_info=None, bw_range=None, freq=None, dest='.', max_iter=10, r_score=False, weighting=False, bfactor=None, symmetries=None, adaptive_res=0.1, fsc_criterion=0.5):
        self.particleList = pl
        if ref is None:
            self.reference = []
        else:
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
    
    def fromXML(self, xmlObj):
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
        
        self.reference = []
        r = jobDescription.xpath('Reference')
        for ref_obj in r:
            ref  = Reference('')
            ref.fromXML(ref_obj)
            self.reference.append(ref)
        
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
    
    def toXML(self):
        from lxml import etree

        jobElement = etree.Element("FRMJob")
        jobElement.append(self.particleList.toXML())
        for ref in self.reference:
            jobElement.append(ref.toXML())
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
        if self.adaptive_res is False:
            jobElement.set("AdaptiveResolution", '+1')
        else:
            jobElement.set("AdaptiveResolution", str(self.adaptive_res))
        jobElement.set("FSC", str(self.fsc_criterion))
        
        return jobElement
    
    def check(self):
        from pytom.tools.files import checkDirExists
        self.particleList.check()
        for ref in self.reference:
            ref.check()
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
        if self.mpi_id == 0:
            from pytom.basic.structures import ParticleList, Reference
            from pytom.basic.resolution import bandToAngstrom
            from pytom.basic.filter import lowpassFilter
            from math import ceil
            
            # randomly split the particle list into 2 half sets
            if len(job.particleList.splitByClass()) != 2:
                import numpy as np
                n = len(job.particleList)
                labels = np.random.randint(2, size=(n,))
                print(self.node_name + ': Number of 1st half set:', n-np.sum(labels), 'Number of 2nd half set:', np.sum(labels))
                for i in range(n):
                    p = job.particleList[i]
                    p.setClass(labels[i])

            self.destination = job.destination
            new_reference = job.reference
            old_freq = job.freq
            new_freq = job.freq
            # main node
            for i in range(job.max_iter):
                if verbose:
                    print(self.node_name + ': starting iteration %d ...' % i)
                
                # construct a new job by updating the reference and the frequency
                new_job = FRMJob(job.particleList, new_reference, job.mask, job.peak_offset, job.sampleInformation, job.bw_range, new_freq, job.destination, job.max_iter-i, job.r_score, job.weighting)
                
                # distribute it
                self.distribute_job(new_job, verbose)
                
                # get the result back
                all_even_pre = None # the 1st set
                all_even_wedge = None
                all_odd_pre = None # the 2nd set
                all_odd_wedge = None
                pl = ParticleList()
                for j in range(self.num_workers):
                    result = self.get_result()
                    pl += result.pl
                    pre, wedge = self.retrieve_res_vols(result.name)
                    
                    if self.assignment[result.worker_id] == 0:
                        if all_even_pre:
                            all_even_pre += pre
                            all_even_wedge += wedge
                        else:
                            all_even_pre = pre
                            all_even_wedge = wedge
                    else:
                        if all_odd_pre:
                            all_odd_pre += pre
                            all_odd_wedge += wedge
                        else:
                            all_odd_pre = pre
                            all_odd_wedge = wedge
                
                # write the new particle list to the disk
                pl.toXMLFile('aligned_pl_iter'+str(i)+'.xml')
                
                # create the averages separately
                if verbose:
                    print(self.node_name + ': determining the resolution ...')
                even = self.create_average(all_even_pre, all_even_wedge)
                odd = self.create_average(all_odd_pre, all_odd_wedge)
                
                # apply symmetries if any
                even = job.symmetries.applyToParticle(even)
                odd = job.symmetries.applyToParticle(odd)
                
                # determine the transformation between even and odd
                # here we assume the wedge from both sets are fully sampled
                from pytom.lib.frm import frm_align
                pos, angle, score = frm_align(odd, None, even, None, job.bw_range, new_freq, job.peak_offset)
                print(self.node_name + 'Transform of even set to match the odd set - shift: '+str(pos)+' rotation: '+str(angle))
                
                # transform the odd set accordingly
                from pytom.lib.pytom_volume import vol, transformSpline
                from pytom.basic.fourier import ftshift
                from pytom.lib.pytom_volume import reducedToFull
                from pytom.lib.pytom_freqweight import weight
                transformed_odd_pre = vol(odd.size_x(), odd.size_y(), odd.size_z())
                full_all_odd_wedge = reducedToFull(all_odd_wedge)
                ftshift(full_all_odd_wedge)
                odd_weight = weight(full_all_odd_wedge) # the funny part of pytom
                transformed_odd = vol(odd.size_x(), odd.size_y(), odd.size_z())
                
                transformSpline(all_odd_pre,transformed_odd_pre,-angle[1],-angle[0],-angle[2],odd.size_x()/2,odd.size_y()/2,odd.size_z()/2,-(pos[0]-odd.size_x()/2),-(pos[1]-odd.size_y()/2),-(pos[2]-odd.size_z()/2),0,0,0)
                odd_weight.rotate(-angle[1],-angle[0],-angle[2])
                transformed_odd_wedge = odd_weight.getWeightVolume(True)
                transformSpline(odd,transformed_odd,-angle[1],-angle[0],-angle[2],odd.size_x()/2,odd.size_y()/2,odd.size_z()/2,-(pos[0]-odd.size_x()/2),-(pos[1]-odd.size_y()/2),-(pos[2]-odd.size_z()/2),0,0,0)
                
                all_odd_pre = transformed_odd_pre
                all_odd_wedge = transformed_odd_wedge
                odd = transformed_odd
                
                # determine resolution
                resNyquist, resolutionBand, number_bands = self.determine_resolution(even, odd, job.fsc_criterion, None, job.mask, verbose)
                
                # write the half set to the disk
                even.write(os.path.join(self.destination, 'fsc_'+str(i)+'_even.em'))
                odd.write(os.path.join(self.destination, 'fsc_'+str(i)+'_odd.em'))
                
                current_resolution = bandToAngstrom(resolutionBand, job.sampleInformation.getPixelSize(), number_bands, 1)
                if verbose:
                    print(self.node_name + ': current resolution ' + str(current_resolution), resNyquist)
                
                # create new average
                all_even_pre += all_odd_pre
                all_even_wedge += all_odd_wedge
                average = self.create_average(all_even_pre, all_even_wedge)
                
                # apply symmetries
                average = job.symmetries.applyToParticle(average)
                
                # filter average to resolution 
                average_name = os.path.join(self.destination, 'average_iter'+str(i)+'.em')
                average.write(average_name)
                
                # update the references
                new_reference = [Reference(os.path.join(self.destination, 'fsc_'+str(i)+'_even.em')),
                                 Reference(os.path.join(self.destination, 'fsc_'+str(i)+'_odd.em'))]
                
                # low pass filter the reference and write it to the disk
                filtered = lowpassFilter(average, ceil(resolutionBand), ceil(resolutionBand)/10)
                filtered_ref_name = os.path.join(self.destination, 'average_iter'+str(i)+'_res'+str(current_resolution)+'.em')
                filtered[0].write(filtered_ref_name)
                
                # if the position/orientation is not improved, break it
                
                # change the frequency to a higher value
                new_freq = int(ceil(resolutionBand))+1
                if new_freq <= old_freq:
                    if job.adaptive_res is not False: # two different strategies
                        print(self.node_name + ': Determined resolution gets worse. Include additional %f percent frequency to be aligned!' % job.adaptive_res)
                        new_freq = int((1+job.adaptive_res)*old_freq)
                    else: # always increase by 1
                        print(self.node_name + ': Determined resolution gets worse. Increase the frequency to be aligned by 1!')
                        new_freq = old_freq+1
                        old_freq = new_freq
                else:
                    old_freq = new_freq
                if new_freq >= number_bands:
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
        from pytom.lib.frm import frm_align
        from pytom.basic.structures import Shift, Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        
        while True:
            # get the job
            try:
                job = self.get_job()
            except:
                if verbose:
                    print(self.node_name + ': end')
                break # get some non-job message, break it
            
            if verbose:
                prog = FixedProgBar(0, len(job.particleList)-1, self.node_name+':')
                i = 0
            
            ref = job.reference[0].getVolume()
            # run the job
            for p in job.particleList:
                if verbose:
                    prog.update(i)
                    i += 1
                v = p.getVolume()
                
                pos, angle, score = frm_align(v, p.getWedge(), ref, None, job.bw_range, job.freq, job.peak_offset, job.mask.getVolume())
                    
                p.setShift(Shift([pos[0]-v.size_x()/2, pos[1]-v.size_y()/2, pos[2]-v.size_z()/2]))
                p.setRotation(Rotation(angle))
                p.setScore(FRMScore(score))
                
            # average the particle list
            name_prefix = os.path.join(self.destination, self.node_name+'_'+str(job.max_iter))
            self.average_sub_pl(job.particleList, name_prefix, job.weighting)
            
            # send back the result
            self.send_result(FRMResult(name_prefix, job.particleList, self.mpi_id))
        
        pytom_mpi.finalise()
    
    def average_sub_pl(self, pl, name_prefix, weight_average):
        """For worker node.
        """
        pl.average(name_prefix+'.em', progressBar=False, createInfoVolumes=False, _mpiParallel=False, weighting=weight_average)
    
    def retrieve_res_vols(self, name_prefix):
        """For master node, retrieve the sub-averages and do the cleaning.
        """
        from pytom.lib.pytom_volume import read
        pre = read(name_prefix+'-PreWedge.em')
        wedge = read(name_prefix+'-WedgeSumUnscaled.em')
        
        # delete the volumes from disk
        import os
        try:
            os.remove(name_prefix+'-PreWedge.em')
            os.remove(name_prefix+'-WedgeSumUnscaled.em')
            os.remove(name_prefix+'.em')
        except:
            pass
        
        return (pre, wedge)
    
    def create_average(self, pre, wedge):
        """For the master node, create the average according to the pre-wedge and wedge volumes.
        """
        from pytom.lib.pytom_volume import complexDiv, limit
        from pytom.basic.fourier import fft,ifft
        
        limit(wedge, 0.1, 0, 0,0,True,False) # set all the values below the specified value to 0
        
        f_pre = fft(pre)
        r = complexDiv(f_pre, wedge)
        average = ifft(r)
        average.shiftscale(0.0,1/float(average.size_x()*average.size_y()*average.size_z()))
        
        return average
    
    def determine_resolution(self, even, odd, criterion, number_bands, mask, verbose=False):
        """For the master node, determine the resolution.
        """
        from pytom.basic.correlation import fsc, determine_resolution
        
        if not number_bands:
            number_bands = even.size_x()/2
        
        calc_fsc = fsc(even, odd, number_bands, mask, verbose=False)
        if verbose:
            print(self.node_name + ': FSC: ' + str(calc_fsc))
        
        return determine_resolution(calc_fsc, criterion, verbose=False)
    
    def send_job(self, job, dest):
        pytom_mpi.send(str(job), dest)
    
    def get_job(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        job = FRMJob()
        job.fromStr(mpi_msgString)
        
        return job
    
    def send_result(self, result):
        pytom_mpi.send(str(result), 0)
    
    def get_result(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        result = FRMResult()
        result.fromStr(mpi_msgString)
        
        return result
    
    def distribute_job(self, job, verbose=False):
        if len(job.reference) != 2:
            raise Exception('Number of provided references must be 2!')
        
        pls = job.particleList.splitByClass()
        assert len(pls) == 2
        
        self.assignment = {} # note the assignment down
        
        # the first half set
        pl = pls[0]
        label = int(pl[0].getClass())
        num_workers = self.num_workers/2
        particlesPerNode = int(len(pl)/num_workers)
        residual = len(pl)-particlesPerNode*num_workers
        start_idx = 0
        for i in range(1, self.num_workers/2+1):
            if i < residual+1: # since i starts from 1
                l = particlesPerNode+1
            else:
                l = particlesPerNode
            
            subPL = pl[start_idx : start_idx+l]
            start_idx += l
            
            subJob = FRMJob(subPL, [job.reference[label]], job.mask, job.peak_offset, job.sampleInformation, job.bw_range, job.freq, job.destination, job.max_iter, job.r_score, job.weighting)
            self.send_job(subJob, i)
            
            self.assignment[i] = label
            
            if verbose:
                print(self.node_name + ': distributed %d particles of label %d to node %d' % (len(subPL), label, i))
        
        # the second half set
        pl = pls[1]
        label = int(pl[0].getClass())
        num_workers = self.num_workers - num_workers # get the rest
        particlesPerNode = int(len(pl)/num_workers)
        residual = len(pl)-particlesPerNode*num_workers
        start_idx = 0
        for i in range(self.num_workers/2+1, self.num_workers+1):
            if i < residual+self.num_workers/2+1: # since i starts from self.num_workers/2+1
                l = particlesPerNode+1
            else:
                l = particlesPerNode
            
            subPL = pl[start_idx : start_idx+l]
            start_idx += l
            
            subJob = FRMJob(subPL, [job.reference[label]], job.mask, job.peak_offset, job.sampleInformation, job.bw_range, job.freq, job.destination, job.max_iter, job.r_score, job.weighting)
            self.send_job(subJob, i)
            
            self.assignment[i] = label
            
            if verbose:
                print(self.node_name + ': distributed %d particles of label %d to node %d' % (len(subPL), label, i))

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Subtomogram alignment by Fast Rotational Matching.',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-j'], 'Job xml file.', True, True),
                                    ScriptOption(['-v'], 'Verbose mode.', False, False),
                                    ScriptOption(['--help'], 'Help info.', False, False)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    
    try:
        job_filename, verbose, bHelp = parse_script_options(sys.argv[1:], helper)     
    except:
        sys.exit()
        
    if bHelp is True:
        print(helper)
        sys.exit()
    
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
