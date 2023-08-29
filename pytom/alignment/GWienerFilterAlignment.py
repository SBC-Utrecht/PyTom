'''
Gold standard Wiener Filter alignment
Created on Nov 6, 2012

@author: yuxiangchen
'''


from pytom.basic.structures import PyTomClass
import pytom.lib.pytom_mpi as pytom_mpi
from pytom.alignment.FRMAlignment import FRMJob, FRMScore, FRMResult, FRMWorker

class ParticleListPair(PyTomClass):
    def __init__(self, phase_flip_pl, ctf_conv_pl, ctf_sqr, snr=10):
        self.phase_flip_pl = phase_flip_pl # only store the file name of this particle list
        self.ctf_conv_pl = ctf_conv_pl # only store the directory name, the names of the particles assumed to be the same with the previous ones
        self.ctf_sqr = ctf_sqr # store the file name of the 3D volume
        self.snr = snr
        
        self.phase_flip_pl_obj = None # the actual ParticleList object is stored here
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "ParticleListPair":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('ParticleListPair')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a ParticleListPair.")
            jobDescription = jobDescription[0]
        
        self.phase_flip_pl = jobDescription.get('PhaseFlippedParticleList')
        self.ctf_conv_pl = jobDescription.get('CTFConvolutedParticleList')
        self.ctf_sqr = jobDescription.get('CTFSquared')
        self.snr = float(jobDescription.get('SNR'))
    
    def toXML(self):
        from lxml import etree

        jobElement = etree.Element("ParticleListPair")
        
        jobElement.set("PhaseFlippedParticleList", str(self.phase_flip_pl))
        jobElement.set("CTFConvolutedParticleList", str(self.ctf_conv_pl))
        jobElement.set("CTFSquared", str(self.ctf_sqr))
        jobElement.set("SNR", str(self.snr))
        
        return jobElement
    
    def check(self):
        from pytom.tools.files import checkFileExists
        self.get_phase_flip_pl().check()
        self.get_ctf_conv_pl().check()
        if not checkFileExists(self.ctf_sqr):
            raise Exception('File: %s does not exist!' % self.ctf_sqr)
    
    def set_phase_flip_pl(self, pl):
        self.phase_flip_pl_obj = pl
    
    def get_phase_flip_pl(self):
        if self.phase_flip_pl_obj is None:
            from pytom.basic.structures import ParticleList
            pl = ParticleList('.')
            pl.fromXMLFile(self.phase_flip_pl)
            self.phase_flip_pl_obj = pl
        return  self.phase_flip_pl_obj
    
    def get_ctf_conv_pl(self):
        pl = self.get_phase_flip_pl().copy() # get a copied version of the phase flipped particle list, otherwise will change the original pl
        for p in pl:
            name = p.getFilename()
            p.setFilename(self.ctf_conv_pl+'/'+name.split('/')[-1]) # remain the name unchanged but append the correct path
        return pl
    
    def get_ctf_sqr_vol(self):
        from pytom.lib.pytom_volume import read
        v = read(self.ctf_sqr)
        return v

class ParticleListSet(PyTomClass):
    def __init__(self, pairs=[]):
        self.pairs = pairs
    
    def append(self, pair):
        self.pairs.append(pair)
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "ParticleListSet":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('ParticleListSet')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a ParticleListSet.")
            jobDescription = jobDescription[0]
        
        self.pairs = []
        pairs = jobDescription.xpath('ParticleListPair')
    
        if len(pairs) > 0:
            for p in pairs:
                pp = ParticleListPair('', '', '')
                pp.fromXML(p)
                self.pairs.append(pp)
    
    def toXML(self):
        from lxml import etree

        jobElement = etree.Element("ParticleListSet")
        
        if len(self.pairs) > 0:
            for p in self.pairs:
                jobElement.append(p.toXML())
        
        return jobElement


class MultiDefocusJob(FRMJob):
    """For the entry of the whole procedure.
    """
    def fromXML(self, xmlObj): # only rewrite this function
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "FRMJob": # the name is not changed here!
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('FRMJob')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a FRMJob.")
            jobDescription = jobDescription[0]
        
        from pytom.basic.structures import Reference, Mask, SampleInformation, MultiSymmetries
        
        particleList_element = jobDescription.xpath('ParticleListSet')[0]
        pl = ParticleListSet()
        pl.fromXML(particleList_element)
        self.particleList = pl # here i still use the original name!
        
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
    
    def check(self):
        for p in self.particleList.pairs:
            p.check()


class CTFCorrectionJob(PyTomClass):
    """For the worker to actually run.
    """
    def __init__(self, pl=None, ctf_conv_pl=None, ref=None, mask=None, peak_offset=0, sample_info=None, bw_range=None, freq=None, dest='.', max_iter=10, r_score=False, weighting=False, bfactor=None, sum_ctf_sqr=None):
        self.particleList = pl
        self.ctf_conv_pl = ctf_conv_pl
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
        self.sum_ctf_sqr = sum_ctf_sqr
    
    def fromXML(self, xmlObj):
        from lxml.etree import _Element
        
        if xmlObj.__class__ != _Element:
            raise Exception('You must provide a valid XML object.')
        
        if xmlObj.tag == "CTFCorrectionJob":
            jobDescription = xmlObj
        else:  
            jobDescription = xmlObj.xpath('CTFCorrectionJob')
            if len(jobDescription) == 0:
                raise Exception("This XML is not a CTFCorrectionJob.")
            jobDescription = jobDescription[0]
        
        from pytom.basic.structures import ParticleList, Reference, Mask, SampleInformation
        
        particleList_element = jobDescription.xpath('ParticleList')[0]
        pl = ParticleList('.')
        pl.fromXML(particleList_element)
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
            self._sampleInformation = SampleInformation()
        
        self.ctf_conv_pl = jobDescription.get('CTFConvolutedParticleList')
        self.peak_offset = int(jobDescription.get('PeakOffset'))
        self.bw_range = [int(i) for i in jobDescription.get('BandwidthRange')[1:-1].split(',')]
        self.freq = int(jobDescription.get('Frequency'))
        self.destination = jobDescription.get('Destination')
        self.max_iter = int(jobDescription.get('MaxIterations'))
        self.r_score = jobDescription.get('RScore')=='True'
        self.weighting = jobDescription.get('WeightedAverage')=='True'
        self.bfactor = jobDescription.get('BFactor')
        self.sum_ctf_sqr = jobDescription.get('CTFSquared')
    
    def toXML(self):
        from lxml import etree

        jobElement = etree.Element("CTFCorrectionJob")
        jobElement.append(self.particleList.toXML())
        for ref in self.reference:
            jobElement.append(ref.toXML())
        jobElement.append(self.mask.toXML())
        jobElement.append(self.sampleInformation.toXML())
        jobElement.set("CTFConvolutedParticleList", self.ctf_conv_pl)
        jobElement.set("PeakOffset", str(self.peak_offset))
        jobElement.set("BandwidthRange", str(self.bw_range))
        jobElement.set("Frequency", str(self.freq))
        jobElement.set("Destination", self.destination)
        jobElement.set("MaxIterations", str(self.max_iter))
        jobElement.set("RScore", str(self.r_score))
        jobElement.set("WeightedAverage", str(self.weighting))
        jobElement.set("BFactor", str(self.bfactor))
        jobElement.set("CTFSquared", str(self.sum_ctf_sqr))
        
        return jobElement

class MultiDefocusWorker(FRMWorker):
    # rewrite several functions listed below
    
    def start(self, job, verbose=False):
        if self.mpi_id == 0:
            from pytom.basic.structures import ParticleList, Reference
            from pytom.basic.resolution import bandToAngstrom
            from pytom.basic.filter import lowpassFilter
            from math import ceil
            from pytom.basic.fourier import convolute
            from pytom.lib.pytom_volume import vol, power, read
            
            # randomly split the particle list into 2 half sets
            import numpy as np
            num_pairs = len(job.particleList.pairs)
            for i in range(num_pairs):
                # randomize the class labels to indicate the two half sets
                pl = job.particleList.pairs[i].get_phase_flip_pl()
                n = len(pl)
                labels = np.random.randint(2, size=(n,))
                print(self.node_name + ': Number of 1st half set:', n-np.sum(labels), 'Number of 2nd half set:', np.sum(labels))
                for j in range(n):
                    p = pl[j]
                    p.setClass(labels[j])
            
            new_reference = job.reference
            old_freq = job.freq
            new_freq = job.freq
            # main node
            for i in range(job.max_iter):
                if verbose:
                    print(self.node_name + ': starting iteration %d ...' % i)
                
                # construct a new job by updating the reference and the frequency
                # here the job.particleList is actually ParticleListSet
                new_job = MultiDefocusJob(job.particleList, new_reference, job.mask, job.peak_offset, job.sampleInformation, job.bw_range, new_freq, job.destination, job.max_iter-i, job.r_score, job.weighting, job.bfactor)
                
                # distribute it
                num_all_particles = self.distribute_job(new_job, verbose)
                
                # calculate the denominator
                sum_ctf_squared = None
                for pair in job.particleList.pairs:
                    if sum_ctf_squared is None:
                        sum_ctf_squared = pair.get_ctf_sqr_vol() * pair.snr
                    else:
                        sum_ctf_squared += pair.get_ctf_sqr_vol() * pair.snr
                
                # get the result back
                all_even_pre = None
                all_even_wedge = None
                all_odd_pre = None
                all_odd_wedge = None
                pls = []
                for j in range(len(job.particleList.pairs)):
                    pls.append(ParticleList())
                
                for j in range(self.num_workers):
                    result = self.get_result()
                    
                    pair_id = self.assignment[result.worker_id]
                    pair = job.particleList.pairs[pair_id]
                    
                    pl = pls[pair_id]
                    pl += result.pl
                    even_pre, even_wedge, odd_pre, odd_wedge = self.retrieve_res_vols(result.name)
                    
                    if all_even_pre:
                        all_even_pre += even_pre * pair.snr
                        all_even_wedge += even_wedge
                        all_odd_pre += odd_pre * pair.snr
                        all_odd_wedge += odd_wedge
                    else:
                        all_even_pre = even_pre * pair.snr
                        all_even_wedge = even_wedge
                        all_odd_pre = odd_pre * pair.snr
                        all_odd_wedge = odd_wedge
                
                # write the new particle list to the disk
                for j in range(len(job.particleList.pairs)):
                    pls[j].toXMLFile('aligned_pl'+str(j)+'_iter'+str(i)+'.xml')
                
                # correct for the number of particles in wiener filter
                sum_ctf_squared = sum_ctf_squared/num_all_particles
#                all_even_pre = all_even_pre/(num_all_particles/2)
#                all_odd_pre = all_odd_pre/(num_all_particles/2)
                
                # bfactor
                if job.bfactor and job.bfactor != 'None':
#                    bfactor_kernel = create_bfactor_vol(sum_ctf_squared.size_x(), job.sampleInformation.getPixelSize(), job.bfactor)
                    bfactor_kernel = read(job.bfactor)
                    bfactor_kernel_sqr = vol(bfactor_kernel)
                    power(bfactor_kernel_sqr, 2)
                    all_even_pre = convolute(all_even_pre, bfactor_kernel, True)
                    all_odd_pre = convolute(all_odd_pre, bfactor_kernel, True)
                    sum_ctf_squared = sum_ctf_squared*bfactor_kernel_sqr

                
                # create averages of two sets
                if verbose:
                    print(self.node_name + ': determining the resolution ...')
                even = self.create_average(all_even_pre, sum_ctf_squared, all_even_wedge) # assume that the CTF sum is the same for the even and odd
                odd = self.create_average(all_odd_pre, sum_ctf_squared, all_odd_wedge)
                
                # determine the transformation between even and odd
                # here we assume the wedge from both sets are fully sampled
                from pytom.lib.frm import frm_align
                pos, angle, score = frm_align(odd, None, even, None, job.bw_range, new_freq, job.peak_offset)
                print(self.node_name + ': transform of even set to match the odd set - shift: '+str(pos)+' rotation: '+str(angle))
                
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
                
                transformSpline(all_odd_pre,transformed_odd_pre,-angle[1],-angle[0],-angle[2],
                                int(odd.size_x()/2), int(odd.size_y()/2), int(odd.size_z()/2),
                                -(pos[0]-odd.size_x()/2),-(pos[1]-odd.size_y()/2),-(pos[2]-odd.size_z()/2),0,0,0)
                odd_weight.rotate(-angle[1],-angle[0],-angle[2])
                transformed_odd_wedge = odd_weight.getWeightVolume(True)
                transformSpline(odd,transformed_odd,-angle[1],-angle[0],-angle[2],
                                int(odd.size_x()/2), int(odd.size_y()/2), int(odd.size_z()/2),
                                -(pos[0]-odd.size_x()/2),-(pos[1]-odd.size_y()/2),-(pos[2]-odd.size_z()/2),0,0,0)
                
                all_odd_pre = transformed_odd_pre
                all_odd_wedge = transformed_odd_wedge
                odd = transformed_odd
                
                # apply symmetries before determine resolution
                # with gold standard you should be careful about applying the symmetry!
                even = job.symmetries.applyToParticle(even)
                odd = job.symmetries.applyToParticle(odd)
                resNyquist, resolutionBand, number_bands = self.determine_resolution(even, odd, job.fsc_criterion, None, job.mask, verbose)
                
                # write the half set to the disk
                even.write('fsc_'+str(i)+'_even.em')
                odd.write('fsc_'+str(i)+'_odd.em')
                
                current_resolution = bandToAngstrom(resolutionBand, job.sampleInformation.getPixelSize(), number_bands, 1)
                if verbose:
                    print(self.node_name + ': current resolution ' + str(current_resolution), resNyquist)
                
                # create new average
                all_even_pre += all_odd_pre
                all_even_wedge += all_odd_wedge
#                all_even_pre = all_even_pre/2 # correct for the number of particles in wiener filter
                average = self.create_average(all_even_pre, sum_ctf_squared, all_even_wedge)
                
                # apply symmetries
                average = job.symmetries.applyToParticle(average)
                
                # filter average to resolution and update the new reference
                average_name = 'average_iter'+str(i)+'.em'
                average.write(average_name)
                
                # update the references
                new_reference = [Reference('fsc_'+str(i)+'_even.em'), Reference('fsc_'+str(i)+'_odd.em')]
                
                # low pass filter the reference and write it to the disk
                filtered = lowpassFilter(average, ceil(resolutionBand), ceil(resolutionBand)/10)
                filtered_ref_name = 'average_iter'+str(i)+'_res'+str(current_resolution)+'.em'
                filtered[0].write(filtered_ref_name)
                
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
                    print(self.node_name + ': Determined frequency too high. Terminate!')
                    break
                
                if verbose:
                    print(self.node_name + ': change the frequency to ' + str(new_freq))
            
            # send end signal to other nodes and terminate itself
            self.end(verbose)
        else:
            # other nodes
            self.run(verbose)
    
    def run(self, verbose=False):
        from pytom.lib.frm import frm_align
        from pytom.basic.structures import Shift, Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.basic.fourier import convolute
        from pytom.lib.pytom_volume import read, power
        
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
            
            ref = []
            ref.append(job.reference[0].getVolume())
            ref.append(job.reference[1].getVolume())
            
            # convolute with the approximation of the CTF
            if job.sum_ctf_sqr:
                ctf = read(job.sum_ctf_sqr)
                power(ctf, 0.5) # the number of CTFs should not matter, should it?
                ref0 = ref[0]
                ref1 = ref[1]
                ref0 = convolute(ref0, ctf, True)
                ref1 = convolute(ref1, ctf, True)
                ref = [ref0, ref1]
            
            if job.bfactor and job.bfactor != 'None':
#                restore_kernel = create_bfactor_restore_vol(ref.size_x(), job.sampleInformation.getPixelSize(), job.bfactor)
                from pytom.lib.pytom_volume import vol, read
                bfactor_kernel = read(job.bfactor)
                unit = vol(bfactor_kernel)
                unit.setAll(1)
                restore_kernel = unit/bfactor_kernel
            
            # run the job
            for p in job.particleList:
                if verbose:
                    prog.update(i)
                    i += 1
                v = p.getVolume()
                
#                if weights is None: # create the weights according to the bfactor
#                    if job.bfactor == 0:
#                        weights = [1 for k in xrange(job.freq)]
#                    else:
#                        restore_fnc = create_bfactor_restore_fnc(ref.size_x(), job.sampleInformation.getPixelSize(), job.bfactor)
#                        # cut out the corresponding part and square it to get the weights!
#                        weights = restore_fnc[1:job.freq+1]**2

                if job.bfactor and job.bfactor != 'None':
                    v = convolute(v, restore_kernel, True) # if bfactor is set, restore it
                
                pos, angle, score = frm_align(v, p.getWedge(), ref[int(p.getClass())], None, job.bw_range, job.freq, job.peak_offset, job.mask.getVolume())
                
                p.setShift(Shift([pos[0]-v.size_x()/2, pos[1]-v.size_y()/2, pos[2]-v.size_z()/2]))
                p.setRotation(Rotation(angle))
                p.setScore(FRMScore(score))
                
            # average the particle list
            name_prefix = self.node_name+'_'+str(job.max_iter)
            pair = ParticleListPair('', job.ctf_conv_pl, None, None)
            pair.set_phase_flip_pl(job.particleList)
            self.average_sub_pl(pair.get_ctf_conv_pl(), name_prefix) # operate on the CTF convoluted projection!
            
            # send back the result
            self.send_result(FRMResult(name_prefix, job.particleList, self.mpi_id))
        
        pytom_mpi.finalise()
    
    def get_job(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        job = CTFCorrectionJob()
        job.fromStr(mpi_msgString)
        
        return job
    
    def average_sub_pl(self, pl, name_prefix):
        """For worker node, this function has been rewritten.
        """
        from pytom.basic.structures import ParticleList
        even = ParticleList('.')
        odd = ParticleList('.')
        
        for p in pl:
            i = int(p.getClass())
            if i%2 == 0:
                even.append(p)
            else:
                odd.append(p)
        
        assert len(even) > 0, "average_sub_pl: length even particle list == 0 :("
        assert len(odd)  > 0, "average_sub_pl: length odd particle list == 0 :("
        
        # in some rare cases this would fail
        # that is: this worker only get one class, in this case, one of the pl will be null
        # just restart the job
        self.sum_sub_pl(even, name_prefix+'even')
        self.sum_sub_pl(odd, name_prefix+'odd')
    
    def sum_sub_pl(self, pl, name_prefix):
        """This is a sub-routine for average_sub_pl.
        """
        from pytom.lib.pytom_volume import vol
        from pytom.lib.pytom_volume import transformSpline as transform
        from pytom.basic.normalise import mean0std1
        
        result = None
        wedgeSum = None
        for p in pl:
            particle = p.getVolume()
            mean0std1(particle)
            wedgeInfo = p.getWedge()
            
            if result is None:
                size_x = particle.size_x() 
                size_y = particle.size_y()
                size_z = particle.size_z()
                
                newParticle = vol(size_x,size_y,size_z)
                # make consistent for python3
                centerX = size_x//2 
                centerY = size_y//2 
                centerZ = size_z//2 
                
                result = vol(size_x,size_y,size_z)
                result.setAll(0)
                
                wedgeSum = wedgeInfo.returnWedgeVolume(size_x,size_y,size_z)
                wedgeSum.setAll(0)
            
            # create wedge weighting
            rotation = p.getRotation()
            
            wedge = wedgeInfo.returnWedgeVolume(size_x,size_y,size_z,False,rotation.invert())
            wedgeSum = wedgeSum + wedge
            
            # shift and rotate particle
            shift = p.getShift()
            newParticle.setAll(0)
            transform(particle,newParticle,-rotation[1],-rotation[0],-rotation[2],
                    centerX,centerY,centerZ,-shift[0],-shift[1],-shift[2],0,0,0)
            
            result = result + newParticle
        
        # write them back to disk
        result.write(name_prefix+'-PreWedge.em')
        wedgeSum.write(name_prefix+'-WedgeSumUnscaled.em')
    
    def create_average(self, sum_ctf_conv, sum_ctf_squared, wedge_weight):
        """For the master node, this function is rewritten.
        """
        from pytom.lib.pytom_volume import complexDiv, fullToReduced, limit
        from pytom.basic.fourier import fft, ifft, ftshift
        from pytom.basic.normalise import mean0std1
        
        limit(wedge_weight, 0.1, 0, 0,0,True,False) # set all the values below the specified value to 0
        
        # for mask out the outside area
#        mask = vol(sum_ctf_conv)
#        mask.setAll(0)
#        initSphere(mask, sum_ctf_conv.size_x()/2-1, 0,0, sum_ctf_conv.size_x()/2, sum_ctf_conv.size_x()/2, sum_ctf_conv.size_x()/2)
#        mask = fullToReduced(ftshift(mask, inplace=False))
        
        # Wiener filter
        numerator = fft(sum_ctf_conv)
        sum_ctf_squared = fullToReduced(ftshift(sum_ctf_squared, inplace=False))
        denominator = (sum_ctf_squared+1)*wedge_weight
        r = complexDiv(numerator, denominator)
#        average = ifft(complexRealMult(r, mask))
        average = ifft(r)
        average.shiftscale(0.0,1/float(average.size_x()*average.size_y()*average.size_z()))
        
        # nomalize the average
        try:
            average = mean0std1(average, True)
        except:
            average *= 1000 # in case the average volume is too small to normalize
            average = mean0std1(average, True)
        
        return average
    
    def distribute_job(self, job, verbose=False):
        # get the number of different particle list pairs
        num_pairs = len(job.particleList.pairs)
        num_all_particles = 0
        for i in range(num_pairs):
            num_all_particles += len(job.particleList.pairs[i].get_phase_flip_pl())
        
        self.assignment = {} # dict for retrieving back the corresponding result
        start_worker_idx = 1
        for i in range(num_pairs):
            # nodes available for this pair
            if i != num_pairs-1:
                num_workers = int(float(len(job.particleList.pairs[i].get_phase_flip_pl()))/num_all_particles*self.num_workers)
            else: # the last one takes all the rest nodes
                num_workers = self.num_workers - start_worker_idx + 1
            
            if verbose:
                print(self.node_name + ': assign pair %d to workers from %d to %d' % (i, start_worker_idx, start_worker_idx+num_workers-1))
            
            pl = job.particleList.pairs[i].get_phase_flip_pl()
            
            sub_job = CTFCorrectionJob(pl, job.particleList.pairs[i].ctf_conv_pl, job.reference, job.mask, job.peak_offset, job.sampleInformation, job.bw_range, job.freq, job.destination, job.max_iter, job.r_score, job.weighting, job.bfactor, job.particleList.pairs[i].ctf_sqr)
            self.distribute_pl(sub_job, i, num_workers, start_worker_idx, verbose)
            
            start_worker_idx += num_workers
        
        return num_all_particles
    
    def distribute_pl(self, job, pair_id, num_workers, start_worker_idx, verbose=False):
        """This is a sub-routine used in distribute_job.
        """
        n = len(job.particleList)
        particlesPerNode = int(n/num_workers)
        residual = n-particlesPerNode*num_workers
        start_idx = 0
        for i in range(start_worker_idx, num_workers+start_worker_idx):
            if i < residual+start_worker_idx:
                l = particlesPerNode+1
            else:
                l = particlesPerNode
            
            subPL = job.particleList[start_idx : start_idx+l]
            start_idx += l
            
            subJob = CTFCorrectionJob(subPL, job.ctf_conv_pl, job.reference, job.mask, job.peak_offset, job.sampleInformation, job.bw_range, job.freq, job.destination, job.max_iter, job.r_score, job.weighting, job.bfactor, job.sum_ctf_sqr)
            self.send_job(subJob, i)
            
            self.assignment[i] = pair_id
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
    job = MultiDefocusJob()
    job.fromXMLFile(job_filename)
    job.check()
    
    worker = MultiDefocusWorker()
    
    if verbose:
        from pytom.tools.timing import Timing
        t = Timing()
        t.start()
    
    # start it
    worker.start(job, verbose)
    
    if verbose:
        print('Overall execution time: %f s.' % t.end())
