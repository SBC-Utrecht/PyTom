'''
Created on May 7, 2012

@author: yuxiangchen
'''

import pytom_mpi
import cPickle

class CMWorker():
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
            import numpy as np
            pl_filename = job["ParticleList"]
            from pytom.basic.structures import ParticleList
            pl = ParticleList('.')
            pl.fromXMLFile(pl_filename)
            n = len(pl)
            correlation_matrix = np.ones((n, n))

            # distribute the job
            self.distribute_job(job, verbose)

            # fill in the diagnal line
            for i in xrange(n):
                correlation_matrix[i][i] = 1

            # gather the results
            for i in xrange(self.num_workers):
                result = self.get_result()
                pairs = result.keys()
                for pair in pairs:
                    correlation_matrix[pair[0]][pair[1]] = result[pair]
                    correlation_matrix[pair[1]][pair[0]] = result[pair]
            
            # write the correlation matrix to the disk
            np.savetxt('correlation_matrix.csv', correlation_matrix, delimiter=',')

            # send end signal to other nodes and terminate itself
            self.end(verbose)
        else:
            # other nodes
            self.run(verbose)
    
    def end(self, verbose=False):
        if verbose == True:
            print self.node_name + ': sending end messages to others'
        
        from pytom.parallel.messages import StatusMessage
        
        mpi_numberNodes = pytom_mpi.size()
        mpi_myid = pytom_mpi.rank()
        
        for i in range(1, mpi_numberNodes):
            msg = StatusMessage(str(mpi_myid),str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg),i)
        
        pytom_mpi.finalise()
    
    def run(self, verbose=False):
        from pytom_volume import read, sum
        from pytom.basic.filter import lowpassFilter
        from pytom.basic.correlation import nxcc
        from pytom.basic.structures import Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        
        while True:
            # get the job
            job = self.get_job()
            
            try:
                pairs = job["Pairs"]
                pl_filename = job["ParticleList"]
            except:
                if verbose:
                    print self.node_name + ': end'
                break # get some non-job message, break it

            from pytom.basic.structures import ParticleList
            pl = ParticleList('.')
            pl.fromXMLFile(pl_filename)

            if verbose:
                prog = FixedProgBar(0, len(pairs)-1, self.node_name+':')
                i = 0

            # run the job
            result = {}
            last_filename = None
            binning = job["Binning"]
            mask = read(job["Mask"], 0,0,0,0,0,0,0,0,0, binning,binning,binning)
            for pair in pairs:
                if verbose:
                    prog.update(i)
                    i += 1
                g = pl[pair[0]]
                f = pl[pair[1]]
                vf = f.getTransformedVolume(binning)
                wf = f.getWedge().getWedgeObject()
                wf_rotation = f.getRotation().invert()
                # wf.setRotation(Rotation(-rotation[1],-rotation[0],-rotation[2]))
                # wf_vol = wf.returnWedgeVolume(vf.sizeX(), vf.sizeY(), vf.sizeZ(), True, -rotation[1],-rotation[0],-rotation[2])
                vf = lowpassFilter(vf, job["Frequency"], 0)[0]

                if g.getFilename() != last_filename:
                    vg = g.getTransformedVolume(binning)
                    wg = g.getWedge().getWedgeObject()
                    wg_rotation = g.getRotation().invert()
                    # wg.setRotation(Rotation(-rotation[1],-rotation[0],-rotation[2]))
                    # wg_vol = wg.returnWedgeVolume(vg.sizeX(), vg.sizeY(), vg.sizeZ(), True, -rotation[1],-rotation[0],-rotation[2])
                    vg = lowpassFilter(vg, job["Frequency"], 0)[0]

                    last_filename = g.getFilename()

                score = nxcc(wg.apply(vf, wg_rotation), wf.apply(vg, wf_rotation), mask)
                # overlapped_wedge_vol = wf_vol * wg_vol
                # scaling = float(overlapped_wedge_vol.numelem())/sum(overlapped_wedge_vol)
                # score *= scaling

                result[pair] = score
            
            # send back the result
            self.send_result(result)
        
        pytom_mpi.finalise()

    
    def send_job(self, job, dest):
        pytom_mpi.send(cPickle.dumps(job), dest)
    
    def get_job(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        try:
            job = cPickle.loads(mpi_msgString)
        except:
            return None
        
        return job
    
    def send_result(self, result):
        pytom_mpi.send(cPickle.dumps(result), 0)
    
    def get_result(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        result = cPickle.loads(mpi_msgString)
        
        return result
    
    def distribute_job(self, job, verbose=False):
        pl_filename = job["ParticleList"]
        from pytom.basic.structures import ParticleList
        pl = ParticleList('.')
        pl.fromXMLFile(pl_filename)
        nn = len(pl)
        all_pairs = []
        for i in xrange(nn):
            for j in xrange(i+1, nn):
                all_pairs.append((i, j))
        n = len(all_pairs)
        particlesPerNode = int(n/self.num_workers)
        residual = n-particlesPerNode*self.num_workers
        start_idx = 0
        for i in xrange(1, self.num_workers+1):
            if i < residual+1: # since i starts from 1
                l = particlesPerNode+1
            else:
                l = particlesPerNode
            
            sub_pairs = all_pairs[start_idx : start_idx+l]
            start_idx += l
            
            # construct the job
            sub_job = {}
            sub_job["Pairs"] = sub_pairs
            sub_job["ParticleList"] = pl_filename
            sub_job["Mask"] = job["Mask"]
            sub_job["Frequency"] = job["Frequency"]
            sub_job["Binning"] = job["Binning"]
            self.send_job(sub_job, i)
            
            if verbose:
                print self.node_name + ': distributed %d particles to node %d' % (len(sub_pairs), i)

if __name__ == '__main__':
    # parse command line arguments
    import sys
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    
    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Calculate the correlation matrix, given the aligned particle list. The result will be written to the disk named correlation_matrix.csv.',
                          authors='Yuxiang Chen',
                          options= [ScriptOption(['-p'], 'Aligned particle list file.', True, False),
                                    ScriptOption(['-m'], 'Mask.', True, False),
                                    ScriptOption(['-f'], 'Frequency (after binning).', True, False),
                                    ScriptOption(['-b'], 'Binning factor.', True, True),
                                    ScriptOption(['-v'], 'Verbose mode.', False, True),
                                    ScriptOption(['--help'], 'Help info.', False, True)])
    
    if len(sys.argv) == 1:
        print helper
        sys.exit()
    
    try:
        pl_filename, mask_filename, freq, binning, verbose, bHelp = parse_script_options(sys.argv[1:], helper)     
    except:
        raise
        sys.exit()
    
    if bHelp is True:
        print helper
        sys.exit()
    
    if verbose:
        from pytom.tools.timing import Timing
        t = Timing()
        t.start()
    
    # construct the job
    freq = int(freq)
    if not binning:
        binning = 1
    else:
        binning = int(binning)
    job = {}
    job["ParticleList"] = pl_filename
    job["Mask"] = mask_filename
    job["Frequency"] = freq
    job["Binning"] = binning

    worker = CMWorker()
    worker.start(job, verbose)
    
    if verbose:
        print 'Overall execution time: %f s.' % t.end()
