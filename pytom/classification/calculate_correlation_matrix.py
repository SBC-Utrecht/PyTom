'''
Created on May 7, 2012

@author: yuxiangchen
'''

import pytom_mpi
import pickle
import os
import numpy
from pytom.gpu.initialize import xp



class CMWorkerGPU():
    def __init__(self):
        if not pytom_mpi.isInitialised():
            pytom_mpi.init()

        self.mpi_id = pytom_mpi.rank()
        self.num_workers = pytom_mpi.size() - 1
        self.node_name = 'node_' + str(self.mpi_id)

        if self.num_workers < 1:
            raise RuntimeError("Not enough nodes to parallelize the job!")

    def start(self, job, verbose=False):
        outdir = job["outdir"]
        if self.mpi_id == 0:
            import numpy as np
            pl_filename = job["ParticleList"]
            from pytom.basic.structures import ParticleList
            pl = ParticleList('.')
            pl.fromXMLFile(pl_filename)
            n = len(pl)
            correlation_matrix = np.zeros((n, n))

            # distribute the job
            self.distribute_job(job, verbose)

            # fill in the diagnal line
            for i in range(n):
                correlation_matrix[i][i] = 1

            # gather the results
            for i in range(self.num_workers):
                aa = np.array(self.get_result())
                correlation_matrix += aa + aa.T

            # write the correlation matrix to the disk
            np.savetxt(os.path.join(outdir, 'correlation_matrix.csv'), correlation_matrix, delimiter=',')
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
            msg = StatusMessage(str(mpi_myid), str(i))
            msg.setStatus("End")
            pytom_mpi.send(str(msg), i)

        pytom_mpi.finalise()

    def run(self, verbose=False):
        from pytom_volume import read, sum
        from pytom.basic.filter import lowpassFilter
        from pytom.basic.correlation import nxcc
        from pytom.basic.structures import Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        from pytom.agnostic.transform import fourier_reduced2full
        from pytom.agnostic.io import read
        from pytom.agnostic.tools import create_sphere
        from pytom.gpu.gpuStructures import CCCPlan

        while True:
            # get the job
            job = self.get_job()

            try:
                pl_filename = job["ParticleList"]
                maskname = job["Mask"]
                slists = job["Pairs"]
                nlists = job["NeededParticles"]
                freq= job["Frequency"]
                max_num_load = job["MaxNumPartPerGPU"]
                device = f'gpu:{job["GpuID"]}'
                binning = job['Binning']

            except:
                if verbose:
                    print(self.node_name + ': end')
                break  # get some non-job message, break it

            from pytom.basic.structures import ParticleList
            pl = ParticleList('.')
            pl.fromXMLFile(pl_filename)

            plan = CCCPlan(pl, maskname, freq, cp=xp, device=device, max_num_part=max_num_load, profile=verbose, binning=binning)

            import time
            num_jobs = len(nlists)
            for nn, (subparts, needed) in enumerate(zip(slists, nlists)):
                tt = time.time()
                m = 0
                for n, curr in enumerate(plan.currently_loaded):
                    if curr in sorted(needed):
                        plan.currently_loaded[m] = plan.currently_loaded[n]
                        plan.formatting[m] = plan.formatting[n]
                        m += 1

                plan.currently_loaded[m:] = -1

                for id in needed:
                    if id in plan.currently_loaded:
                        continue
                    plan.currently_loaded[m] = int(id)
                    plan.store_particle(m, pl[int(id)].getFilename())
                    m += 1

                elapsed = (time.time() - tt) * 1000
                print(f'finished reading for job {nn+1}/{num_jobs} on {device} in \t\t{elapsed/1000:10.3f} sec ({elapsed/max(1,len(subparts)):7.3f})')
                tt = time.time()
                plan.ccc_go(subparts)
                elapsed2 = (time.time()-tt)*1000
                elapsed += elapsed2
                # print(f'finished {len(subparts)} comparisons for job {nn+1}/{num_jobs} on {device} in {elapsed2/1000:10.3f} sec ({elapsed2/max(1,len(subparts)):7.3f})')
                print(f'finished {len(subparts)} comparisons for job {nn+1}/{num_jobs} on {device} in {elapsed/1000:10.3f} sec ({elapsed/max(1,len(subparts)):7.3f})')
            # send back the result

            self.send_result(plan.results.get().tolist())
            del plan

        pytom_mpi.finalise()

    def send_job(self, job, dest):
        pickled = pickle.dumps(job, protocol=0, fix_imports=True).decode('ISO-8859-1')
        pytom_mpi.send(pickled, dest)

    def get_job(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        try:
            job = pickle.loads(mpi_msgString.encode('ISO-8859-1'))
        except:
            return None

        return job

    def send_result(self, result):
        pickled = pickle.dumps(result, protocol=0, fix_imports=True).decode('ISO-8859-1')
        pytom_mpi.send(pickled, 0)

    def get_result(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        result = pickle.loads(mpi_msgString.encode('ISO-8859-1'))

        return result

    def distribute_job(self, job, verbose=False, max_num_load=150):
        pl_filename = job["ParticleList"]
        from pytom.basic.structures import ParticleList
        pl = ParticleList('.')
        pl.fromXMLFile(pl_filename)
        num_particles = len(pl)

        max_num_load = min(max_num_load, num_particles//(len(job['gpuIDs'])+1) if len(job['gpuIDs']) > 1 else max_num_load)

        #max_num_load = 50  # min(int(10 * 1024**3 / (x*y*z*8)),num_particles//((len(gpus)+7)**2))

        needed_list, subparts_list = [], []

        bx, by = numpy.meshgrid(numpy.arange(num_particles), numpy.arange(num_particles))

        all_pairs = numpy.zeros((num_particles, num_particles, 2), dtype=int)
        all_pairs[:, :, 0] = bx
        all_pairs[:, :, 1] = by
        all_pairs[by >= bx] = [-1, -1]

        cy, cx = numpy.meshgrid(numpy.arange(0, num_particles, max_num_load),
                                numpy.arange(0, num_particles, max_num_load))

        coords = numpy.zeros((cx.shape[0], cx.shape[1], 2), dtype=int)
        coords[:, :, 0] = cx
        coords[:, :, 1] = cy
        coords = coords[cx <= cy]

        for i, j in coords:
            bb = all_pairs[i:i + max_num_load, j:j + max_num_load, ::-1]
            pairs = bb[bb[:, :] != numpy.array([-1, -1])]
            pairs = pairs.reshape(pairs.shape[0] // 2, 2)
            sp = pairs.copy()

            needed_list.append(sorted(list(numpy.unique(sp))))
            sp[:, 0] -= i
            if j - i >= max_num_load:
                sp[:, 1] -= j - max_num_load
            if i == j:
                sp[:, 1] -= j
            subparts_list.append(sp.tolist())

        num_comps = numpy.zeros((len(job['gpuIDs'])),dtype=int)
        subpartsGPU = []
        neededGPU = []
        for g in job['gpuIDs']:
            subpartsGPU.append([])
            neededGPU.append([])

        for n, subpart in enumerate(subparts_list):
            subpartsGPU[num_comps.argmin()] += [subpart]
            neededGPU[num_comps.argmin()]+= [needed_list[n]]
            num_comps[num_comps.argmin()] += len(subpart)


        for n, gpuID in enumerate(job['gpuIDs']):
            # construct the job
            sub_job = {}
            sub_job["Pairs"] = subpartsGPU[n]
            sub_job["ParticleList"] = pl_filename
            sub_job["Mask"] = job["Mask"]
            sub_job["Frequency"] = job["Frequency"]
            sub_job["Binning"] = job["Binning"]
            sub_job['GpuID'] = gpuID
            sub_job['NeededParticles'] = neededGPU[n]
            sub_job['MaxNumPartPerGPU'] = max_num_load

            self.send_job(sub_job, n+1)

            if verbose:
                print(self.node_name + ': distributed %d particles to node %d' % (len(subpartsGPU[n]), n))


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
        outdir = job["outdir"]
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
            for i in range(n):
                correlation_matrix[i][i] = 1

            # gather the results
            for i in range(self.num_workers):
                result = self.get_result()
                pairs = list(result.keys())
                for pair in pairs:
                    correlation_matrix[pair[0]][pair[1]] = result[pair]
                    correlation_matrix[pair[1]][pair[0]] = result[pair]
            
            # write the correlation matrix to the disk
            np.savetxt(os.path.join(outdir, 'correlation_matrix.csv'), correlation_matrix, delimiter=',')

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
        from pytom_volume import read, sum
        from pytom.basic.filter import lowpassFilter
        from pytom.basic.correlation import nxcc
        from pytom.basic.structures import Rotation
        from pytom.tools.ProgressBar import FixedProgBar
        import time

        t = time.time()
        while True:
            # get the job
            job = self.get_job()
            
            try:
                pairs = job["Pairs"]
                pl_filename = job["ParticleList"]
            except:
                if verbose:
                    print(self.node_name + ': end')
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
            binning = int(job["Binning"])

            mask = read(job["Mask"], 0, 0, 0, 0, 0, 0, 0, 0, 0, binning, binning, binning)
            
            
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


                score = nxcc( wg.apply(vf, wg_rotation), wf.apply(vg, wf_rotation), mask)
                # overlapped_wedge_vol = wf_vol * wg_vol
                # scaling = float(overlapped_wedge_vol.numelem())/sum(overlapped_wedge_vol)
                # score *= scaling

                result[pair] = score
            
            # send back the result
            self.send_result(result)
        print((time.time()-t)*1000/len(pairs))
        pytom_mpi.finalise()

    
    def send_job(self, job, dest):
        pickled = pickle.dumps(job, protocol=0, fix_imports=True).decode('utf-8')
        pytom_mpi.send(pickled, dest)
    
    def get_job(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        try:
            job = pickle.loads(mpi_msgString.encode('utf-8'))
        except:
            return None
        
        return job
    
    def send_result(self, result):
        pickled = pickle.dumps(result, protocol=0, fix_imports=True).decode('utf-8')
        pytom_mpi.send(pickled, 0)
    
    def get_result(self):
        from pytom.localization.parallel_extract_peaks import getMsgStr
        mpi_msgString = getMsgStr()
        result = pickle.loads(mpi_msgString.encode('utf-8'))
        
        return result
    
    def distribute_job(self, job, verbose=False):
        pl_filename = job["ParticleList"]
        from pytom.basic.structures import ParticleList
        pl = ParticleList('.')
        pl.fromXMLFile(pl_filename)
        nn = len(pl)
        all_pairs = []
        for i in range(nn):
            for j in range(i+1, nn):
                all_pairs.append((i, j))
        n = len(all_pairs)
        particlesPerNode = int(n/self.num_workers)
        residual = n-particlesPerNode*self.num_workers
        start_idx = 0
        for i in range(1, self.num_workers+1):
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
                print(self.node_name + ': distributed %d particles to node %d' % (len(sub_pairs), i))

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
                                    ScriptOption(['-o'], 'Output directory.', True, True),
                                    ScriptOption(['-g', '--gpuID'], "Index or indices of the gpu's one wants to use. CCC can run on multiple gpu's simultaneously. The indices "
                                                                    "of multiple gpu's are separated by a comma (no space). For example 0,2,3,5 **Please note that the number "
                                                                    "of mpi cores should be one more than the number of GPUs you are using.**", True, True),
                                    ScriptOption(['--help'], 'Help info.', False, True)])
    
    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    
    try:
        pl_filename, mask_filename, freq, binning, verbose, outdir, gpuIDs, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(e)
        print(helper)
        sys.exit()
    
    if help is True:
        print(helper)
        sys.exit()
    
    if 1:
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
    job["outdir"] = outdir if outdir else './'
    job['gpuIDs'] = [] if gpuIDs is None else list(map(int, gpuIDs.split(',')))

    if not job['gpuIDs']:
        worker = CMWorker()
        worker.start(job, verbose)
    else:
        worker = CMWorkerGPU()
        worker.start(job, verbose)

    if worker.mpi_id == 0:
        print('Overall execution time: %f s.' % t.end())
