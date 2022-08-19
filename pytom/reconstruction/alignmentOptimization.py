#!/usr/bin/env pytom
'''
Created on jan 2021
@author: Gvd
'''

class FiducialLessAlignment():
    def __init__(self):
        import sys, os, re, subprocess
        from scipy import misc
        import numpy as np
        import scipy.optimize
        from pytom_volume import transformSpline, vol
        from pytom.basic.correlation import xcc, nxcc
        from pytom.basic.transformations import scale, general_transform2d
        from pytom.reconstruction.ccAlign import modelProjCrop, vol_projection, punchHole, cylindricalMask, create_cylindrical_mask, reAlignWeightReconstruct, reAlignWeight
        from pytom.basic.files import read, read_em, write_em, read_em_header, write_em_header
        from pytom_volume import vol, read, subvolume
        from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
        from pytom.reconstruction.TiltAlignmentStructures import TiltSeries
        from pytom.reconstruction.writeAlignedProjections import writeAlignedProjections
        from pytom.gpu.initialize import xp

        self.cp = xp
        self.score = self.singleNXCC

        self.num_threads = 1024

        self.sumMeanStdv = xp.RawKernel(r'''
            __device__ void warpReduce(volatile float* sdata, volatile float* stdv, int tid, int blockSize) {
                if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; stdv[tid] += stdv[tid+32];}
                if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; stdv[tid] += stdv[tid+16];}
                if (blockSize >= 16) {sdata[tid] += sdata[tid +  8]; stdv[tid] += stdv[tid+ 8];}
                if (blockSize >=  8) {sdata[tid] += sdata[tid +  4]; stdv[tid] += stdv[tid+ 4];}
                if (blockSize >=  4) {sdata[tid] += sdata[tid +  2]; stdv[tid] += stdv[tid+ 2];}
                if (blockSize >=  2) {sdata[tid] += sdata[tid +  1]; stdv[tid] += stdv[tid+ 1];} }



            extern "C" __global__ 
            void sumMeanStdv(float *g_idata, float *mask, float *g_mean, float * g_stdv,  int n) {

                __shared__ float mean[1024];
                __shared__ float stdv[1024];
                int blockSize = blockDim.x; 
                unsigned int tid = threadIdx.x;                                                                                                                                        
                int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
                int gridSize = blockSize*gridDim.x*2;                                                                                                                         

                mean[tid] = 0.;                                                                                                                                                       
                stdv[tid] = 0.;

                while (i < n) {
                    mean[tid] += g_idata[i] * mask[i] ; 
                    stdv[tid] += g_idata[i] * g_idata[i] * mask[i];
                    if (i + blockSize < n){
                        mean[tid] += g_idata[i + blockSize] * mask[i+blockSize];
                        stdv[tid] += g_idata[i + blockSize] * g_idata[i + blockSize] * mask[i+blockSize];}
                    i += gridSize;}
                 __syncthreads();                                                                                                                                                       

                if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512]; stdv[tid] += stdv[tid+512];} __syncthreads(); }                                                                                                    
                if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256]; stdv[tid] += stdv[tid+256];} __syncthreads(); }                                                                          
                if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128]; stdv[tid] += stdv[tid+128];} __syncthreads(); }                                                                          
                if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64]; stdv[tid] += stdv[tid+ 64];} __syncthreads(); }                                                                          
                if (tid < 32){ warpReduce(mean, stdv, tid, blockSize);}                                                                                                                                                                      
                if (tid == 0) {g_mean[blockIdx.x] = mean[0]; g_stdv[blockIdx.x] = stdv[0];} 
               __syncthreads(); 

            }''', 'sumMeanStdv')

        self.sumMultNorm = xp.RawKernel(r'''
            __device__ void warpReduce(volatile float* sdata, int tid, int blockSize) {
                if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
                if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
                if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
                if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
                if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
                if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



            extern "C" __global__ 
            void sumMultNorm(float *g_idata, float *h_idata, float g_mean, float g_stdv, float *g_result, int n) {

                __shared__ float mean[1024];
                __shared__ float stdv[1024];
                int blockSize = blockDim.x; 
                unsigned int tid = threadIdx.x;                                                                                                                                        
                int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
                int gridSize = blockSize*gridDim.x*2;                                                                                                                         

                mean[tid] = 0.;                                                                                                                                                       
                

                while (i < n) {
                    mean[tid] += (g_idata[i] - g_mean)/g_stdv * h_idata[i]; 
                    
                    if (i + blockSize < n){
                        mean[tid] += (g_idata[i + blockSize] - g_mean) / g_stdv * h_idata[i+blockSize];
                        }
                    i += gridSize;}
                 __syncthreads();                                                                                                                                                       

                if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
                if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
                if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
                if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
                if (tid < 32){ warpReduce(mean, tid, blockSize);}                                                                                                                                                                      
                if (tid == 0) {g_result[blockIdx.x] = mean[0];} 
               __syncthreads(); 

            }''', 'sumMultNorm')

        self.sumKernel = xp.RawKernel(r'''
                __device__ void warpReduceSum(volatile float* sdata, int tid, int blockSize) {
                    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
                    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
                    if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
                    if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
                    if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
                    if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



                extern "C" __global__ 
                void sumKernel(float *g_idata, float *g_mean,  int n) {

                    __shared__ float mean[1024];
                    int blockSize = blockDim.x; 
                    unsigned int tid = threadIdx.x;                                                                                                                                        
                    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
                    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

                    mean[tid] = 0.;                                                                                                                                                       

                    while (i < n) {
                        mean[tid] += g_idata[i]; 
                        if (i + blockSize < n){mean[tid] += g_idata[i + blockSize];}
                        i += gridSize;}
                     __syncthreads();                                                                                                                                                       

                    if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
                    if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
                    if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
                    if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
                    if (tid < 32){ warpReduceSum(mean, tid, blockSize);}                                                                                                                                                                      
                    if (tid == 0) {g_mean[blockIdx.x] = mean[0];} 
                   __syncthreads(); 

                }''', 'sumKernel')

    def normaliseUnderMask(self, volume):

        self.fast_sum_stdv *= 0.
        self.fast_sum_mean *= 0.

        self.sumMeanStdv((self.nblocks, 1, 1), (self.num_threads,),
                         (volume, self.mask, self.fast_sum_mean, self.fast_sum_stdv, volume.size),
                         shared_mem=16 * self.num_threads)
        # if self.profile:
        #     t_end = self.stream.record()
        #     t_end.synchronize()
        #
        #     time_took = self.cp.cuda.get_elapsed_time(t_start, t_end)
        #     print(f'nxcc norm s: \t{time_took:.3f}ms')
        #     t_start = self.stream.record()

        # print(meanT[meanT.argmax()-4:meanT.argmax()+2])

        meanT = self.fast_sum_mean.sum() / self.p
        stdT = self.cp.sqrt(self.fast_sum_stdv.sum() / self.p - meanT * meanT)

        # print(meanT, (volume*self.mask).sum()/self.p, stdT, self.cp.sqrt((volume*volume*self.mask).sum()/self.p - meanT**2))

        return ((volume - meanT) / stdT)

    def singleNXCC(self, k, shiftX, shiftY, rot):
        from pytom.agnostic.correlation import nxcc
        import time

        cr_projTrans = self.craw_images[k]
        cr_simProj = self.csim_images[k]

        # Shift raw image and calculated the cross-correlation and score
        optTrans = cr_projTrans.transform(translation=[shiftX, shiftY,0], rotation=[0,0,rot], rotation_order='rzxz')
        nccVal = nxcc(cr_simProj.squeeze(), optTrans.squeeze(), volumeIsNormalized=True)

        return (1 - nccVal)


        optTransNorm = self.normaliseUnderMask(optTrans)
        result = cr_simProj * optTransNorm
        ncc = result.sum() / self.p
        return (1 - ncc.get())

        # self.fast_sum_mean *= 0.
        # self.sumKernel((self.nblocks, 1, 1), (self.num_threads,),
        #                (result, self.fast_sum_mean, result.size),
        #                shared_mem=16 * self.num_threads)
        #
        # ncc = self.fast_sum_mean.sum() / self.p





        # ss = time.time()
        # self.fast_sum_mean *= 0
        # self.fast_sum_stdv *= 0
        # self.sumMeanStdv((self.nblocks, 1, 1), (self.num_threads,),
        #                  (optTrans, self.mask, self.fast_sum_mean, self.fast_sum_stdv, optTrans.size),
        #                  shared_mem=16 * self.num_threads)
        # meanT = self.fast_sum_mean.sum() / self.p
        # stdT = self.cp.sqrt(self.fast_sum_stdv.sum() / self.p - meanT * meanT)
        #
        # self.fast_sum_mean *= 0
        # self.sumMultNorm((self.nblocks, 1, 1), (self.num_threads,),
        #                  (optTrans, cr_simProj, meanT, stdT, self.fast_sum_mean, optTrans.size),
        #                  shared_mem=16 * self.num_threads)
        #
        # nccValExp = self.fast_sum_mean.sum() / self.p
        # t3 = time.time() - ss
        #
        # print(t1, t2, t3)
        # print(nccVal, ncc, nccValExp)





    def calcScore(self, param):
        """
        evaluate the scoring function for given rotation and translation and isotropic isoMagnification
        @param param: rotIP for RotationInPlane, shiftX and shiftY are traslations in X and Y axis, isoMag is the isotropic isoMagnification
        @type param: L{pytom.basic.transformations.general_transform2d}, L{pytom_volume.vol}
        @return: Score
        @rtype: L{float}
        @author: GvdS
        """
        from pytom.gpu.initialize import xp
        score = 0

        for k in range(len(param)//3):

            score += self.score(k, param[k*3], param[k*3+1], param[k*3+2])

        if score < self.low_score:
            self.low_score = score
            self.params = param
            print(f'score: {score:.3f}, {xp.abs(xp.array(param)).sum().get():.1f}')
        return score


    def calcBestFitXY(self, params, show=False):
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.correlation import nxcc
        import pytom.voltools as vt
        from pytom.agnostic.filter import ramp_filter
        cos,sin = xp.cos, xp.sin
        Z1,Y,Z2 = self.angs

        if len(params) == 3:
            Z2 = 0

        center_proj = xp.zeros([2],dtype=xp.float32)
        center_proj[0] += self.src.shape[0] // 2
        center_proj[1] += self.src.shape[1] // 2

        tr11 = cos(Y) * cos(Z1) * cos(Z2) - sin(Z1) * sin(Z2)
        tr31 = -sin(Y) * cos(Z2)
        tr22 = -cos(Y) * sin(Z1) * sin(Z2) + cos(Z1) * cos(Z2)

        tr = xp.array([tr11, tr22, tr31], dtype=xp.float32)

        rec = xp.zeros_like(self.backProjRef,dtype=xp.float32)

        #rec = rec[:,:,:int(50+sin(Y)*rec.shape[2]//2)]


        rec2 = xp.zeros_like(self.backProjRef, dtype=xp.float32)
        cent_rec = xp.array([rec.shape[0]//2, rec.shape[1]//2, rec.shape[2]//2], dtype=xp.int32)
        dims = xp.array(list(rec.shape),dtype=xp.int32)
        outputImage = xp.zeros_like(self.src, dtype=xp.float32)

        vt.transform(self.src.astype(xp.float32), rotation=[0, 0, params[2]], rotation_order='rxyz', output=outputImage,
                     device=device, translation=[params[0], params[1], 0], scale=[params[3],params[3], 1], interpolation='filt_bspline')

        image = outputImage.squeeze()

        #image = bandpass_circle(image, high=lowpassFilter * imdim // 2, sigma=lowpassFilter /5. * imdim)
        image *= self.taper_mask
        image = xp.fft.ifftn(xp.fft.fftn(image) * self.weightSlice ).real


        self.reconstruction_wbp((self.Nblocks, 1, 1,), (self.Nthreads, 1, 1), (image.astype(xp.float32), center_proj, xp.array(list(self.src.shape),dtype=xp.int32), rec, cent_rec, dims, tr, rec.size))

        sx,ex,sy,ey = 10, 60,240,-20

        nccVal = nxcc((self.maskWBP*self.refImage)[sx-1:ex-1,sy-1:ey-1], (self.maskWBP*(self.reconstruction.sum(axis=2)+rec.sum(axis=2)))[sx:ex,sy:ey], volumeIsNormalized=False)

        if show:
            print(params)

            import matplotlib
            try: matplotlib.use('Qt5Agg')
            except:
                pass
            from pylab import imshow, show, subplots
            fig, ax = subplots(2, 1, figsize=(5, 10))
            ax[0].set_title(nccVal)
            #ax[0].imshow( (reconstruction).sum(axis=2)[50:-50,50:-50].get())
            ax[0].imshow((self.reconstruction.sum(axis=2)+rec.sum(axis=2)*1).get()[sx:ex,sy:ey])
            #ax[0].imshow(src.get()[50:-50,50:-50])
            ax[1].imshow(self.refImage.get()[sx-1:ex-1,sy-1:ey-1])
            show()

        # print(nccVal, float(params[0]), float(params[1]))
        return 1-nccVal

    def iterAlign(self,tiltSeriesName='', markerFileName='', lastProj=None, volumeName='', volumeFileType='em',
                  voldims=[464,464,464], recCent=[0,0], tiltSeriesFormat='.mrc', firstProj=0,
                  irefmark=0, ireftilt=0, handflip=False, alignedTiltSeriesName='', weightingType=-1,
                  alignResultsFile='', lowpassFilter=0.9, projBinning=8, verbose=True, volsize=100,
                  numIter=1, startCoor=[0,0,0], CRawDir='CRaw/', CSimDir='CSim/', simulDir='simulDir/', outdir='./',
                  interp='filt_bspline', profile=True, gpuID=-1):
        import matplotlib
        try: matplotlib.use('Qt5Agg')
        except:pass
        import sys, os, re, subprocess
        from scipy import misc
        import numpy as np
        import scipy.optimize
        from pytom_volume import transformSpline, vol
        from pytom.reconstruction.ccAlign import modelProjCrop
        from pytom.basic.files import read_em, write_em, read_em_header, write_em_header
        from pytom.basic.files import read as read_c
        from pytom.agnostic.tools import taper_edges, paste_in_center
        from pytom.gpu.kernels import reconstruction_wbp_text
        from pytom.agnostic.transform import resize
        from pytom.agnostic.io import read, write
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.filter import ramp_filter
        from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
        from pytom.reconstruction.TiltAlignmentStructures import TiltSeries
        from pytom.reconstruction.writeAlignedProjections import writeAlignedProjections
        from pytom.reconstruction.ccAlign import modelProjCrop
        from pytom.gui.guiFunctions import fmtAR, headerAlignmentResults, datatypeAR, loadstar
        import time
        import pytom.voltools as vt


        if gpuID is None: gpuID=-1

        for dirname in (simulDir, CRawDir, CSimDir, outdir):
            if not os.path.exists(dirname): os.mkdir(dirname)

        filenameVolume = volumeName
        #print(f'Generate tomogram: {filenameVolume)')

        ########### 3D reconstruciton with WBP using either ramp or exact filter ##############
        if 0 or not os.path.exists(filenameVolume):
            alignWeightReconstruct(tiltSeriesName=tiltSeriesName, markerFileName=markerFileName, lastProj=lastProj,
                                   volumeName=volumeName, volumeFileType=volumeFileType, voldims=voldims, recCent=recCent,
                                   tiltSeriesFormat=tiltSeriesFormat, firstProj=firstProj, irefmark=irefmark,
                                   ireftilt=ireftilt, handflip=handflip, alignedTiltSeriesName=alignedTiltSeriesName,
                                   weightingType=weightingType, lowpassFilter=lowpassFilter, projBinning=projBinning,
                                   verbose=False, alignResultFile=alignResultsFile, profile=verbose, gpuID=gpuID)

        recVolume = read_c(filenameVolume)

        tiltSeries = self.genTiltSeries(tiltSeriesName=tiltSeriesName, alignedTiltSeriesName=alignedTiltSeriesName,
                                        firstProj=firstProj, lastProj=lastProj, tiltSeriesFormat=tiltSeriesFormat, markerFileName=markerFileName)


        if not (alignResultsFile is None) and os.path.exists(alignResultsFile):
            alignmentResults = loadstar(alignResultsFile, dtype=datatypeAR)
            indices = list(range(len(alignmentResults['TiltAngle'])))

        else:

            indices = list(range(len(tiltSeries._ProjectionList)))
            angs = np.array([projection.getTiltAngle() for projection in tiltSeries._ProjectionList])

            alignmentResults = np.zeros((len(indices)),dtype=datatypeAR)
            alignmentResults['TiltAngle'] = angs
            alignmentResults['Magnification'] = 1
            alignResultsFile = os.path.join(outdir, f'alignmentResultsFile_it.txt')

        outname=alignResultsFile

        startCoor = [ 104,0, 104]

        imdim = max(recVolume.sizeX(), recVolume.sizeY())

        self.low_score = 9999

        self.nblocks = int(xp.ceil(recVolume.sizeX()*recVolume.sizeY() / self.num_threads / 2))
        self.fast_sum_mean = xp.zeros((self.nblocks ), dtype=xp.float32)
        self.fast_sum_stdv = xp.zeros((self.nblocks ), dtype=xp.float32)
        self.mask = xp.ones((recVolume.sizeX(),recVolume.sizeY(),1),dtype=xp.float32)
        self.p = self.mask.size

        reconstruction = xp.zeros(voldims, dtype=xp.float32)
        theta_angles = xp.deg2rad(xp.array(alignmentResults['TiltAngle']))
        ntilts = len(alignmentResults['TiltAngle'])

        sx, sy = resize(read(alignmentResults['FileName'][0]).squeeze(),1/projBinning).shape
        projections = xp.zeros((sx, sx, ntilts))


        cfreq = abs(alignmentResults['TiltAngle'][1:]-alignmentResults['TiltAngle'][:-1])
        cfreq = float(1/xp.sin(cfreq.min()*xp.pi/180))//1
        self.weightSlice = xp.fft.fftshift(ramp_filter(sx, sx, cfreq, 2))


        x,y = xp.meshgrid(xp.arange(voldims[0]), xp.arange(voldims[2]))
        x -= voldims[0]//2
        y -= voldims[1]//2

        r = xp.sqrt(x**2+y**2)


        self.maskWBP = xp.ones_like(reconstruction[:,:,0])
        #self.maskWBP[r <= voldims[0]//2-2] = 1

        for i in range(ntilts):
            image = resize(read(alignmentResults['FileName'][i]),1/projBinning).squeeze()
            immean = image.mean()
            image = (image - immean) / immean
            image, self.taper_mask = taper_edges(image, sx // 30)
            out = xp.zeros((sx,sx),dtype=xp.float32)
            projections[:, :, i] = paste_in_center(image, out)
            image, self.taper_mask = taper_edges(out, sx // 30)
        dims = xp.array(reconstruction.shape, dtype=xp.int32)

        refid = np.abs(alignmentResults['TiltAngle']).argmin()
        self.reconstruction_wbp = xp.RawKernel(reconstruction_wbp_text, 'reconstruction_wbp')
        Z1, Z2 = 0, 0
        Y = theta_angles[refid]

        cos = xp.cos
        sin = xp.sin

        # preparation


        center_recon = xp.zeros((3), dtype=xp.int32)
        center_recon[0] = int(dims[0] // 2 )
        center_recon[1] = int(dims[1] // 2 )
        center_recon[2] = int(dims[2] // 2 )

        dims_proj = xp.array(projections.shape, dtype=xp.int32)

        self.Nthreads = 1024
        self.Nblocks = int(xp.ceil(reconstruction.size / self.Nthreads).get())


        vol_offsetProjections = xp.zeros((projections.shape[2],2),dtype=xp.float32)
        center_proj = xp.zeros_like(vol_offsetProjections)
        center_proj[:, 0] += dims_proj[0] // 2 + vol_offsetProjections[:, 0]
        center_proj[:, 1] += dims_proj[1] // 2 + vol_offsetProjections[:, 1]

        tr11 = cos(Y)*cos(Z1)*cos(Z2)-sin(Z1)*sin(Z2)
        # tr21 = cos(Y)*sin(Z1)*cos(Z2)+cos(Z1)*sin(Z2)
        tr31 = -sin(Y)*cos(Z2)
        # tr12 = -cos(Y)*cos(Z1)*sin(Z2)-sin(Z1)*cos(Z2)
        tr22 = -cos(Y)*sin(Z1)*sin(Z2)+cos(Z1)*cos(Z2)

        tr = xp.array([tr11, tr22, tr31],dtype=xp.float32)
        # tr32 = sin(Y)*sin(Z2)

        # if not recPosVol is None:
        #     cx, cy, cz = recPosVol[n, :]
        #     print(cx,cy,cz)
        #     center_proj[n,0] += xp.cos(xp.arctan2(cz,cx) + Y)* xp.sqrt(cx**2+cz**2) - cx
        #center_proj = cp.array([dims[0] // 2 + vol_offsetProjections[0, 0, n], dims[1] // 2 + vol_offsetProjections[0, 1, n]], dtype=cp.float32)

        src = xp.array(projections[:,:,refid], dtype=xp.float32)

        outputImage = xp.zeros_like(src, dtype=xp.float32)

        vt.transform(src.astype(xp.float32), rotation=[0, 0, alignmentResults['InPlaneRotation'][refid]], rotation_order='rxyz', output=outputImage, device=device,
                      translation=[alignmentResults['AlignmentTransX'][refid],alignmentResults['AlignmentTransY'][refid], 0],
                     scale=[1/alignmentResults['Magnification'][refid],1/alignmentResults['Magnification'][refid], 1], interpolation='filt_bspline')

        self.outputImage = outputImage.copy()

        self.refImage = self.outputImage[104:-104,:]

        src = xp.fft.ifftn(xp.fft.fftn(outputImage.squeeze()) * self.weightSlice ).real

        # imshow(src.get())
        # show()

        print(self.outputImage.shape, src.shape)
        self.reconstruction_wbp((self.Nblocks,1,1,), (self.Nthreads,1,1), (src.astype(xp.float32), center_proj[refid,:], dims_proj, reconstruction, center_recon, dims, tr, reconstruction.size))

        self.backProjRef = reconstruction.copy()

        # imshow(self.backProjRef.get()[:,:,0])
        # show()

        self.params = np.zeros([ntilts*4],dtype=xp.float32)
        dims = xp.array(list(reconstruction.shape), dtype=xp.int32)
        cent_rec= xp.zeros_like(dims,dtype=xp.int32)
        cent_rec[0] = dims[0]//2
        cent_rec[1] = dims[1] // 2
        cent_rec[2] = dims[2] // 2

        for t in range(numIter):
            if 0:
                # Create simulated projection images from reconstruction
                modelProjCrop(tiltSeries, recVolume, tomogramSizeX, tomogramSizeY, tomogramSizeZ, startCoor,vol_size=voldims,
                              verbose=False, binning=projBinning, alignmentResultsFile=outname, weighting=weightingType,
                              profile=verbose, CRawDir=CRawDir, CSimDir=CSimDir, simulDir=simulDir)

                # Set init values of parmaters to current alignment values
                param = np.zeros((len(alignmentResults['TiltAngle'])*3), dtype=np.float32)

                pp = param.reshape(len(alignmentResults['TiltAngle']),3)
                pp[:,2] = -264
                param = pp.flatten()
                self.params = param

                # Update threshold
                ### if t < 5: eps = 1e01 (works good)
                if t < 1:
                   eps = 1
                elif t < 3:
                   eps = 1
                else:
                   eps = 1

                print(f"eps in iteration {t}: {eps}")

                # Reload images

                from scipy.ndimage.filters import gaussian_laplace

                # self.craw_images = [StaticVolume(gaussian_laplace(read(f'{CRawDir}CRaw_{k:02d}.mrc',keepnumpy=True).squeeze(),3), interpolation=interp, device=device) for k in indices]
                # self.csim_images = [mean0std1(xp.array(gaussian_laplace(read(f'{CSimDir}CSim_{k:02d}.mrc', keepnumpy=True).squeeze(),3)),True) for k in indices]

                self.csim_images = []
                self.craw_images = []
                for k in indices:
                    a = xp.array((read(f'{CRawDir}CRaw_{k:02d}.mrc', keepnumpy=True).squeeze()))
                    b = xp.array((read(f'{CSimDir}CSim_{k:02d}.mrc', keepnumpy=True).squeeze()))
                    # a[a<a.mean()+3*a.std()] = 0
                    # b[b<b.mean()+3*b.std()] = 0

                    self.craw_images.append(a)
                    self.csim_images.append(b)

                sx,sy = self.craw_images[0].shape

                # self.low_score = self.calcScore(param)
                #
                # if verbose:
                #     print( f'start_score it{t:03d}: {self.low_score:.3f}')

                # Optimize alignment parameters
                #optParam = scipy.optimize.fmin_bfgs(self.calcScore, param, maxiter=30, epsilon=1e-01)

                if verbose:
                    start_time = time.time()
                #optParam = scipy.optimize.fmin_cg(self.calcScore, param, maxiter=5, epsilon=eps)

                from pytom.agnostic.correlation import FLCF
                for i in range(len(alignmentResults['TiltAngle'])):
                    res = FLCF(self.craw_images[i], self.csim_images[i])
                    px, py = xp.unravel_index(res.argmax(),res.shape)
                    self.params[i*3:i*3+2] = [px-sx//2,py-sy//2]
                    print(i, res.max(), px-sx//2, py-sy//2)
                    import matplotlib
                    try:
                        matplotlib.use('Qt5Agg')
                    except:
                        pass
                    from pylab import imshow, show, subplots
                    fig,ax = subplots(1,3,figsize=(15,5))
                    ax[0].imshow(self.craw_images[i].get())
                    ax[1].imshow(self.csim_images[i].get())
                    ax[2].imshow(res.get())
                    show()

                # shifts = self.params.reshape((len(self.csim_images),3))
                # shifts[:26, :2] -= shifts[26, :2]
                # shifts[26:, :2] = shifts[26, :2] - shifts[26:, :2]
                # self.params = shifts.flatten()

                # if verbose:
                #     print(f'optimization of alignment parameters took {(time.time()-start_time):.1f} sec')
                #     self.final_score = self.calcScore(self.params)
                #     print(f'final_score it{t:03d}: {self.final_score:.3f}')

            for i in range(0,ntilts):
                cfreq = abs(alignmentResults['TiltAngle'][refid] - alignmentResults['TiltAngle'][i])
                self.cfreq = float(1 / xp.sin(cfreq.min() * xp.pi / 180)) // 1
                self.weightSlice = xp.fft.fftshift(ramp_filter(imdim, imdim, self.cfreq, 2))
                self.reconstruction = xp.zeros_like(self.backProjRef, dtype=xp.float32)

                src = xp.fft.ifftn(xp.fft.fftn(self.outputImage.squeeze()) * self.weightSlice).real
                self.reconstruction_wbp((self.Nblocks, 1, 1,), (self.Nthreads, 1, 1), (
                    src.astype(xp.float32), center_proj, xp.array(list(src.shape), dtype=xp.int32), self.reconstruction,
                    cent_rec, dims, xp.array([1, 1, 0], dtype=xp.float32),
                    self.reconstruction.size))


                scoreMin = 1000
                self.src = xp.array(projections[:,:,i], dtype=xp.float32)
                self.angs = [0, theta_angles[i], 0]
                params = [alignmentResults['AlignmentTransX'][i]/projBinning,
                          alignmentResults['AlignmentTransY'][i]/projBinning,
                          alignmentResults['InPlaneRotation'][i],
                          1/alignmentResults['Magnification'][i]]

                print(params)
                if i == refid:
                    optParam = params


                else:

                    x,y = np.meshgrid(np.arange(-12,12.05,12), np.arange(-12,12.05,12))
                    ks = x.flatten()
                    js = y.flatten()
                    for jj in range(len(js)):
                        j,k = js[jj], ks[jj]
                        score = self.calcBestFitXY([params[0]+k, params[1]+j, params[2], params[3]], 0)
                        if score < scoreMin:
                            scoreMin = score
                            #print(score, params[0]+k, params[1]+j)
                            optParam = [k+params[0],j+params[1],params[2],0]


                # optParam = scipy.optimize.fmin(self.calcBestFitXY, params, maxiter=100)
                # params = [optParam[0], optParam[1]]
                # optParam = scipy.optimize.fmin(self.calcBestFitXY, params, maxiter=100)
                # params = [optParam[0], optParam[1], 90]
                # optParam = scipy.optimize.fmin(self.calcBestFitXY, params, maxiter=100)
                # optParam = [optParam[0], optParam[1], optParam[2], 1]
                # print(i, optParam)
                # print(self.calcBestFitXY( params, 0))
                # optParam = scipy.optimize.fmin_powell(self.calcBestFitXY, params, maxiter=100)
                    print(i, self.calcBestFitXY( optParam, 1), optParam)
                self.params[i*4:i*4+4] = np.array(optParam, dtype=xp.float32)

            print(params[0::4]*projBinning)

            # Update alignment parameters and save results
            alignmentResults['AlignmentTransX'] = (self.params[0::4]+1)*projBinning
            alignmentResults['AlignmentTransY'] = (self.params[1::4]+1)*projBinning
            alignmentResults['InPlaneRotation'] = (self.params[2::4])
            # alignmentResults['Magnification']   = self.params[3::4]

            outname = os.path.join(outdir, f'alignmentResultsFile_it{t:03d}.txt')
            np.savetxt(outname, alignmentResults, fmt=fmtAR, header=headerAlignmentResults)
            tiltSeries.updateAlignmentParams(outname)


            # Update reconstruction
            alignWeightReconstruct(tiltSeriesName=tiltSeriesName, markerFileName='', lastProj=lastProj,
                                   volumeName=f'{volumeName}_it{t:03d}.em', volumeFileType='em', alignResultFile=outname,
                                   voldims=voldims, firstProj=firstProj, weightingType=weightingType,
                                   projBinning=projBinning, profile=verbose, gpuID=gpuID)

            recVolume, header = read_em(f'{volumeName}_it{t:03d}.em')

    def genTiltSeries(self, tiltSeriesName, alignedTiltSeriesName, firstProj, lastProj, tiltSeriesFormat, markerFileName):
        from pytom.reconstruction.TiltAlignmentStructures import TiltSeries, TiltAlignmentParameters

        tiltParas = TiltAlignmentParameters(dmag=True, drot=True, dbeam=False, finealig=True,
                                            finealigfile='xxx.txt', grad=False,
                                            irefmark=0, ireftilt=0, r=None, cent=[2049, 2049],
                                            handflip=False, optimizer='leastsq', maxIter=1000)

        tiltSeries = TiltSeries(tiltSeriesName=tiltSeriesName, TiltAlignmentParas=tiltParas,
                                alignedTiltSeriesName=alignedTiltSeriesName,
                                markerFileName=markerFileName, firstProj=firstProj, lastProj=lastProj,
                                tiltSeriesFormat=tiltSeriesFormat)

        return tiltSeries


if __name__ == '__main__':
    import sys
    #from pytom.reconstruction.TiltAlignmentStructures import TiltAlignmentParameters, TiltSeries, TiltAlignment
    from pytom.tools.script_helper import ScriptHelper, ScriptOption
    from pytom.tools.parse_script_options import parse_script_options
    from pytom.reconstruction.reconstructionFunctions import alignWeightReconstruct
    import numpy
    import os

    options=[ScriptOption(['--tiltSeriesName'], 'Name tilt series - either prefix of sequential tilt series files \
             expected as "tiltSeriesName_index.em/mrc" or full name of stack "tiltSeriesName.st"',
                          arg=True, optional=False),
             ScriptOption(['--tiltSeriesFormat'], 'Format of tilt series (series of "em" or "mrc" images or "st" stack).',
                          arg=True, optional=True),
             ScriptOption(['--firstIndex'], 'Index of first projection.', arg=True, optional=True),
             ScriptOption(['--lastIndex'], 'Index of last projection.', arg=True, optional=True),
             ScriptOption(['--projIndices'], 'Use numbering in filename as index', arg=False, optional=True),
             ScriptOption(['--tltFile'], 'tltFile containing tilt angles.', arg=True, optional=True),
             ScriptOption(['--prexgFile'], 'prexgFile containing pre-shifts from IMOD.', arg=True, optional=True),
             ScriptOption(['--preBin'], 'pre-Binning in IMOD prior to marker determination.', arg=True, optional=True),
             ScriptOption(['--referenceIndex'], 'Index of reference projection used for alignment.', arg=True,
                          optional=True),
             ScriptOption(['--markerFile'], 'Name of EM markerfile or IMOD wimp File containing marker coordinates.',
                          arg=True, optional=True),
             ScriptOption(['--alignmentResultsFile'], 'Name of alignmentResults file markerfile containing per image shift, '
                                                      'tilt angle, in-plane rotation, magnification and filename.',
                          arg=True, optional=True),
             ScriptOption(['--referenceMarkerIndex'], 'Index of reference marker to set up coordinate system.',
                          arg=True, optional=True),
             ScriptOption(['--expectedRotationAngle'], 'Is your tilt series outside of 0-180deg (Specify if yes).',
                          arg=True, optional=True),
             ScriptOption(['--projectionTargets'],
                          'Relative or absolute path to the aligned projections that will be generated + file prefix.\
                          default: "align/myTilt"', arg=True, optional=True),
             ScriptOption(['--fineAlignFile'],
                          'Relative or absolute path to the file with fineAlign parameters (type should be *.dat).',
                          arg=True, optional=True),
             ScriptOption(['--projectionBinning'], 'Binning of projections during read - default: 1.', arg=True,
                          optional=True),
             ScriptOption(['--lowpassFilter'], 'Lowpass filter in Nyquist after binning.', arg=True, optional=True),
             ScriptOption(['--tomogramFile'],
                          'Relative or absolute path to final tomogram (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--fileType'], 'File type (can be em or mrc - no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeX'], 'Size of tomogram in x (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeY'], 'Size of tomogram in y (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--tomogramSizeZ'], 'Size of tomogram in z (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterX'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterY'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--reconstructionCenterZ'],
                          'Center where tomogram will be reconstructed (no tomogram written if not specified).',
                          arg=True, optional=True),
             ScriptOption(['--weightingType'], 'Type of weighting (-1 default r-weighting, 0 no weighting)', arg=True,
                          optional=True),
             ScriptOption(['--noOutputImages'], 'When specified, not output images are saved.', arg=False, optional=True),
             ScriptOption(['--verbose'], 'Enable verbose mode', arg=False, optional=True),
             ScriptOption(['-g', '--gpuID'], 'provide a gpu for running this algorithm', arg=True, optional=False),
             ScriptOption(['-h', '--help'], 'Help.', False, True)]

    helper = ScriptHelper(sys.argv[0].split('/')[-1],
                          description='Align and weight projections, save them and reconstruct tomogram (optional). \n\
                                      See http://pytom.org/doc/pytom/reconstructTomograms.html for documentation.',
                          authors='Friedrich Foerster',
                          options = options)

    if len(sys.argv) == 1:
        print(helper)
        sys.exit()
    try:
        tiltSeriesName, tiltSeriesFormat, firstProj, lastProj, projIndices,\
        tltFile, prexgFile, preBin, referenceIndex, markerFileName, alignResultsFile, referenceMarkerIndex, expectedRotationAngle, \
        projectionTargets, fineAlignFile, projBinning, lowpassFilter, \
        volumeName, filetype, \
        tomogramSizeX, tomogramSizeY, tomogramSizeZ, \
        reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ, \
        weightingType, noOutputImages, verbose, gpuID, help = parse_script_options(sys.argv[1:], helper)
    except Exception as e:
        print(sys.version_info)
        print(e)
        print()
        print(helper)
        sys.exit()

    if help is True:
        print(helper)
        sys.exit()

    if (markerFileName is None) and (alignResultFile is None):
        raise Exception('Please provide either a markerfile or an alignmentResults.txt file')

    gpuID = int(gpuID)

    irefmark, ireftilt = 1, 1

    if not markerFileName is None and os.path.exists(markerFileName):

        # input parameters
        #tiltSeriesName = tiltSeriesPath + tiltSeriesPrefix  # "../projections/tomo01_sorted" # ending is supposed to be tiltSeriesName_index.em (or mrc)
        if not tiltSeriesFormat:
            tiltSeriesFormat = 'em'
        if firstProj:
            firstProj = int(firstProj)  # 1 # index of first projection
        else:
            firstProj = 1

        lastProj = int(lastProj)  # 41 # index of last projection
        ireftilt = int(referenceIndex)  # 21 # reference projection (used for alignment)
        if referenceMarkerIndex:
            irefmark = int(referenceMarkerIndex)  # reference marker (defines 3D coordinate system)
        else:
            irefmark = 1

        try:
            expectedRotationAngle = int(expectedRotationAngle)
        except:
            expectedRotationAngle = 0
        #handflip = handflip is not None  # False # is your tilt axis outside 0-180 deg?
        # output parameters
    if projectionTargets:
        # weighted and aligned projections are stored as alignedTiltSeriesName_index.em
        alignedTiltSeriesName = projectionTargets
        if not os.path.exists(os.path.dirname(projectionTargets)):
            os.mkdir(os.path.dirname(projectionTargets))
    else:
        alignedTiltSeriesName = 'align/myTilt'

    if projBinning:
        projBinning = int(projBinning)  # binning factor
    else:
        projBinning = 1
    if lowpassFilter:
        lowpassFilter = float(lowpassFilter)  # lowpass filter in Nyquist (post-binning)
    else:
        lowpassFilter = 1.

    if weightingType is None:
        weightingType = -1
    else:
        weightingType = int(weightingType)

    # only write projections and do NOT reconstruct tomogram (following parameters would be obsolete)
    onlyWeightedProjections = False
    if not volumeName:
        onlyWeightedProjections = True

    if filetype is None:
        if volumeName:
            filetype = volumeName.split('.')[-1]
        else:
            filetype = tiltSeriesFormat

    if tomogramSizeX is not None or tomogramSizeY is not None or tomogramSizeZ is not None:
        # dimensions of reconstructed tomogram
        voldims = [int(tomogramSizeX), int(tomogramSizeY), int(tomogramSizeZ)]
    else:
        onlyWeightedProjections = True
        voldims = [0, 0, 0]       # offset from center of volume - for example choose z!=0 to shift in z (post-binning coordinates)
    if reconstructionCenterX:
        reconstructionCenterX = int(reconstructionCenterX)
    else:
        reconstructionCenterX = 0
    if reconstructionCenterY:
        reconstructionCenterY = int(reconstructionCenterY)
    else:
        reconstructionCenterY = 0
    if reconstructionCenterZ:
        reconstructionCenterZ = int(reconstructionCenterZ)
    else:
        reconstructionCenterZ = 0
    reconstructionPosition = [reconstructionCenterX, reconstructionCenterY, reconstructionCenterZ]
    if preBin:
        preBin=int(preBin)

    if noOutputImages is None:
        write_images = True
    else:
        write_images = False

    outMarkerFileName = 'MyMarkerFile.em'

    alignResultsFile = '' if alignResultsFile is None else alignResultsFile

    if not alignResultsFile and 'alignmentResults.txt' in os.listdir(os.path.dirname(tiltSeriesName)):
        alignResultsFile = os.path.join(os.path.dirname(tiltSeriesName), 'alignmentResults.txt')

    outfile = ''
    outfolder = os.path.dirname(projectionTargets)

    if 'reconstruction/WBP' in projectionTargets or 'reconstruction/INFR' in projectionTargets:
        outfolder = os.path.dirname(outfolder)
    reconstructionAlgorithm = os.path.basename(outfolder)
    tomogramID = os.path.basename(os.getcwd())
    outfile = os.path.join(outfolder, 'markerLocations_{}_irefmark_{}.txt'.format(tomogramID, referenceMarkerIndex))

    if verbose:
        print("Tilt Series: "+str(tiltSeriesName)+", "+str(firstProj)+"-"+str(lastProj))
        print("Index of Reference Projection: "+str(referenceIndex))
        print("Marker Filename: "+str(markerFileName))
        print("TltFile: "+str(tltFile))
        print("prexgFile: "+str(prexgFile))
        print("Index of Reference Marker: "+str(referenceMarkerIndex))
        print("Expected Rotation Angle: "+str(expectedRotationAngle))
        print("Projection Targets: "+str(projectionTargets))
        print("FineAlignmentFile: "+str(fineAlignFile))
        print("Binning Factor of Projections: "+str(projBinning)+", lowpass filter (in Ny): "+str(lowpassFilter))
        print("Name of Reconstruction Volume: "+str(volumeName)+" of Filetype: "+str(filetype))
        print("Reconstruction size: "+str(voldims))
        print("Reconstruction center: "+str(reconstructionPosition))
        print("write only aligned projections out: "+str(onlyWeightedProjections))
        print(f"Marker locations are written to: {outfile}")


    fa = FiducialLessAlignment()

    fa.iterAlign(tiltSeriesName=tiltSeriesName, markerFileName=markerFileName, lastProj=lastProj, volumeName=volumeName,
                 volumeFileType=filetype, voldims=voldims, recCent=reconstructionPosition,
                 tiltSeriesFormat=tiltSeriesFormat, firstProj=firstProj, irefmark=irefmark, ireftilt=ireftilt,
                 handflip=float(expectedRotationAngle)*numpy.pi/180, alignedTiltSeriesName=alignedTiltSeriesName,
                 weightingType=weightingType, alignResultsFile=alignResultsFile,lowpassFilter=lowpassFilter,
                 projBinning=projBinning, gpuID=-1)



