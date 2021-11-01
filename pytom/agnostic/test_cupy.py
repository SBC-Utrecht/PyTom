



def plot_central_sections(vol):
    import matplotlib
    try:
        matplotlib.use('Qt5Agg')
    except:
        pass
    from pylab import subplots, show

    fig,ax = subplots(1,3,figsize=(15,5))
    dx,dy,dz = vol.shape
    ax[0].imshow(vol[dx//2,:,:])
    ax[1].imshow(vol[:,dy//2,:])
    ax[2].imshow(vol.sum(axis=2))
    show()

if __name__=='__main__': 
    import sys
    from scipy.ndimage import rotate as ROTATE
    import cupy as xp
    import numpy as np
    from pytom.agnostic.correlation import *
    from pytom.agnostic.tools import paste_in_center
    from cupyx.scipy.ndimage import rotate
    num_angles, size = map(int, sys.argv[1:3])
    
    template = xp.zeros((size, size, size), dtype=xp.float32)

    temp = xp.zeros((64,64,64))
    temp[16:-16,16:-16,16:-16] = 1.
    
    for i in range(3):
        tr = ROTATE(temp,0*i,axes=(2,0),reshape=False,order=5)
        ii = i+1
        z,y,x = 150, 150, ii*75
        template[z-32:z+32, y-32:y+32, x-32:x+32] = temp

    volume = xp.random.rand(size, size, size)
    volume = volume.astype(xp.float32)
  
    import time

    #s = time.time()

    import cupy as xp

    gpu=True

    if gpu:
        vcp = xp.array(template+volume)
        tcp = xp.array(temp)
    else:
        vcp = template
        tcp = temp
    #fvol = xp.fft.fftn(vcp)
    #temp = xp.conj(xp.fft.fftn(tcp))
 
    #m = xp.abs(xp.fft.ifftn(fvol * temp)).get()
  
    #print(time.time()-s)
    mask = tcp.copy()

    meanT = meanUnderMask(tcp, mask, gpu=gpu)
    stdT = stdUnderMask(tcp, mask, meanT, gpu=gpu)

    if stdT > 1E-09:
        temp2 = (tcp - meanT) / stdT
        temp2 = temp * mask
    else:
        temp2 = tcp
        mask2 = mask
    if vcp.shape[0] != temp2.shape[0] or vcp.shape[1] != temp2.shape[1] or vcp.shape[2] != temp2.shape[2]:
        tempV = xp.zeros(vcp.shape)
        temp2 = paste_in_center(temp2, tempV,gpu=gpu)

    if vcp.shape[0] != mask.shape[0] or vcp.shape[1] != mask.shape[1] or vcp.shape[2] != mask.shape[2]:
        maskV = xp.zeros(vcp.shape)
        mask2 = paste_in_center(mask, maskV, gpu=gpu)

    meanV = meanVolUnderMask(vcp, temp2, gpu=gpu)
    stdV = stdVolUnderMask(vcp, mask2, meanV, gpu=gpu)




    from pytom.agnostic.transform import rotate3d

    s = time.time()
    for i in range(num_angles):
        tcp2 = xp.array(rotate3d(temp, 10, 10, 10))
        m = FLCF(vcp, temp2, mask=mask2, stdV=stdV, gpu=gpu)

    print((time.time()-s))


    

