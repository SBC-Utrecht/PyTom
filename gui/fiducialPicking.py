
import sys
import os
import glob
import numpy
import pytom.gui.mrcOperations
from time import time
#from pytomgui import mrcOperations
from numpy.fft  import fft2, ifft2, fftshift
from numpy import *
from numpy import int32, arange, conj, zeros, ones, bool8, sqrt, newaxis, tan, median
from pytom.gui.mrcOperations import read_mrc as read

from numpy import tan

from scipy.ndimage import gaussian_filter, laplace, label, median_filter, maximum_filter, convolve
from scipy.ndimage import binary_closing, binary_fill_holes, binary_erosion, binary_dilation, center_of_mass
from scipy.ndimage.filters import minimum_filter
from scipy.ndimage.morphology import generate_binary_structure
from scipy.signal import wiener, argrelextrema



from skimage.morphology import remove_small_objects, watershed
from skimage.feature import canny, peak_local_max

def diff(vec0, vec1):
    return numpy.linalg.norm(vec0-vec1)

def closest_point(ref_vec, points,cutoff=3.4, it=0, alpha=270):

    if alpha < 0.: alpha += 360
    if alpha < 45 or alpha > 315 or (alpha > 135 and alpha < 225):
        vertical = True
    else:
        vertical = False

    t90 = tan((alpha - 90) * numpy.pi / 180.)
    t00 = tan((alpha) * numpy.pi / 180.)

    c = numpy.zeros(len(points),dtype=float)

    if not vertical:
        for n, p in enumerate(points):
            c[n] = diff(ref_vec,p)
            if abs(ref_vec[1] - (p[1] - t90 * (p[0]-ref_vec[0]))) > cutoff:
                c[n] = 1000
    else:
        for n, p in enumerate(points):
            if vertical and abs(ref_vec[0] - (p[0] - t00 * (p[1]-ref_vec[1]))) > cutoff:
                c[n] = 1000

    try:
        for i in range(it): c[c.argmin()] = 10000
        return c.argmin(), c.min()
    except:
        return -1


def dist_z(y0,y1,a1):
    z1 = (y1-y0*numpy.cos(a1) )/ (numpy.sin(a1))
    return z1
    
def compare(cur, ref, imnr, cutoff=3.4,v=False,tiltangles=[],y0=0,a0=0,pos=1,cutoff_dist=10,tiltaxis=0):
    angle = tiltaxis
    dist_map = [1000]*(ref.shape[0])
    index_map = numpy.zeros((ref.shape[0]),dtype=int)
    coordinates = numpy.zeros_like(ref)
    for i,p in enumerate(cur):
        if p.sum() == 0: continue
        index,dist = closest_point(p,ref,cutoff,alpha=angle)
        
        if dist< 4*cutoff_dist:
            index_map[i] = index+1

            tgt = 4
            #if index == tgt :
            #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
            #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]
                #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(numpy.round(p[0])), int(numpy.round(p[1])), dist_z(ref[index][0], p[0], d_theta*numpy.pi/180.) ), p[0]
                #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+2, imnr, int(numpy.round(d_theta)), int(numpy.round(tiltangles[imnr])), dist_z(ref[tgt+1][0], p[0], d_theta*numpy.pi/180.) ), p[0]
            
            #if index == tgt+1 :
            #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
            #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]   
                #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+1, imnr, int(numpy.round(p[0])), int(numpy.round(p[1])), dist_z(ref[tgt][0], p[0], d_theta*numpy.pi/180.) ), p[0]
                #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(numpy.round(p[0])), int(numpy.round(p[1])), dist_z(ref[index][0], p[0], d_theta*numpy.pi/180.) ), p[0]
            
    count = 0
    #print
    
    while count < 10 and len(index_map) != len(numpy.unique(index_map))+(1*(index_map==0)).sum()-1:
        #if 1: print imnr, index_map
        
        for ind in numpy.unique(index_map):
            if ind == 0: continue
            
            closest, second = [], []
            for n,index in enumerate(index_map):
                if index==0 or index != ind or cur[n].sum() == 0: continue
                for ii in range(len(index_map)):
                    ind1,dis1 = closest_point(cur[n],    ref,it=ii, alpha=angle)
                    if ind1+1==ind: break
                ind2,dis2 = closest_point(cur[n],    ref,it=ii+1, alpha=angle)
                closest.append([n,ind1+1,dis1,ind2+1,dis2])
            
            candidates = []
            m = [0,0,9999,0,9999]
            for c in closest:
                if c[2] < cutoff_dist:
                    candidates.append(c)
                if c[2] < m[2]:
                    m = c
            
            #if v and len(candidates) > 1: print candidates
            if len(candidates) > 2:
                
                for c in closest:
                    if m[0]==c[0]: continue

                    if c[4] < 100:
                        if v: print ('changed', index_map[c[0]], index_map[c[1]])
                        index_map[c[0]]==c[3]
                    else: index_map[c[0]] = 0
            
            elif len(candidates) ==2:
                todo=True
                c1,c2 =candidates
                dist,dis,d1,d2 = c1[2],c1[4],c2[2],c2[4]

                std2,std1 = numpy.std( [d1,dis] ), numpy.std( [d2,dist] )

                #if imnr == 15:
                #    print c1
                #    print c2

                if std2<4 and std2 <4:
                    p1,p2 = cur[c1[0]],cur[c2[0]]
                    
                    diff = abs(p1[1]*4.-p2[1]*4.)
                    #print diff
                    if diff>3.:
                        #print p1, p2
                        if p1[1] > p2[1]:
                            m = c1
                            index_map[c2[0]] = c2[3]
                        else:
                            index_map[c1[0]] = c1[3]
                            m = c2
                        #print imnr, index_map
                    else:
                        if std1 < std2:                                                                                                                       
                            m = c1                                                                                                                              
                            index_map[c2[0]] = c2[3]                                                                                                           
                        else:                                     
                            index_map[c1[0]] = c1[3]
                            m = c2
                        
                elif numpy.std( [d1,dis] ) < numpy.std( [d2,dist] ) and d2 < cutoff_dist and dis < cutoff_dist:
                    index_map[c1[0]] = c1[3]
                    m = c2
                    if v: print (count,'option a')
                        
                elif numpy.std( [d1,dis] ) < numpy.std( [d2,dist] ) and d2 < cutoff_dist and dis < cutoff_dist:
                    index_map[c2[0]] = c2[3]
                    m = c1
                    if v: print (count, 'option b')
                    todo=False
                elif d1> dist and dist < 100:
                    m = c1
                    if d1< 100:
                        index_map[c2[0]] = c2[3]
                    else: index_map[c2[0]] = 0
                    if v: print (count, 'option c')
                    
                elif d1 <= dist and d1 < 100:
                    m = c2
                    if dist< 10 and dis <20:
                        index_map[c1[0]] = c1[3]
                    else: index_map[c1[0]] = 0
                    if v: print (count, 'option d', c1, c2, todo)
                        
                else:
                    if v: print (count, 'no_change' )

                    
            if len(closest) > 1:
                for c in closest:
                    if m[0]==c[0]:continue
                    if c[4] < cutoff_dist: index_map[c[0]]==c[3]
                    else: index_map[c[0]] = 0

        count +=1
        if v: print (imnr, index_map)
    
    for n,p in enumerate(cur):
        if index_map[n]: coordinates[index_map[n]-1] = p

    return index_map, coordinates

def sort( obj, nrcol ):
    obj.sort(key=lambda i: float(i[nrcol]))


def detect_shifts_many(mark_frames, max_shiftx=5.5,max_shifty=2.,sd=1,diag=False, image=[]):
    from time import time


    if diag: s=time()

    tilt_nr,mark_nr, k = mark_frames.shape
    dimx,dimy = image.shape

    out = numpy.zeros((dimx,dimy,tilt_nr))

    frame_shifts = numpy.zeros((tilt_nr,mark_nr,2),dtype=float)
    fs = numpy.zeros((tilt_nr,2),dtype=float)
    outline_detect_shifts = ''
    num_shifts = 0

    for itilt in range(tilt_nr):
        for imark in range(mark_nr):
            x,y = mark_frames[itilt][imark]
            #if x> dimx-1 or y> dimy-1: print x,y
            if x > 0 and y > 0 and x < dimx-1 and y< dimy-1: out[int(round(x))][int(round(y))][itilt] = 1
    for itilt in range(tilt_nr-1):
        c, shifty, shiftx, m = detect_shift(gaussian_filter(out[:,:,itilt],sd), gaussian_filter(out[:,:,itilt+1],sd),image=image)
        
        if m > 0.001:#abs(shiftx) > max_shiftx or abs(shifty) > max_shifty: 
            print ("\t{:2d}-{:2d} {:4.0f} {:4.0f} {:4.3f}".format(itilt,itilt+1,shiftx,shifty,m))
            frame_shifts[itilt+1:] += [ shiftx, shifty]
            fs[itilt+1:] += [ shiftx, shifty]
            num_shifts +=1
            outline_detect_shifts += "\t{:2d}-{:2d} {:4.0f} {:4.0f} {:4.3f}\n".format(itilt,itilt+1,shiftx,shifty,m)
    
    frame_shifts *= (mark_frames.sum(axis=-1) > 0)[:,:,numpy.newaxis]
    if diag: print ('detect frame shifts: {:.0f} msec'.format(1000*(time()-s)))
    print
    return frame_shifts, num_shifts, outline_detect_shifts, fs


def detect_shift(arr0,arr1,image=[]):
    x,y = image.shape
    cross = abs(fftshift( ifft2(fftshift(fftshift(fft2(arr0))*conj(fftshift(fft2(arr1)))))))**2
    locx,locy =  (abs((cross))**2).flatten().argmax()%y, (abs((cross))**2).flatten().argmax()/y
    return cross, locx-y/2, locy-x/2, cross[int(locy)][int(locx)]



def detect_shifts_few(mark_frames, maxshift=9, max_shiftx=5.5, max_shifty=2.,diag=False, image=[]):
    from time import time
    if diag: s = time()

    x,y,z = mark_frames.shape
           
    distance_matrix = 1000*numpy.ones((x-1,y),dtype=float)
    shift = numpy.zeros_like(mark_frames)
    size_shift = numpy.ones((x-1,y),dtype=float)

    for ii in range(x-1):
        distance_matrix[ii,:], shift[ii,:,:], size_shift[ii,:] = pair(mark_frames[ii],mark_frames[ii+1],pr=(ii==8))
        #print ii, distance_matrix[ii]
    frame_shifts = numpy.zeros((x,y,2),dtype=float)
    fs = numpy.zeros((x,2),dtype=float)
    num_shifts = 0

    outline_detect_shifts = ''
    for iii in range(x-1):     
        sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]<150)
        d_shift = size_shift[iii][sel]
        shifty, shiftx = numpy.median(shift[iii,:,1][sel]), numpy.median(shift[iii,:,0][sel])
        #if iii == 31: print sorted(d_shift)
        #print iii, numpy.std(d_shift), d_shift

        n, h = numpy.histogram(d_shift, bins=10, range=None, normed=False, weights=None, density=None)
        #if iii == 31:
        #    print n, h
            #print numpy.std(h[n.argmax():n.argmax()+1]
        #distance_matrix[iii], d_shift, shift[iii]                                                                   
        if abs(shifty) > max_shifty or abs(shiftx) > max_shiftx :
            if numpy.std(d_shift) > 5:
                sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]<85)
                d_shift = size_shift[iii][sel]
                shifty, shiftx = numpy.median(shift[iii,:,1][sel]), numpy.median(shift[iii,:,0][sel])
                #print iii, numpy.std(d_shift), d_shift
                #print

            if numpy.std(d_shift) > 5:
                #diff = h[-1]-h[0]
                #total = len(d_shift)
                sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]< h[min(len(n),n.argmax()+1)] )*(size_shift[iii]>h[n.argmax()])
                d_shift = size_shift[iii][sel]
                shifty, shiftx = numpy.median(shift[iii,:,1][sel]), numpy.median(shift[iii,:,0][sel])
                
                
            if numpy.std(d_shift) > 5:
                continue
            outline_detect_shifts += "\t{:2d}-{:2d} {:4.1f} {:4.1f} {:4.1f} | {:4.0f} {:4.0f}\n".format(iii,iii+1, numpy.median(d_shift ), numpy.std(d_shift), numpy.mean(d_shift), shiftx,shifty)
            frame_shifts[iii+1:] += [ shiftx, shifty]
            fs[iii+1:] += [ shiftx, shifty]
            num_shifts +=1

    frame_shifts *= (mark_frames.sum(axis=-1) > 0)[:,:,numpy.newaxis]
    print ("encountered {} frame shifts.".format(num_shifts))
    print (outline_detect_shifts   )
 
    if diag: print ('detect frame shifts: {:.0f} msec'.format(1000*(time()-s)))

    return frame_shifts, num_shifts, outline_detect_shifts,fs

def pair(l1,l2,pr=0):
    len_l1,len_l2 = len(l1), len(l2)
    com1 = numpy.zeros((len_l1,len_l1,2))
    com2 = numpy.zeros((len_l2,len_l2,2))
    
    for i in range(len_l1):
        com1[i,:,:] = l1-l1[i]
    for j in range(len_l2):
        com2[j,:,:] = l2-l2[j]

    
    
    
    dl = numpy.ones((len_l1))*1000
    sl = numpy.ones_like(l1)*1000
    nl = dl.copy()
    
    for n, el in enumerate(com1.reshape(len_l1**2,2)):
        if el.sum() == 0:
            continue
        diff = numpy.linalg.norm(com2-el,axis=-1)

        n2,distance=diff.argmin(),diff.min()

        x1,y1=n/len_l1,n%len_l1
        x2,y2=n2/len_l1,n2%len_l1

        shift = (l1[x1]+l1[y1])/2-(l2[x2]+l2[y2])/2
        norm = numpy.linalg.norm(shift)

            
        if dl[x1] > distance and l1[x1].sum() and l1[y1].sum() and l2[x2].sum() and l2[y2].sum():
            dl[x1] = distance
            sl[x1] = shift
            nl[x1] = norm
        
    return dl,sl,nl


def frames_without_indexed_fiducials(c,tiltangle):
    empty_frames = []
    try:
        a, b = int((c[:,0][c[:,0]>0]).mean()), int((c[:,1][c[:,1]>0]).mean())
    except:
        empty_frames.append(tiltangle)
    return empty_frames


def index_potential_fiducials(frames, mark_frames, frame_shifts, plot=True, cut=3, zero_angle=20, tiltangles=[],
                              add_marker=False, excluded=[], user_coords=[],diag=False,tiltaxis=0):
    from time import time
    s = time()
    c = 300
    r = len(frames)
    cmap = ['r','orange','lightgreen','yellow','b','cyan','magenta','pink','violet','black']


    mark_frames += frame_shifts.astype(float)
    
    coordinates_sorted = numpy.zeros_like(mark_frames)
    index_map = -1*numpy.ones((r,len(mark_frames[0])),dtype=int)
    
            
    cur = 1
    for i in range(len(mark_frames[zero_angle])):
        if mark_frames[zero_angle][i][0] > frame_shifts[zero_angle][i][0]-1.5:
            index_map[zero_angle][i] = cur
            coordinates_sorted[zero_angle][i] = mark_frames[zero_angle][i]
            cur +=1

    for tiltangle in range(zero_angle+1,r):
        cur_fid = mark_frames[tiltangle,:].copy()
        ref_fid = coordinates_sorted[tiltangle-1,:].copy()
        c =ref_fid[:,:]
    
        for temp_i in numpy.arange(tiltangle-1,zero_angle,-1):
            markers = coordinates_sorted[temp_i-1,:].copy()
            absent = (((ref_fid == [0,0]).sum(axis=1) ==2)*1)[:,numpy.newaxis]*markers
            ref_fid += absent
     
        #print cur_fid.T
        if len(user_coords):
            for n,(x,y) in enumerate(user_coords[tiltangle]):
                if x>0 and y>0:
                    print ('found user selected marker: {:6.1f},{:6.1f}'.format(x,y))
                    ref_fid[n] = [x+1,y+1]
             

        index_cur_fid, coordinates_cur_sorted = compare(cur_fid, ref_fid, tiltangle, cutoff=cut,tiltangles=tiltangles,pos=1,tiltaxis=tiltaxis)
        index_map[tiltangle,:] = index_cur_fid
        #print coordinates_cur_sorted
        coordinates_sorted[tiltangle,:] = coordinates_cur_sorted
    #frames_without_indexed_fiducials(coordinates_sorted[tiltangle,:],tiltangle)
    #for j, [mx,my] in enumerate(coordinates_cur_sorted):
    #    if mx+my> 0: self.ax.scatter(tiltangle+(10+len(self.fnames))*j, mx,c=cmap[j%9], cmap=plt.cm.jet)
       
    for tiltangle in numpy.arange(zero_angle-1,-1,-1):
        cur_fid = mark_frames[tiltangle,:].copy()
        ref_fid = coordinates_sorted[tiltangle+1,:].copy()

        for temp_i in numpy.arange(tiltangle+1,zero_angle,+1):
            markers = coordinates_sorted[temp_i+1,:]
            absent = (((ref_fid == [0,0]).sum(axis=1) ==2)*1)[:,numpy.newaxis]*markers
            ref_fid += absent

        index_cur_fid, coordinates_cur_sorted = compare(cur_fid, ref_fid, tiltangle, cutoff=cut, tiltangles=tiltangles, pos=0, tiltaxis=tiltaxis)
        index_map[tiltangle,:] = index_cur_fid
        coordinates_sorted[tiltangle,:] = coordinates_cur_sorted
    #frames_without_indexed_fiducials(coordinates_sorted[tiltangle,:], tiltangle)

    #print coordinates_sorted
    printline = ''
    cntr = 0
    tot = 0
    listdx= []
    #for nr in numpy.unique(index_map)[:]:

    outc = numpy.zeros_like(coordinates_sorted)
    imark = 0
    for i in range(len(frames)):
        if excluded[i]:
            coordinates_sorted[i,:,:] *= 0
    
    for nr in range(coordinates_sorted.shape[1]):
        if nr+1 > 0.5:# and (index_map == nr).sum() > 2:
            cntr+=1
            cur = 0
            for j, [mx,my] in enumerate(coordinates_sorted[:,nr,:]):
                #if mx >0 and my >0: cur += 1
                tempn = 0
                for mcnt, (markx,marky) in enumerate(mark_frames[j,:,:]):
                    
                    if abs(mx-markx) < 0.1 and abs(marky-my) < 0.1:
                        
                        if (markx-frame_shifts[j,mcnt,0]) > 0.01 and (marky-frame_shifts[j,mcnt,1] ) > 0.01:
                            #print (nr,j, marky, marky-frame_shifts[j,mcnt,1])
                            cur += 1
                            coordinates_sorted[j,nr,:] -= frame_shifts[j,mcnt,:]
                            tempn+=1
                
            if cur > 19:
                tot +=cur

                index = coordinates_sorted[:,nr,0].argmax()
                mx,my = coordinates_sorted[index,nr,:]
                printline +=  "\tindex {:4d}: {}/{}\n".format(int(cntr), cur, r)
                listdx.append('{:02d}/{}'.format(cur,r-sum(excluded)))
                outc[:,imark,:] = coordinates_sorted[:,nr,:]
                try: add_marker(['Marker_{:03d}'.format(imark), '{:02d}/{:02d}'.format(cur,r)])
                except: pass
                imark += 1
                
            elif cur>5:
                if nr > 0.5 and (index_map == nr).sum() > 2:
                    index = coordinates_sorted[:,nr,0].argmax()
                    mx,my = coordinates_sorted[index,nr,:]
                    if (mx+my) > 0:
                        printline += "\tindex {:4d}: {}/{}\n".format(int(cntr), cur, r)
                        listdx.append('{:02d}/{}'.format(cur,r-sum(excluded)))

                    outc[:,imark,:] = coordinates_sorted[:,nr,:]

                    try: add_marker(['Marker_{:03d}'.format(imark), '{:02d}/{:02d}'.format(cur,r)])
                    except: pass
                    imark += 1
            else:
                coordinates_sorted[:,nr,:] = 0.*coordinates_sorted[:,nr,:] -1


    print ("indexed {}/{} potential fiducials.".format(tot, (mark_frames.sum(axis=-1) > 0).sum()) )
    print (printline)


    frame_shifts_sorted = numpy.zeros_like(frame_shifts)
    
    '''
    for j in range(frame_shifts.shape[0]):
        for i in range(frame_shifts.shape[1]):
            if numpy.abs(coordinates_sorted[j][i]).sum() > 0.1:
                frame_shifts_sorted[j][i] = frame_shifts[j][0]
    '''
    if diag:
        print ('index frames: {:.0f} msec'.format(1000*(time()-s)) )
        s =time()
    mark_frames -= frame_shifts
    return outc, index_map, frame_shifts_sorted, listdx
    
def find_potential_fiducials(frames, frames_full, bin, bin_full, cropsize=50,target=0,fid_list=[],proc_id=0,num_procs=1,sen=False):
    from time import time

    s = time()
    
    list_cx_cy_imnr = []
    
    #fiducials = numpy.zeros_like(frames_full[0])
    #cf = numpy.zeros_like(frames_full[0])
    x,y = frames[0].shape
    sum2 = 0
    fact = (bin*1.)/bin_full
    ss = cropsize//2
    
    for imnr in range(len(frames)):
        #if len(fid_list) and not imnr==target: continue
        st = time()
        image =  gaussian_filter(frames[imnr],1)
        #cf = numpy.zeros_like(frames_full[0])
        lap = laplace(image.copy())
        lap -= lap.min()
        lap /= lap.max()
        lap[lap< lap.mean()*1.8]=lap.min()
        label_im, nb_labels = label(lap)
        laps = remove_small_objects(label_im,3)
        label_im, nr_labels = label(laps)
        for n in range(1, nr_labels+1):
            cx,cy = center_of_mass(label_im==n)
               
            #refine center by checking its position on large frame
            ccx, ccy =  cx,cy#int(numpy.round(cx)),int(numpy.round(cy) )


            xmin,xmax = max(0,int(ccx*fact)-ss), int(ccx*fact)+ss
            ymin,ymax = max(0,int(ccy*fact)-ss), int(ccy*fact)+ss
            cutout2 = (frames_full[imnr][xmin:xmax,ymin:ymax]).copy()
            cutout = cutout2-cutout2.min()
            cutout /= cutout.max()
            cutout = gaussian_filter(cutout,2)
            l = laplace((cutout.max()-cutout)/cutout.max())
           
            gfl = gaussian_filter((l.max()-l ),5)
            gfl[gfl<numpy.median(gfl)*1.1] = numpy.median(gfl)*1.1
            gfl -= gfl.min()
            if gfl.max(): gfl /= gfl.max()
            local_maxi = peak_local_max(gfl, indices=False, footprint=numpy.ones((5, 5)) )
            lm = local_maxi.copy()
            #imshow(lm+cutout2*(lm==0))
            #show()
            tot = local_maxi.sum()
            if tot ==0:
                sum2+=1
                #if abs(cy-226) < 8: print imnr,cx,cy
                add = True
                for a,b,c in list_cx_cy_imnr:
                    if abs(cx- a)+abs(cy-b) < 4 and imnr*num_procs+proc_id ==c:
                        add = False
                        break
                if add and cx*fact+16 < frames_full[0].shape[1] and cy*fact+16 < frames_full[0].shape[0]:
                    list_cx_cy_imnr.append([float(cx),float(cy),imnr*num_procs+proc_id])
                

                
            else:
                for i in range(tot):
                    x,y = cutout.shape

                    cx,cy = (local_maxi.argmax() )/y,(local_maxi.argmax())%y
                    rx,ry = cx+ccx*fact-ss,cy+ccy*fact-ss
                    if cropsize > cutout.shape[1] and ccy+ss < cutout.shape[1]: ry += cropsize - cutout.shape[1]
                    if cropsize > cutout.shape[0] and ccx+ss < cutout.shape[0]: rx += cropsize - cutout.shape[0]
                    
                    # sometimes the algorithm gives two adjacent peaks. We do not keep the second.
                    add = True
                    for a,b,c in list_cx_cy_imnr:
                        if abs((rx)/fact - a)+abs((ry)/fact-b) < 4 and imnr*num_procs+proc_id ==c:
                            add = False
                            break
                    if add and cx+16 < frames_full[0].shape[1] and cy+16 < frames_full[0].shape[0]:
                        dcx,dcy = (rx)/max(fact,1),(ry)/max(fact,1)
                        
                        #if abs(dcy-226) < 8:print imnr, dcx,dcy
                        list_cx_cy_imnr.append([float(dcx),float(dcy),imnr*num_procs+proc_id])
                        #frames_full[imnr][ int((cx+(ccx*fact)-ss))][int((cy+(ccy*fact)-ss))] = frames_full[imnr].max()
                        #fiducials[min(int(numpy.round(dcx)),y),min(int(numpy.round(dcy)-1),x)] = 1
                        #if dcx < 20:
                        #print dcx,dcy
                    else:
                        lm[int(cx)][int(cy)] = 0
                    local_maxi[int(cx)][int(cy)] = False
                    #print local_maxi
               
                #fiducials[max(0,int(ccx*fact)-ss):int(ccx*fact)+ss,max(0,int(ccy*fact)-ss):int(ccy*fact)+ss] += lm
                #cf[max(0,int(ccx*fact)-ss):int(ccx*fact)+ss,max(0,int(ccy*fact)-ss):int(ccy*fact)+ss] += lm
        
        
                # append the peak to the peaklist
        
    #print len(list_cx_cy_imnr), (fiducials>0).sum()
    #print sum2
    if sen==False: fid_list += list_cx_cy_imnr
    return numpy.zeros_like(frames_full[0]), list_cx_cy_imnr

#def find_potential_fiducials          (frames, frames_full, bin, bin_full, cropsize=50,target=0,fid_list=[],proc_id=0,num_procs=1):
def find_potential_fiducials_sensitive(frames, frames_full,  bin, bin_full, cropsize=32,target=0,fid_list=[],proc_id=0,num_procs=1,aa=20):
    
    a, l1 = find_potential_fiducials(frames, frames_full, bin, bin_full, cropsize, target, fid_list, proc_id, num_procs,sen=True)
    
    list_cx_cy_imnr = []
    fact = 1.*bin/bin_full
    for imnr in range(len(frames_full)):

        ds = frames[imnr]  #downsample(wiener(full_image),n)**0.7
        ds2 = ds.copy()
        for xx in range(0,len(ds)):
            for yy in range(0,len(ds[0])):
                dd = median(ds[max(0,xx-10):xx+10,max(0,yy-10):yy+10])
                           

                ds2[xx,yy] = ds[xx,yy]-dd
        ds = ds2

        minf = minimum_filter( gaussian_filter(ds - gaussian_filter(ds, 5),.1) , footprint=ones((1,1)) )

        hpf = minf - gaussian_filter(minf,4)
        hpf -= hpf.min()
        hpf /= median(hpf)
        hpf = 1-hpf

        gg = ((gaussian_filter(hpf,.1) )  )
        gg -= gg.min()
        gg /= median(gg)

        rr = gaussian_filter( laplace( gaussian_filter( (1- gg/gg.max()), .1) ) , 1)
        rr -= rr.min()
        rr /= rr.max()
        rr += 0.0001

        tot = zeros_like(rr)
        for ii in arange(4.1,3.59,-0.05):
            rrr = rr.copy()
            tot += testing_done(rrr,aa=aa,ff=ii)

            ll, nn = label(tot)

            
            for N in range(1,nn+1):
                cx,cy = center_of_mass(ll==N)
                #TOT[int(round(cx)),int(round(cy))] = 1                                                                                                                                                            
                add = True
                for a,b,c in list_cx_cy_imnr:
                    if abs(cx- a)+abs(cy-b) < 4 and imnr*num_procs+proc_id ==c:
                        add = False
                        break
                if add:
                    #dcx,dcy = (cx+(ccx*fact)-ss)/max(fact,1),(cy+(ccy*fact)-ss)/max(fact,1)                                                                                                                       
                    list_cx_cy_imnr.append([cx,cy,imnr*num_procs+proc_id])

            '''
            list_cx_cy = []
            for N in range(1,nn+1):
                cx,cy = center_of_mass(ll==N)

                list_cx_cy += refine(cx,cy,frames_full[imnr],fact,cropsize)
                #TOT[int(round(cx)),int(round(cy))] = 1
                add = True

                for cx,cy in list_cx_cy:
                    add = True
                    for a,b,c in list_cx_cy_imnr:
                        if abs(cx- a)+abs(cy-b) < 4 and imnr*num_procs+proc_id == c:
                            add = False
                            break
                    if add:
                    #dcx,dcy = (cx+(ccx*fact)-ss)/max(fact,1),(cy+(ccy*fact)-ss)/max(fact,1)
                        list_cx_cy_imnr.append([cx,cy,imnr*num_procs+proc_id])
            '''


    out = l1
    for rx,ry,im in list_cx_cy_imnr:
        #print rx,ry,imnr
        points = refine(rx, ry, frames_full[imnr], imnr*num_procs+proc_id, fact, cropsize)
        #print 'points: ', points
        for cx,cy,ims in points:
            add = True
            for a,b,c in out:
                if abs(cx- a)+abs(cy-b) < 4 and imnr*num_procs+proc_id == c:
                    add = False
                    break
            if add:
                out.append([cx,cy,ims])

    fid_list += out#list_cx_cy_imnr


    out

    return list_cx_cy_imnr, tot, ds

def refine(cx, cy, frame, imnr, fact, cropsize):
    '''This function finds the maxima in a 2D frame'''
    if len(frame.shape) != 2: raise Exception('wrong shape of input frame')
    try: 
        float(cx)
        float(cy)
    except:
        raise Exception('Please provide float or int as input coordinates')
    list_cx_cy = []

    ss = cropsize//2
    
    ccx, ccy =  cx,cy#int(numpy.round(cx)),int(numpy.round(cy) )                                                                                                                                         

    cutout2 = (frame[max(0,int(ccx*fact)-ss):int(ccx*fact)+ss,max(0,int(ccy*fact)-ss):int(ccy*fact)+ss]).copy()
    cutout = cutout2-cutout2.min()
    cutout /= cutout.max()
    cutout = gaussian_filter(cutout,2)
    l = laplace((cutout.max()-cutout)/cutout.max())

    gfl = gaussian_filter((l.max()-l ),5)
    gfl[gfl<numpy.median(gfl)*1.1] = numpy.median(gfl)*1.1
    gfl -= gfl.min()
    if gfl.max(): gfl /= gfl.max()
    local_maxi = peak_local_max(gfl, indices=False, footprint=numpy.ones((5, 5)) )
    lm = local_maxi.copy()
    #fig,ax = subplots(1,1)
    #ax.imshow(lm+cutout2*(lm==0))                                                                                                                                                                          
    #show()                                                                                                                                                                                              
    tot = local_maxi.sum()
    if tot ==0:
        list_cx_cy.append( [cx,cy,imnr] )
    else:
        
        for i in range(tot):
            x, y = cutout.shape

            cx, cy = (local_maxi.argmax() )/y,(local_maxi.argmax())%y
            rx, ry = cx+ccx*fact-ss,cy+ccy*fact-ss
            if cropsize > cutout.shape[1] and ccy+ss < cutout.shape[1]: ry += cropsize - cutout.shape[1]
            if cropsize > cutout.shape[0] and ccx+ss < cutout.shape[0]: rx += cropsize - cutout.shape[0]

            # sometimes the algorithm gives two adjacent peaks. We do not keep the second.                                                                                                               
            add = True
            for a,b,c in list_cx_cy:
                if abs((rx)/fact - a)+abs((ry)/fact - b) < 4:
                    add = False
                    break
            if add:
                dcx,dcy = (rx)/max(fact,1),(ry)/max(fact,1)
                list_cx_cy.append([dcx,dcy,imnr])
                
    return list_cx_cy




def testing_done(rr,aa=20,ff=4):

    d1 = rr.copy()
    d2 = rr.copy()

    edge = zeros_like(rr)

    for xx in range(0,len(rr),aa):
        for yy in range(0,len(rr[0]),aa):

            dd = rr[xx:xx+aa,yy:yy+aa].copy()
            if dd.std() >0.7:
                #print 'add', std(dd), dd.min(), dd.max()                                                                                                                                  
                edge[xx:xx+aa,yy:yy+aa] = 1
            #else:                                                                                                                                                                         
                #print 'tast', std(dd), dd.min(), dd.max()                                                                                                                                 

            dd[dd < median(dd)+std(dd)*ff] = 0
            d1[xx:xx+aa,yy:yy+aa] = dd

    ll1,nn1 = label(remove_small_objects(d1>0,2))

    for xx in range(-aa//2,len(rr),aa):
        for yy in range(-aa//2,len(rr[0]),aa):

            dd = rr[max(0,xx):xx+aa,max(0,yy):yy+aa]
            if dd.min() > 0.2:
                #print dd.min()                                                                                                                                                            
                edge[xx:xx+aa,yy:yy+aa] = 1
                #print median(dd),median(dd)+std(dd)*2, median(dd)*1.2                                                                                                                     
            dd[dd < median(dd)+std(dd)*ff] = 0
            d2[max(0,xx):xx+aa,max(0,yy):yy+aa] = dd

    ll2,nn2 = label(remove_small_objects(d2>0,2))



    return ((ll2+ll1) >0)


def read_data(fnames,bin,bin_full):
    
    frames = []
    frames_full = []
    for fname in fnames:
            #data  = read('{}'.format(str(fname)),binning=[self.binning.get(),self.binning.get(),1])
        dataf  = read('{}'.format(str(fname)),binning=[bin_full,bin_full,1])

        w = wiener(dataf)
        hpf = dataf - gaussian_filter(dataf,10)*0
        data = mrcOperations.downsample( w+hpf, int(round(bin*1./bin_full)))
            #data = vol2npy(f).copy()
            #data = read_mrc(fname)
        data -= data.min()
        data /= data.max()

        dataf -= dataf.min()
        dataf /= dataf.max()
            
        frames.append( (data)**0.75 )
        frames_full.append( wiener(dataf) )
        
        
    return frames, frames_full


if __name__ == '__main__':
    import glob
    from matplotlib.patches import Circle
    from matplotlib.collections import PatchCollection
    from time import time 
    from find_fiducials_refined import get_positions_potential_fiducials
    
    cmap = ['r','purple','orange','lightgreen','yellow','cyan','magenta','b','pink','lightgreen','violet','brown']
    
    exp = "{}/03_Tomographic_Reconstruction/tomogram_{:03d}/sorted/".format(sys.argv[1].strip('/'), int(sys.argv[2]))
    bin, bin_full = int(sys.argv[3]),int(sys.argv[4])
    fnames = numpy.sort( numpy.array( [ exp+line.split('/')[-1] for line in os.listdir( exp ) if line.endswith('mrc') ] ) )

    s = time()
    frames, frames_full = read_data(fnames[:],bin,bin_full)
    print ("Read time: ", time()-s)
    s = time()

    fiducials, list_cx_cy_imnr = find_potential_fiducials(frames, frames_full, bin, bin_full)

    #print list_cx_cy_imnr
    #list_cx_cy_imnr = []
    #for i in range(len(frames)):
    #    list_cx_cy_imnr,fiducials,ds2 = get_positions_potential_fiducials(frames_full[i],i,list_cx_cy_imnr,n=bin/bin_full,aa=20)
        
    print ("Find Fid time: ", time()-s)
    s = time()           
    coordinates, mark_frames, index_map, frame_shifts = index_potential_fiducials(fnames, list_cx_cy_imnr,plot=True)
    print ("Index time: ", time() -s)
    
    coord_orig = coordinates-frame_shifts
    mark_frames_orig = mark_frames-frame_shifts
    
    t = len(frames)
    c =8
    r= max(2,int(numpy.ceil(t/(c*1.))))
    fig,ax = subplots(r,c,figsize=(c*1.8,r*1.8))
    for nr in range(r):
        for nc in range(c):
            ax[nr][nc].get_xaxis().set_visible(False)
            ax[nr,nc].get_yaxis().set_visible(False)
    
    fig.tight_layout(pad=0)
    subplots_adjust(left=0, right=1, top=1, bottom=0)

    for imnr in range(len(frames)):
        rr,cc = imnr/c,imnr%c
        ax[rr][cc].imshow(frames_full[imnr],cmap='gray')
        for n in range(len(mark_frames_orig[0])):
            mp = mark_frames_orig[imnr][n]
            #print imnr, n, index_map[imnr][n],mark_frames[imnr][n][0],mark_frames[imnr][n][1]
            #for index in range(len(coordinates[0])):
            #p = coord_orig[imnr][index]
            if mp.sum() != 0 and mp.astype(int) in coord_orig.astype(int):
                if cmap[(index_map[imnr][n]%9)] == 'r': continue
                circle = Circle((int( (mp[1]*bin/bin_full) ), int( (mp[0]*bin/bin_full)) ), radius=5, facecolor='none',edgecolor=cmap[(index_map[imnr][n]%9)])
                ax[rr][cc].add_patch(circle)
                #frames_full[imnr][round(p[0]*4.),round(p[1]*4.)] = frames_full[imnr].max()
    

            elif mp.sum() !=0:
                
                circle = Circle((int( (mp[1]*bin/bin_full) ), int( (mp[0]*bin/bin_full)) ), radius=5, facecolor='none',edgecolor='red')
                #ax[rr][cc].add_patch(circle)
    '''
    for imnr in range(len(frames)):
        rr,cc = imnr/c,imnr%c
        for index in range(len(mark_frames[0])):
            p = mark_frames_orig[imnr][index]

            if abs(mark_frames[imnr][index][1] - 226) < 5:
                print p[1],p[0]
            
            if p.sum() != 0:
                
                if not p in coord_orig:
                    circle = Circle((int( (p[1]*4.) ), int( (p[0]*4.)) ), radius=5, facecolor='none',edgecolor='red')
                    ax[rr][cc].add_patch(circle)
                #frames_full[imnr][round(p[0]*4.)-1:round(p[0]*4)+2,round(p[1]*4.)-1:round(p[1]*4.)+2] = frames_full[imnr].max()
            
        ax[rr][cc].imshow(frames_full[imnr],cmap='gray')
    '''        
    #plt.plot((fiducials).sum(axis=0))
    fig.canvas.draw()
    show()
