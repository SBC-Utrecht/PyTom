import numpy as np
from pytom.gui.mrcOperations import downsample, read_mrc as read
from scipy.ndimage import gaussian_filter, laplace, label, gaussian_laplace, center_of_mass
from scipy.ndimage.filters import minimum_filter
from scipy.signal import wiener
from skimage.morphology import remove_small_objects
from skimage.feature import peak_local_max


class PickingFunctions():
    def __init__(self):
        pass

    def diff(self, vec0, vec1):
        return np.linalg.norm(vec0-vec1)

    def closest_point(self, ref_vec, points, cutoff=3.4, it=0, alpha=270):

        if alpha < 0.: alpha += 360

        c = np.zeros(len(points), dtype=float)
        for n, p in enumerate(points):
            c[n] = self.diff(ref_vec, p)
            if abs(self.distance(p, ref_vec, alpha))> cutoff:
                c[n] = 1000
        try:
            for i in range(it): c[c.argmin()] = 10000
            return c.argmin(), c.min()
        except:
            return -1

    def dist_z(self, y0,y1,a1):
        z1 = (y1-y0*np.cos(a1) )/ (np.sin(a1))
        return z1


    def calc_matrix(self, cur, ref, cutoff, alpha):
        alpha = alpha % 360
        distance = np.ones((len(cur), len(ref))) * 9999
        for n, p in enumerate(cur):
            for m, ref_vec in enumerate(ref):
                if abs(self.distance2(p, ref_vec, alpha, cutoff)) <= cutoff:
                    distance[n][m] = self.diff(p, ref_vec)

        return distance

    def cp_it(self, dist, it):
        dd = dist.copy()
        for i in range(it):
            dd[dd.argmin()] = 9999
        return dd.argmin(), dd.min()

    def compare(self, cur, ref, imnr, cutoff=3.4, v=False, tiltangles=[], y0=0, a0=0, pos=1, cutoff_dist=10,
                tiltaxis=0):
        angle = tiltaxis
        dist_map = [1000] * (ref.shape[0])
        index_map = np.zeros((ref.shape[0]), dtype=int)
        coordinates = np.zeros_like(ref)
        #distance_matrix = calc_matrix(cur, ref, cutoff, tiltaxis)
        for i, p in enumerate(cur):
            if p.sum() == 0: continue
            index, dist = self.closest_point(p, ref, cutoff, alpha=angle)

            if dist < 4 * cutoff_dist:
                index_map[i] = index + 1

                tgt = 4
                # if index == tgt :
                #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
                #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]
                # print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[index][0], p[0], d_theta*np.pi/180.) ), p[0]
                # print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+2, imnr, int(np.round(d_theta)), int(np.round(tiltangles[imnr])), dist_z(ref[tgt+1][0], p[0], d_theta*np.pi/180.) ), p[0]

                # if index == tgt+1 :
                #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
                #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]
                # print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[tgt][0], p[0], d_theta*np.pi/180.) ), p[0]
                # print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[index][0], p[0], d_theta*np.pi/180.) ), p[0]

        count = 0
        # print

        while count < 10 and len(index_map) != len(np.unique(index_map)) + (1 * (index_map == 0)).sum() - 1:
            # if 1: print imnr, index_map

            for ind in np.unique(index_map):
                if ind == 0: continue

                closest, second = [], []
                for n, index in enumerate(index_map):
                    if index == 0 or index != ind or cur[n].sum() == 0: continue
                    for ii in range(len(index_map)):
                        ind1, dis1 = self.closest_point(cur[n], ref, it=ii, alpha=angle)
                        if ind1 + 1 == ind: break
                    ind2, dis2 = self.closest_point(cur[n], ref, it=ii + 1, alpha=angle)
                    closest.append([n, ind1 + 1, dis1, ind2 + 1, dis2])

                candidates = []
                m = [0, 0, 9999, 0, 9999]
                for c in closest:
                    if c[2] < cutoff_dist:
                        candidates.append(c)
                    if c[2] < m[2]:
                        m = c

                # if v and len(candidates) > 1: print candidates
                if len(candidates) > 2:

                    for c in closest:
                        if m[0] == c[0]: continue

                        if c[4] < 100:
                            if v: print ('changed', index_map[c[0]], index_map[c[1]])
                            index_map[c[0]] = c[3]
                        else:
                            index_map[c[0]] = 0

                elif len(candidates) == 2:
                    todo = True
                    c1, c2 = candidates
                    dist, dis, d1, d2 = c1[2], c1[4], c2[2], c2[4]

                    std2, std1 = np.std([d1, dis]), np.std([d2, dist])

                    # if imnr == 15:
                    #    print c1
                    #    print c2

                    if std2 < 4 and std2 < 4:
                        p1, p2 = cur[c1[0]], cur[c2[0]]

                        diff = self.distance(p1, p2, angle)

                        # diff = abs(p1[1]*4.-p2[1]*4.)
                        # print diff
                        if abs(diff) > 3.:
                            # print p1, p2
                            if diff < 0.:  # p1[1] > p2[1]:
                                m = c1
                                index_map[c2[0]] = c2[3]
                            else:
                                index_map[c1[0]] = c1[3]
                                m = c2
                            # print imnr, index_map
                        else:
                            if std1 < std2:
                                m = c1
                                index_map[c2[0]] = c2[3]
                            else:
                                index_map[c1[0]] = c1[3]
                                m = c2

                    elif np.std([d1, dis]) < np.std([d2, dist]) and d2 < cutoff_dist and dis < cutoff_dist:
                        index_map[c1[0]] = c1[3]
                        m = c2
                        if v: print (count, 'option a')

                    elif np.std([d1, dis]) < np.std([d2, dist]) and d2 < cutoff_dist and dis < cutoff_dist:
                        index_map[c2[0]] = c2[3]
                        m = c1
                        if v: print (count, 'option b')
                        todo = False
                    elif d1 > dist and dist < 100:
                        m = c1
                        if d1 < 100:
                            index_map[c2[0]] = c2[3]
                        else:
                            index_map[c2[0]] = 0
                        if v: print (count, 'option c')

                    elif d1 <= dist and d1 < 100:
                        m = c2
                        if dist < 10 and dis < 20:
                            index_map[c1[0]] = c1[3]
                        else:
                            index_map[c1[0]] = 0
                        if v: print (count, 'option d', c1, c2, todo)

                    else:
                        if v: print (count, 'no_change')

                if len(closest) > 1:
                    for c in closest:
                        if m[0] == c[0]: continue
                        if c[4] < cutoff_dist:
                            index_map[c[0]] == c[3]
                        else:
                            index_map[c[0]] = 0

            count += 1
            if v: print (imnr, index_map)

        for n, p in enumerate(cur):
            if index_map[n]: coordinates[index_map[n] - 1] = p

        return index_map, coordinates

    def compare_fast(self, cur, ref, imnr, cutoff=3.4,v=False,tiltangles=[],y0=0,a0=0,pos=1,cutoff_dist=10,tiltaxis=0):
        angle = tiltaxis
        dist_map = [1000]*(ref.shape[0])
        index_map = np.zeros((ref.shape[0]),dtype=int)
        coordinates = np.zeros_like(ref)
        distance_matrix = self.calc_matrix(cur, ref, cutoff, tiltaxis)
        for i,p in enumerate(cur):
            if p.sum() == 0: continue
            # index,dist = self.closest_point(p,ref,cutoff,alpha=angle)

            if distance_matrix[i,:].min()< cutoff_dist:

                index_map[i] = distance_matrix[i, :].argmin()+1

                tgt = 4
                #if index == tgt :
                #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
                #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]
                    #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[index][0], p[0], d_theta*np.pi/180.) ), p[0]
                    #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+2, imnr, int(np.round(d_theta)), int(np.round(tiltangles[imnr])), dist_z(ref[tgt+1][0], p[0], d_theta*np.pi/180.) ), p[0]

                #if index == tgt+1 :
                #    if pos == 1: d_theta = tiltangles[imnr] - tiltangles[imnr-1]
                #    else: d_theta = tiltangles[imnr]-tiltangles[imnr+1]
                    #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(tgt+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[tgt][0], p[0], d_theta*np.pi/180.) ), p[0]
                    #print '{:4d} {:4d} {:4d} {:4d} {:6.0f}'.format(index+1, imnr, int(np.round(p[0])), int(np.round(p[1])), dist_z(ref[index][0], p[0], d_theta*np.pi/180.) ), p[0]

        count = 0
        #print

        while count < 10 and len(index_map) != len(np.unique(index_map))+(1*(index_map==0)).sum()-1:
            #if 1: print imnr, index_map

            for ind in np.unique(index_map):
                if ind == 0: continue

                closest, second = [], []
                for n,index in enumerate(index_map):
                    if index==0 or index != ind or cur[n].sum() == 0: continue
                    ii=0
                    ind1, dis1 = 0, 9999
                    for ii in range(len(index_map)):
                        ind1,dis1 = self.cp_it(distance_matrix[n,:], ii)
                        if ind1+1==ind: break
                    ind2,dis2 = self.cp_it(distance_matrix[n,:], ii+1)
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
                            index_map[c[0]]=c[3]
                        else: index_map[c[0]] = 0

                elif len(candidates) ==2:
                    todo=True
                    c1,c2 =candidates
                    dist,dis,d1,d2 = c1[2],c1[4],c2[2],c2[4]

                    std2,std1 = np.std( [d1,dis] ), np.std( [d2,dist] )

                    #if imnr == 15:
                    #    print c1
                    #    print c2

                    if std2<4 and std2 <4:
                        p1,p2 = cur[c1[0]],cur[c2[0]]

                        diff = self.distance(p1,p2,angle)

                        #diff = abs(p1[1]*4.-p2[1]*4.)
                        #print diff
                        if abs(diff)>3.:
                            #print p1, p2
                            if diff < 0. : #p1[1] > p2[1]:
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

                    elif np.std( [d1,dis] ) < np.std( [d2,dist] ) and d2 < cutoff_dist and dis < cutoff_dist:
                        index_map[c1[0]] = c1[3]
                        m = c2
                        if v: print (count,'option a')

                    elif np.std( [d1,dis] ) < np.std( [d2,dist] ) and d2 < cutoff_dist and dis < cutoff_dist:
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
                        if c[4] < cutoff_dist:
                            index_map[c[0]]==c[3]
                        else: index_map[c[0]] = 0

            count +=1
            if v: print (imnr, index_map)

        for n,p in enumerate(cur):
            if index_map[n]:
                coordinates[index_map[n]-1] = p

        return index_map, coordinates

    def distance2(self, p0, p1, angle, cutoff):
        if p0.sum() == 0: return 9999
        dx = p1[0]-p0[0]
        dy = p1[1]-p0[1]
        a = np.arctan2(dy, dx)
        v = np.sqrt(dx**2+dy**2)
        diff_angle = a-angle*np.pi/180.
        distx = v*np.cos(diff_angle)
        disty = v * np.sin(diff_angle)
        if distx > cutoff: return 9999
        return distx

    def distance(self, p0, p1, angle):
        dx = p1[0]-p0[0]
        dy = p1[1]-p0[1]
        a = np.arctan2(dy, dx)
        v = np.sqrt(dx**2+dy**2)
        diff_angle = a-angle*np.pi/180.
        dist = v*np.cos(diff_angle)
        return dist

    def sort(self, obj, nrcol ):
        obj.sort(key=lambda i: float(i[nrcol]))


    def cross_correlate(self, im0, im1):
        for i in range(len(im0.shape)):
            assert im0.shape[i] == im1.shape[i]
            assert im0.shape[i] > 0

        ft_im0 = (np.fft.fftn(im0))
        ft_im1 = (np.fft.fftn(im1))
        cc = np.fft.fftshift(abs(np.fft.ifftn((ft_im0 * ft_im1))))
        return cc

    def detect_shifts_many(self, mark_frames, max_shiftx=5.5,max_shifty=2.,sd=1,diag=False, image=[]):
        from time import time

        diag=False
        if diag: s=time()

        tilt_nr,mark_nr, k = mark_frames.shape
        dimx,dimy = image.shape
        factor = max(1,np.ceil(dimx / 1024.))
        factor = 1
        out = np.zeros((dimx//factor,dimy//factor,tilt_nr))

        frame_shifts = np.zeros((tilt_nr,mark_nr,2),dtype=float)
        fs = np.zeros((tilt_nr,2),dtype=float)
        outline_detect_shifts = ''
        num_shifts = 0

        for itilt in range(tilt_nr):
            for imark in range(mark_nr):
                x,y = mark_frames[itilt][imark]/factor
                #if x> dimx-1 or y> dimy-1: print x,y
                if x > 0 and y > 0 and x < dimx-1 and y< dimy-1: out[int(round(x))][int(round(y))][itilt] = 1
        for itilt in range(tilt_nr-1):

            c, shifty, shiftx, m = self.detect_shift(gaussian_filter(out[:,:,itilt],sd), gaussian_filter(out[:,:,itilt+1],sd),image=image)
            shiftx *= factor
            shifty *= factor
            if m > 0.001:#abs(shiftx) > max_shiftx or abs(shifty) > max_shifty:
                print ("\t{:2d}-{:2d} {:4.0f} {:4.0f} {:4.3f}".format(itilt,itilt+1,shiftx,shifty,m))
                frame_shifts[itilt+1:] += [ shiftx, shifty]
                fs[itilt+1:] += [ shiftx, shifty]
                num_shifts +=1
                outline_detect_shifts += "\t{:2d}-{:2d} {:4.0f} {:4.0f} {:4.3f}\n".format(itilt,itilt+1,shiftx,shifty,m)

        frame_shifts *= (mark_frames.sum(axis=-1) > 0)[:,:,np.newaxis]
        if diag: print ('detect frame shifts: {:.0f} msec'.format(1000*(time()-s)))

        return frame_shifts, num_shifts, outline_detect_shifts, fs

    def detect_shift(self, arr0,arr1,image=[]):
        x,y = image.shape
        cross = abs(np.fft.fftshift( np.fft.ifft2(np.fft.fftshift(np.fft.fftshift(
            np.fft.fft2(arr0))*np.conj(np.fft.fftshift(np.fft.fft2(arr1)))))))**2
        locx,locy =  (abs((cross))**2).flatten().argmax()%y, (abs((cross))**2).flatten().argmax()/y
        return cross, locx-y/2, locy-x/2, cross[int(locy)][int(locx)]

    def detect_shifts_few(self, mark_frames, maxshift=9, max_shiftx=5.5, max_shifty=2.,diag=False, image=[]):
        from time import time
        if diag: s = time()

        x,y,z = mark_frames.shape

        distance_matrix = 1000*np.ones((x-1,y),dtype=float)
        shift = np.zeros_like(mark_frames)
        size_shift = np.ones((x-1,y),dtype=float)

        for ii in range(x-1):
            distance_matrix[ii,:], shift[ii,:,:], size_shift[ii,:] = self.pair(mark_frames[ii],mark_frames[ii+1],pr=(ii==8))
            #print ii, distance_matrix[ii]
        frame_shifts = np.zeros((x,y,2),dtype=float)
        fs = np.zeros((x,2),dtype=float)
        num_shifts = 0

        outline_detect_shifts = ''
        for iii in range(x-1):
            sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]<150)
            d_shift = size_shift[iii][sel]
            shifty, shiftx = np.median(shift[iii,:,1][sel]), np.median(shift[iii,:,0][sel])
            #if iii == 31: print sorted(d_shift)
            #print iii, np.std(d_shift), d_shift

            n, h = np.histogram(d_shift, bins=10, range=None, normed=False, weights=None, density=None)
            #if iii == 31:
            #    print n, h
                #print np.std(h[n.argmax():n.argmax()+1]
            #distance_matrix[iii], d_shift, shift[iii]
            if abs(shifty) > max_shifty or abs(shiftx) > max_shiftx :
                if np.std(d_shift) > 5:
                    sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]<85)
                    d_shift = size_shift[iii][sel]
                    shifty, shiftx = np.median(shift[iii,:,1][sel]), np.median(shift[iii,:,0][sel])
                    #print iii, np.std(d_shift), d_shift
                    #print

                if np.std(d_shift) > 5:
                    #diff = h[-1]-h[0]
                    #total = len(d_shift)
                    sel = (distance_matrix[iii]<maxshift)*(size_shift[iii]< h[min(len(n),n.argmax()+1)] )*(size_shift[iii]>h[n.argmax()])
                    d_shift = size_shift[iii][sel]
                    shifty, shiftx = np.median(shift[iii,:,1][sel]), np.median(shift[iii,:,0][sel])


                if np.std(d_shift) > 5:
                    continue
                outline_detect_shifts += "\t{:2d}-{:2d} {:4.1f} {:4.1f} {:4.1f} | {:4.0f} {:4.0f}\n".format(iii,iii+1, np.median(d_shift ), np.std(d_shift), np.mean(d_shift), shiftx,shifty)
                frame_shifts[iii+1:] += [ shiftx, shifty]
                fs[iii+1:] += [ shiftx, shifty]
                num_shifts +=1

        frame_shifts *= (mark_frames.sum(axis=-1) > 0)[:,:,np.newaxis]
        print ("encountered {} frame shifts.".format(num_shifts))
        print (outline_detect_shifts   )

        if diag: print ('detect frame shifts: {:.0f} msec'.format(1000*(time()-s)))

        return frame_shifts, num_shifts, outline_detect_shifts,fs

    def pair(self, l1,l2,pr=0):
        len_l1,len_l2 = len(l1), len(l2)
        com1 = np.zeros((len_l1,len_l1,2))
        com2 = np.zeros((len_l2,len_l2,2))

        for i in range(len_l1):
            com1[i,:,:] = l1-l1[i]
        for j in range(len_l2):
            com2[j,:,:] = l2-l2[j]


        dl = np.ones((len_l1))*1000
        sl = np.ones_like(l1)*1000
        nl = dl.copy()

        for n, el in enumerate(com1.reshape(len_l1**2,2)):
            if el.sum() == 0:
                continue
            diff = np.linalg.norm(com2-el,axis=-1)

            n2,distance=diff.argmin(),diff.min()

            x1,y1=n//len_l1, n%len_l1
            x2,y2=n2//len_l1, n2%len_l1

            shift = (l1[x1]+l1[y1])/2-(l2[x2]+l2[y2])/2
            norm = np.linalg.norm(shift)


            if dl[x1] > distance and l1[x1].sum() and l1[y1].sum() and l2[x2].sum() and l2[y2].sum():
                dl[x1] = distance
                sl[x1] = shift
                nl[x1] = norm

        return dl,sl,nl


    def frames_without_indexed_fiducials(self, c,tiltangle):
        empty_frames = []
        try:
            a, b = int((c[:,0][c[:,0]>0]).mean()), int((c[:,1][c[:,1]>0]).mean())
        except:
            empty_frames.append(tiltangle)
        return empty_frames


    def index_potential_fiducials(self, frames, mark_frames, frame_shifts, plot=True, cut=3, zero_angle=20, tiltangles=[],
                                  add_marker=False, excluded=[], user_coords=[],diag=False,tiltaxis=0):
        from time import time
        s = time()
        c = 300
        r = len(frames)
        cmap = ['r','orange','lightgreen','yellow','b','cyan','magenta','pink','violet','black']


        mark_frames += frame_shifts.astype(float)

        coordinates_sorted = np.zeros_like(mark_frames)
        index_map = -1*np.ones((r,len(mark_frames[0])),dtype=int)

        num_images = 100
        cur = 1
        for i in range(len(mark_frames[zero_angle])):
            if mark_frames[zero_angle][i][0] > frame_shifts[zero_angle][i][0]-1.5:
                index_map[zero_angle][i] = cur
                coordinates_sorted[zero_angle][i] = mark_frames[zero_angle][i]
                cur +=1

        for tiltangle in range(zero_angle+1,min(len(frames),zero_angle+1+num_images)):
            cur_fid = mark_frames[tiltangle,:].copy()
            ref_fid = coordinates_sorted[tiltangle-1,:].copy()
            c =ref_fid[:,:]

            for temp_i in np.arange(tiltangle-1,zero_angle,-1):
                markers = coordinates_sorted[temp_i-1,:].copy()
                absent = (((ref_fid == [0,0]).sum(axis=1) ==2)*1)[:,np.newaxis]*markers
                ref_fid += absent

            #print cur_fid.T
            if len(user_coords):
                for n,(x,y) in enumerate(user_coords[tiltangle]):
                    if x>0 and y>0:
                        print ('found user selected marker: {:6.1f},{:6.1f}'.format(x,y))
                        ref_fid[n] = [x+.1,y+.1]


            index_cur_fid, coordinates_cur_sorted = self.compare_fast(cur_fid, ref_fid, tiltangle, cutoff=cut,tiltangles=tiltangles,pos=1,tiltaxis=tiltaxis)
            index_map[tiltangle,:] = index_cur_fid
            coordinates_sorted[tiltangle,:] = coordinates_cur_sorted


        for tiltangle in np.arange(zero_angle-1,max(-1,zero_angle-1-num_images),-1):
            cur_fid = mark_frames[tiltangle,:].copy()
            ref_fid = coordinates_sorted[tiltangle+1,:].copy()

            for temp_i in np.arange(tiltangle+1,zero_angle,+1):
                markers = coordinates_sorted[temp_i+1,:]
                absent = (((ref_fid == [0,0]).sum(axis=1) ==2)*1)[:,np.newaxis]*markers
                ref_fid += absent

            index_cur_fid, coordinates_cur_sorted = self.compare_fast(cur_fid, ref_fid, tiltangle, cutoff=cut, tiltangles=tiltangles, pos=0, tiltaxis=tiltaxis)
            index_map[tiltangle,:] = index_cur_fid
            coordinates_sorted[tiltangle,:] = coordinates_cur_sorted
        #frames_without_indexed_fiducials(coordinates_sorted[tiltangle,:], tiltangle)

        #print coordinates_sorted
        printline = ''
        cntr = 0
        tot = 0
        listdx= []
        #for nr in np.unique(index_map)[:]:

        outc = np.zeros_like(coordinates_sorted)
        imark = 0
        for i in range(len(frames)):
            if excluded[i]:
                coordinates_sorted[i,:,:] *= 0

        # self.assignedFiducials[self.assignedFiducials > -0.5] = 0

        for nr in range(coordinates_sorted.shape[1]):
            if nr+1 > 0.5:# and (index_map == nr).sum() > 2:
                cntr+=1
                cur = 0

                ids_assigned = []
                for j, [mx,my] in enumerate(coordinates_sorted[:,nr,:]):
                    #if mx >0 and my >0: cur += 1
                    tempn = 0
                    for mcnt, (markx,marky) in enumerate(mark_frames[j,:,:]):

                        if abs(mx-markx) < 0.1 and abs(marky-my) < 0.1:

                            if (markx-frame_shifts[j,mcnt,0]) > 0.01 and (marky-frame_shifts[j,mcnt,1] ) > 0.01:
                                #print (nr,j, marky, marky-frame_shifts[j,mcnt,1])
                                cur += 1
                                coordinates_sorted[j,nr,:] -= frame_shifts[j,mcnt,:]
                                ids_assigned.append([j,mcnt])
                                tempn+=1

                if cur > 19:

                    cc = 0
                    for x,y in coordinates_sorted[:,nr,:]:
                        if x > 0 and y>0: cc+=1
                    cur = cc
                    tot += cur

                    index = coordinates_sorted[:,nr,0].argmax()
                    mx,my = coordinates_sorted[index,nr,:]
                    printline +=  "\tindex {:4d}: {}/{}\n".format(int(cntr), cur, r-sum(excluded))
                    listdx.append('{:02d}/{}'.format(cur,r-sum(excluded)))
                    outc[:,imark,:] = coordinates_sorted[:,nr,:]
                    # for j,mcnt in ids_assigned: self.assignedFiducials[j][mcnt] = 1
                    try: add_marker(['Marker_{:03d}'.format(imark), '{:02d}/{:02d}'.format(cur,r-sum(excluded))])
                    except: pass
                    imark += 1

                elif cur>5:
                    if nr > 0.5 and (index_map == nr).sum() > 2:
                        index = coordinates_sorted[:,nr,0].argmax()
                        mx,my = coordinates_sorted[index,nr,:]
                        if (mx+my) > 0:
                            printline += "\tindex {:4d}: {}/{}\n".format(int(cntr), cur, r-sum(excluded))
                            listdx.append('{:02d}/{}'.format(cur,r-sum(excluded)))

                        outc[:,imark,:] = coordinates_sorted[:,nr,:]
                        # for j, mcnt in ids_assigned: self.assignedFiducials[j][mcnt] = 2
                        try: add_marker(['Marker_{:03d}'.format(imark), '{:02d}/{:02d}'.format(cur,r-sum(excluded))])
                        except: pass
                        imark += 1
                else:
                    coordinates_sorted[:,nr,:] = 0.*coordinates_sorted[:,nr,:] -1


        print ("indexed {}/{} potential fiducials.".format(tot, (mark_frames.sum(axis=-1) > 0).sum()) )
        print (printline)

        outc = outc[:,:imark,:]

        frame_shifts_sorted = np.zeros_like(frame_shifts)

        '''
        for j in range(frame_shifts.shape[0]):
            for i in range(frame_shifts.shape[1]):
                if np.abs(coordinates_sorted[j][i]).sum() > 0.1:
                    frame_shifts_sorted[j][i] = frame_shifts[j][0]
        '''
        if diag:
            print ('index frames: {:.0f} sec'.format((time()-s)) )
            s =time()
        mark_frames -= frame_shifts
        return outc, index_map, frame_shifts_sorted, listdx

    def find_potential_fiducials(self, frames, frames_full, raw, bin, bin_full, cropsize=50,target=0,fid_list=[],proc_id=0,num_procs=1,
                                 average_marker=None,threshold=1.7,mrcdata=[], sen=False, radius=10):
        from time import time

        s = time()

        list_cx_cy_imnr = []


        x,y = frames[0,:,:].shape
        sum2 = 0
        fact = (bin*1.)/bin_full
        ss = cropsize//2

        for imnr in range( len(frames) ):
            #if len(fid_list) and not imnr==target: continue
            st = time()
            image =  gaussian_filter(frames[imnr,:,:],1)
            #cf = np.zeros_like(frames_full[0])
            lap = laplace(image.copy())
            
            lap -= lap.min()
            lap /= lap.max()
            lap[lap< lap.mean()*threshold]=lap.min()
            label_im, nb_labels = label(lap)
            laps = remove_small_objects(label_im,3)
            label_im, nr_labels = label(laps)
            for n in range(1, nr_labels+1):
                cx,cy = center_of_mass(label_im==n)

                #refine center by checking its position on large frame
                ccx, ccy =  cx,cy#int(np.round(cx)),int(np.round(cy) )


                xmin,xmax = max(0,int(ccx*fact)-ss), int(ccx*fact)+ss
                ymin,ymax = max(0,int(ccy*fact)-ss), int(ccy*fact)+ss
                cutout2 = (frames_full[imnr,xmin:xmax,ymin:ymax]).copy()
                cutout = cutout2-cutout2.min()
                cutout /= cutout.max()
                cutout = gaussian_filter(cutout,2)
                l = laplace((cutout.max()-cutout)/cutout.max())

                gfl = gaussian_filter((l.max()-l ),5)
                gfl[gfl<np.median(gfl)*1.1] = np.median(gfl)*1.1
                gfl -= gfl.min()
                if gfl.max(): gfl /= gfl.max()
                local_maxi = peak_local_max(gfl, indices=False, footprint=np.ones((5, 5)) )
                lm = local_maxi.copy()
                #imshow(lm+cutout2*(lm==0))
                #show()
                tot = local_maxi.sum()
                if tot == 0:
                    sum2+=1
                    #if abs(cy-226) < 8: print imnr,cx,cy
                    add = True
                    for a,b,c in list_cx_cy_imnr:
                        if np.sqrt(abs(cx- a)**2+abs(cy-b)**2) < 4 and imnr*num_procs+proc_id ==c:
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
                            if np.sqrt(abs((rx)/fact - a)**2+abs((ry)/fact-b)**2) < 4 and imnr*num_procs+proc_id == c:
                                add = False
                                break
                        if add and cx+16 < y and cy+16 < x:
                            dcx,dcy = (rx)/max(fact,1),(ry)/max(fact,1)

                            #if abs(dcy-226) < 8:print imnr, dcx,dcy
                            list_cx_cy_imnr.append([float(dcx),float(dcy),imnr*num_procs+proc_id])
                            #frames_full[imnr][ int((cx+(ccx*fact)-ss))][int((cy+(ccy*fact)-ss))] = frames_full[imnr].max()

                        else:
                            lm[int(cx)][int(cy)] = 0
                        local_maxi[int(cx)][int(cy)] = False

        if sen==False: fid_list += list_cx_cy_imnr
        return np.zeros_like(self.frames_full[0]), list_cx_cy_imnr

    def find_potential_fiducials_sensitive(self, frames, frames_full, raw, bin, bin_full, cropsize=32,target=0,fid_list=[],proc_id=0,num_procs=1,
                                           average_marker=None, threshold=1.7, mrcdata=[], aa=20, radius=10):

        a, l1 = self.find_potential_fiducials(frames, frames_full, bin, bin_full, cropsize, target, fid_list, proc_id, num_procs, average_marker, sen=True)

        list_cx_cy_imnr = []
        fact = 1.*bin/bin_full
        for imnr in range(proc_id, len(frames_full), num_procs):

            ds = frames[imnr]  #downsample(wiener(full_image),n)**0.7
            ds2 = ds.copy()
            for xx in range(0,len(ds)):
                for yy in range(0,len(ds[0])):
                    dd = np.median(ds[max(0,xx-10):xx+10,max(0,yy-10):yy+10])


                    ds2[xx,yy] = ds[xx,yy]-dd
            ds = ds2

            minf = minimum_filter( gaussian_filter(ds - gaussian_filter(ds, 5),.1) , footprint=np.ones((1,1)) )

            hpf = minf - gaussian_filter(minf,4)
            hpf -= hpf.min()
            hpf /= np.median(hpf)
            hpf = 1-hpf

            gg = ((gaussian_filter(hpf,.1) )  )
            gg -= gg.min()
            gg /= np.median(gg)

            rr = gaussian_filter( laplace( gaussian_filter( (1- gg/gg.max()), .1) ) , 1)
            rr -= rr.min()
            rr /= rr.max()
            rr += 0.0001

            tot = np.zeros_like(rr)
            for ii in np.arange(4.1,3.59,-0.05):
                rrr = rr.copy()
                tot += self.testing_done(rrr,aa=aa,ff=ii)

                ll, nn = label(tot)


                for N in range(1,nn+1):
                    cx,cy = center_of_mass(ll==N)
                    #TOT[int(round(cx)),int(round(cy))] = 1
                    add = True
                    for a,b,c in list_cx_cy_imnr:
                        if abs(cx- a)+abs(cy-b) < 4 and imnr == c:
                            add = False
                            break
                    if add:
                        #dcx,dcy = (cx+(ccx*fact)-ss)/max(fact,1),(cy+(ccy*fact)-ss)/max(fact,1)
                        list_cx_cy_imnr.append([cx,cy,imnr])




        out = l1
        for rx,ry,im in list_cx_cy_imnr:
            #print rx,ry,imnr
            points = self.refine(rx, ry, frames_full[im], im, fact, cropsize)
            #print 'points: ', points
            for cx,cy,ims in points:
                add = True
                for a,b,c in out:
                    if abs(cx- a)+abs(cy-b) < 4 and ims == c:
                        add = False
                        break
                if add:
                    out.append([cx,cy,ims])

        fid_list += out #list_cx_cy_imnr




        return np.zeros_like(frames_full[0]), list_cx_cy_imnr

    def find_potential_fiducials_crosscorr(self, frames, frames_full, raw, bin, bin_full, cropsize=50, target=0, fid_list=[],
                                           proc_id=0, num_procs=1, average_marker=None, threshold=1.7, radius=10):
        list_cx_cy_imnr = []
        radius = ss = cropsize // 2
        sum2 = 0
        import time
        fact = bin
        ccmaps = []
        for imnr in range(len(self.fnames)):

            startTime = time.time()
            #image = frames_full[imnr]
            image = self.frames_adj[imnr*num_procs+proc_id,:,:]
            cc_map = self.cross_correlate(image, average_marker)
            ccmaps.append(cc_map)
            gol = gaussian_laplace(downsample(cc_map, fact), 1)
            gol -= gol.min()
            gol /= gol.max()
            gol = 1 - gol
            gol[gol < gol.mean() * threshold] = 0

            label_img, nb_labelsg = label(gol > 0)
            lapsg = remove_small_objects(label_img, 0)
            label_imgg, nr_labelsgg = label(lapsg)


            for n in range(1, nr_labelsgg + 1):
                cx, cy = center_of_mass(label_imgg == n)

                # refine center by checking its position on large frame
                ccx, ccy = int(round(cx)), int(round(cy))  # int(np.round(cx)),int(np.round(cy) )

                xmin, xmax = max(0, int(ccx * fact) - ss), int(ccx * fact) + ss
                ymin, ymax = max(0, int(ccy * fact) - ss), int(ccy * fact) + ss
                cutout2 = (image[xmin:xmax, ymin:ymax]).copy()
                cutout = cutout2 - cutout2.min()
                cutout /= cutout.max()
                cutout = gaussian_filter(cutout, 2)
                l = laplace((cutout.max() - cutout) / cutout.max())

                gfl = gaussian_filter((l.max() - l), 5)
                gfl[gfl < np.median(gfl) * 1.1] = np.median(gfl) * 1.1
                gfl -= gfl.min()
                if gfl.max(): gfl /= gfl.max()

                ly = cc_map[xmin:xmax, ymin:ymax].flatten().argmax()
                fx, fy = ly // cutout2.shape[1], ly % cutout2.shape[1]

                cutout2[fx - 2:fx + 3, fy - 2:fy + 3] = 2 * cutout2.max()

                tempy = cc_map[xmin:xmax, ymin:ymax].copy()
                temp = tempy - tempy.min()
                temp /= temp.max()
                temp = gaussian_filter(tempy, 3)
                temp -= temp.min()
                temp /= temp.max()

                # local_maxi = peak_local_max(temp, indices=False, threshold_rel=0.75 )
                local_maxi = peak_local_max(temp, indices=False, footprint=np.ones((5, 5)), threshold_rel=0.6,
                                            exclude_border=12)
                lm = local_maxi.copy()

                tot = local_maxi.sum()
                if tot > 6:
                    continue

                for i in range(tot):
                    x, y = cutout.shape

                    cx, cy = (local_maxi.argmax()) / y, (local_maxi.argmax()) % y
                    rx, ry = cx + ccx * fact - ss, cy + ccy * fact - ss

                    if cropsize > cutout.shape[1] and ccy + ss < cutout.shape[1]: ry += cropsize - cutout.shape[1]
                    if cropsize > cutout.shape[0] and ccx + ss < cutout.shape[0]: rx += cropsize - cutout.shape[0]

                    # If two peaks are close to each other, second peak is not considered.
                    add = True
                    for a, b, c in list_cx_cy_imnr:
                        if abs((rx) - a * fact) + abs((ry) - b * fact) < 5 and imnr  == c:
                            add = False
                            break
                    if add and cx + 16 < image.shape[1] and cy + 16 < image.shape[0]:
                        dcx, dcy = (rx) / max(fact, 1), (ry) / max(fact, 1)
                        rx, ry = int(rx), int(ry)
                        list_cx_cy_imnr.append([float(dcx), float(dcy), imnr*num_procs+proc_id ])

                    else:
                        lm[int(cx)][int(cy)] = 0

                    local_maxi[int(cx)][int(cy)] = False


        fid_list += list_cx_cy_imnr
        return ccmaps, list_cx_cy_imnr

    def find_potential_fiducial_stdv(self, frames, frames_full, raw, bin, bin_full, cropsize=50, target=0, fid_list=[],
                                           proc_id=0, num_procs=1, average_marker=None, threshold=1.7, radius=10):


        import matplotlib
        try:
            matplotlib.use('Qt5Agg')
        except:
            pass
        from pylab import imshow, show, subplots, savefig
        from pytom.agnostic.tools import create_sphere, create_circle
        from pytom.agnostic.io import read, write
        from pytom.voltools import transform
        from pytom.gpu.initialize import xp, device
        from pytom.agnostic.correlation import meanVolUnderMask, stdVolUnderMask
        from pytom.agnostic.transform import resize
        import sys
        import os
        from skimage import restoration
        from skimage.morphology import watershed, label
        import scipy, skimage
        from skimage.feature import peak_local_max

        refine_points = True
        image = np.zeros_like(self.frames_full[0, :, :])
        list_cx_cy_imnr = []

        binning = bin_full
        fact = 1.*bin/bin_full

        for imnr in range(len(frames_full)):
            ID = imnr * num_procs + proc_id

            data = raw[imnr,:,:]
            data[data > data.mean() + 5 * data.std()] = np.median(data)

            if binning != 1:
                data = resize(data, 1/binning).squeeze()


            #mask2 = create_sphere(data.shape, radius=data.shape[0] * 2 // 5)
            mask = create_circle(data.shape, radius=radius)

            volume = data  # * mask2

            meanV = meanVolUnderMask(volume, mask)
            stdV = stdVolUnderMask(volume, mask, meanV)



            r = 8


            stdV[:, :r] = 0
            stdV[:r, :] = 0
            stdV[-r:, :] = 0
            stdV[:, -r:] = 0
            stdV[stdV < stdV.mean() + stdV.std() * threshold] = 0
            stdV /= stdV.max() / 100

            footprint = np.ones((int(radius), int(radius)))
            image = restoration.denoise_tv_chambolle(stdV.squeeze(), weight=0.1)
            # distance = scipy.ndimage.distance_transform_edt(image)
            # local_maxi = peak_local_max(image, indices=False, footprint=footprint)
            #
            # markers = skimage.morphology.label(local_maxi)
            # labels_ws = watershed(-distance, markers, mask=image)
            labels_ws = label(image > 0)
            particles = np.unique(labels_ws)




            for particle in particles[particles > 0]:
                x1, y1 = scipy.ndimage.center_of_mass(labels_ws * (labels_ws == particle))

                if refine_points:
                    points = self.refineStdv(x1, y1, data, ID, 1, cropsize=30)
                    for cx, cy, ims in points:
                        add = True
                        for a, b, c in list_cx_cy_imnr:
                            if np.sqrt((cx - a)**2 + (cy - b)**2) < 4 and ims == c:
                                add = False
                                break
                        if add:
                            list_cx_cy_imnr.append([cx/fact, cy/fact, ims])
                else:
                    list_cx_cy_imnr.append([float(x1)*bin_full/bin, float(y1)*bin_full/bin, ID])

        fid_list += list_cx_cy_imnr
        return image, list_cx_cy_imnr

    def find_potential_fiducial_log(self, frames, frames_full, raw, bin, bin_full, cropsize=50, target=0, fid_list=[],
                                           proc_id=0, num_procs=1, average_marker=None, threshold=1.7, radius=10):

        from pytom.agnostic.transform import resize
        from skimage.morphology import watershed, label
        import scipy, skimage
        from skimage.feature import peak_local_max

        binning = bin_full
        fact = 1.*bin/bin_full
        list_cx_cy_imnr = []

        for imnr in range(len(frames_full)):
            ID = imnr * num_procs + proc_id

            data = raw[imnr,:,:]
            data[data > data.mean() + 5 * data.std()] = np.median(data)

            if binning != 1:
                data = resize(data, 1/binning).squeeze()

            vol = scipy.ndimage.gaussian_laplace(data, radius)
            vol[vol < vol.mean()+vol.std()*threshold] = 0

            labels_log = label((vol>0))
            particles = np.unique(labels_log)




            for particle in particles[particles > 0]:
                x1, y1 = scipy.ndimage.center_of_mass(labels_log * (labels_log == particle))

                list_cx_cy_imnr.append([float(x1)*bin_full/bin, float(y1)*bin_full/bin, ID])

        fid_list += list_cx_cy_imnr
        return data, list_cx_cy_imnr


    def refineStdv(self, cx, cy, frame, imnr, fact, cropsize):
        '''This function finds the maxima in a 2D frame'''
        
        if len(frame.shape) != 2: raise Exception('wrong shape of input frame')
        try:
            float(cx)
            float(cy)
        except:
            raise Exception('Please provide float or int as input coordinates')
        list_cx_cy = []

        ss = cropsize // 2

        ccx, ccy = cx, cy  # int(numpy.round(cx)),int(numpy.round(cy) )

        cutout2 = (frame[max(0, int(ccx * fact) - ss):int(ccx * fact) + ss,
                   max(0, int(ccy * fact) - ss):int(ccy * fact) + ss]).copy()
        cutout = cutout2 - cutout2.min()
        cutout /= cutout.max()

        cutout = gaussian_filter(cutout, 2)
        l = laplace((cutout.max() - cutout) / cutout.max())

        gfl = gaussian_filter((l.max() - l), 5)
        gfl[gfl < np.median(gfl) * 1.1] = np.median(gfl) * 1.1
        gfl -= gfl.min()
        if gfl.max(): gfl /= gfl.max()

        local_maxi = peak_local_max(gfl, indices=False, footprint=np.ones((5, 5)))
        lm = local_maxi.copy()

        tot = local_maxi.sum()

        if tot == 0:
            if gfl.sum():
                cx, cy = np.unravel_index(gfl.argmax(), gfl.shape)
                list_cx_cy.append([cx/fact, cy/fact, imnr])
            else:
                list_cx_cy.append([cx/fact, cy/fact, imnr])

        else:

            for i in range(tot):
                x, y = cutout.shape

                cx, cy = (local_maxi.argmax()) / y, (local_maxi.argmax()) % y
                rx, ry = cx + ccx * fact - ss, cy + ccy * fact - ss
                if cropsize > cutout.shape[1] and ccy + ss < cutout.shape[1]: ry += cropsize - cutout.shape[1]
                if cropsize > cutout.shape[0] and ccx + ss < cutout.shape[0]: rx += cropsize - cutout.shape[0]

                # sometimes the algorithm gives two adjacent peaks. We do not keep the second.
                add = True
                for a, b, c in list_cx_cy:
                    if abs((rx) / fact - a) + abs((ry) / fact - b) < 4:
                        add = False
                        break
                if add:
                    dcx, dcy = (rx) / max(fact, 1), (ry) / max(fact, 1)
                    list_cx_cy.append([dcx, dcy, imnr])

        return list_cx_cy

    def refine(self, cx, cy, frame, imnr, fact, cropsize):
        '''This function finds the maxima in a 2D frame'''
        import matplotlib
        try:
            matplotlib.use('Qt5Agg')
        except:
            pass
        from pylab import subplots, imshow, show

        if len(frame.shape) != 2: raise Exception('wrong shape of input frame')
        try:
            float(cx)
            float(cy)
        except:
            raise Exception('Please provide float or int as input coordinates')
        list_cx_cy = []

        ss = cropsize//2

        ccx, ccy =  cx,cy#int(np.round(cx)),int(np.round(cy) )

        cutout2 = (frame[max(0,int(ccx*fact)-ss):int(ccx*fact)+ss,max(0,int(ccy*fact)-ss):int(ccy*fact)+ss]).copy()
        cutout = cutout2-cutout2.min()
        cutout /= cutout.max()
        cutout = gaussian_filter(cutout,2)
        l = laplace((cutout.max()-cutout)/cutout.max())

        gfl = gaussian_filter((l.max()-l ),5)
        gfl[gfl<np.median(gfl)*1.1] = np.median(gfl)*1.1
        gfl -= gfl.min()
        if gfl.max(): gfl /= gfl.max()
        local_maxi = peak_local_max(gfl, indices=False, footprint=np.ones((5, 5)) )
        lm = local_maxi.copy()
        # fig,ax = subplots(1,1)
        # ax.imshow(lm+cutout2*(lm==0))
        # show()
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


    def testing_done(self, rr,aa=20,ff=4):

        d1 = rr.copy()
        d2 = rr.copy()

        edge = np.zeros_like(rr)

        for xx in range(0,len(rr),aa):
            for yy in range(0,len(rr[0]),aa):

                dd = rr[xx:xx+aa,yy:yy+aa].copy()
                if dd.std() >0.7:
                    #print 'add', std(dd), dd.min(), dd.max()
                    edge[xx:xx+aa,yy:yy+aa] = 1
                #else:
                    #print 'tast', std(dd), dd.min(), dd.max()

                dd[dd < np.median(dd)+np.std(dd)*ff] = 0
                d1[xx:xx+aa,yy:yy+aa] = dd

        ll1,nn1 = label(remove_small_objects(d1>0,2))

        for xx in range(-aa//2,len(rr),aa):
            for yy in range(-aa//2,len(rr[0]),aa):

                dd = rr[max(0,xx):xx+aa,max(0,yy):yy+aa]
                if dd.min() > 0.2:
                    #print dd.min()
                    edge[xx:xx+aa,yy:yy+aa] = 1
                    #print median(dd),median(dd)+std(dd)*2, median(dd)*1.2
                dd[dd < np.median(dd)+np.std(dd)*ff] = 0
                d2[max(0,xx):xx+aa,max(0,yy):yy+aa] = dd

        ll2,nn2 = label(remove_small_objects(d2>0,2))



        return ((ll2+ll1) >0)


    def read_data(self, fnames,bin,bin_full):

        frames = []
        frames_full = []
        for fname in fnames:
                #data  = read('{}'.format(str(fname)),binning=[self.binning.get(),self.binning.get(),1])
            dataf  = read('{}'.format(str(fname)),binning=[bin_full,bin_full,1])

            w = wiener(dataf)
            hpf = dataf - gaussian_filter(dataf,10)*0
            data = downsample( w+hpf, int(round(bin*1./bin_full)))
                #data = vol2npy(f).copy()
                #data = read_mrc(fname)
            data -= data.min()
            data /= data.max()

            dataf -= dataf.min()
            dataf /= dataf.max()

            frames.append( (data)**0.75 )
            frames_full.append( wiener(dataf) )


        return frames, frames_full


    def calculate_error(self, coordinates, tilt_angles, projIndices, imdimX, imdimY, fname='alignmentErrors.txt', reference_marker=1):
        from pytom.reconstruction.tiltAlignmentFunctions import alignmentFixMagRot
        from pytom.reconstruction.TiltAlignmentStructures import Marker
        from pytom.gui.guiFunctions import ALIGNMENT_ERRORS, fmtAE as fmt, headerAlignmentErrors

        Markers = []
        for markerID in range(coordinates.shape[1]):
            coords = coordinates[:,markerID,:]
            marker = Marker(projIndices)
            num_hits = 0
            for n, (cy,cx) in enumerate(coords[projIndices]):

                if cx > 0.005 and cy > 0.005:
                    marker.set_xProj(n, cx)
                    marker.set_yProj(n, cy)
                    num_hits += 1
                else:
                    marker.set_xProj(n, -1)
                    marker.set_yProj(n, -1)

            if num_hits: Markers.append(marker)



        cTilt = np.cos(tilt_angles[projIndices]/ 180. * np.pi)
        sTilt = np.sin(tilt_angles[projIndices]/ 180. * np.pi)
        ireftilt = np.abs(tilt_angles[projIndices]).argmin()


        irefmark =  reference_marker if reference_marker < len(Markers) else 0

        psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY, shiftVarX, shiftVarY =\
                alignmentFixMagRot(Markers, cTilt, sTilt, ireftilt, irefmark=irefmark, imdimX=imdimX, imdimY=imdimY)


        psi = [psiindeg,]*len(Markers)

        nangles = len(projIndices)
        nmarkers = len(Markers)

        data = np.zeros((nmarkers*nangles),dtype=ALIGNMENT_ERRORS)

        errors = np.zeros(nmarkers)
        excluded = [0,]*nmarkers
        for imark in range(nmarkers):
            dx, dy = diffX[:,imark]-shiftX[:], diffY[:,imark]-shiftY[:]



            dist = np.sqrt(dx**2 + dy**2)

            dist[Markers[imark].xProj < -0.001] = -1
            dist[Markers[imark].yProj < -0.001] = -1

            errors[imark] = dist[dist>-0.5].mean()


            print(f'Mean residual marker_{imark:03d}: {dist[dist>-0.5].mean()}')

            data['AlignmentError'][imark*nangles: (imark+1)*nangles] = dist

            data['TiltAngle'][     imark*nangles: (imark+1)*nangles] = tilt_angles[projIndices]
            data['MarkerIndex'][   imark*nangles: (imark+1)*nangles] = imark

        if fname: np.savetxt(fname, data, fmt=fmt, header=headerAlignmentErrors)


        for threshold in (100, 50, 40, 30, 20, 20, 15, 10, 7, 5, 3.5):
            nm, excluded, irefmark = self.new_markerset(Markers, errors, threshold, excluded, irefmark)
            if len(nm) < 3: break
            tpsi, terrors, shiftX, shiftY, tdiffX, tdiffY, tx, ty, tz = self.calcscore(nm, cTilt, sTilt, ireftilt, imdimX, imdimY, projIndices, irefmark=irefmark)
            m = 0
            for n, excl in enumerate(excluded):
                if excl == 0:
                    errors[n] = terrors[m]
                    diffX[:,n] = tdiffX[:,m]
                    diffY[:,n] = tdiffY[:,m]
                    x[n] = tx[m]
                    y[n] = ty[m]
                    z[n] = tz[m]
                    psi[n] = tpsi
                    m += 1

            if nm:
                m=0

                for n, imark in enumerate(list(range(nmarkers))):
                    if excluded[n]: continue
                    dx, dy = diffX[:,n]-shiftX[:], diffY[:,n]-shiftY[:]

                    dist = np.sqrt(dx**2 + dy**2)

                    dist[Markers[imark].xProj < -0.001] = -1
                    dist[Markers[imark].yProj < -0.001] = -1

                    errors[imark] = dist[dist>-0.5].mean()


                    #print(f'Mean residual marker_{imark:03d}: {dist[dist>-0.5].mean()}')

                    data['AlignmentError'][imark*nangles: (imark+1)*nangles] = dist

                    data['TiltAngle'][     imark*nangles: (imark+1)*nangles] = tilt_angles[projIndices]
                    data['MarkerIndex'][   imark*nangles: (imark+1)*nangles] = imark
                    m += 1

                if fname: np.savetxt(fname, data, fmt=fmt, header=headerAlignmentErrors)


        self.errors_coordinates = errors



        #centers = self.redoCenters(coordinates, tilt_angles, imdimX, imdimY, max_angle=10)

        # for ang in np.arange(10, abs(tilt_angles).max(), 10):
        #     centers = self.redoCenters(coordinates, tilt_angles, imdimX, imdimY, max_angle=ang)
        #
        #     shift = list(zip(shiftX, shiftY))
        #     Markers = self.remove_outliers(Markers, centers, tilt_angles, shift, psiindeg, excluded, ireftilt)
        #     # centers = self.redoCenters(coordinates, tilt_angles, imdimX, imdimY, max_angle=90)
        #
        #     Markers = self.find_coordinates(Markers, centers, tilt_angles, shift, ireftilt, psiindeg, self.coordinates,
        #                                     self.mark_frames * self.bin_read, self.dataRaw, excluded)
        #     break
        #
        # psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY, shiftVarX, shiftVarY =\
        #         alignmentFixMagRot(Markers, cTilt, sTilt, ireftilt, irefmark=1,imdimX=imdimX, imdimY=imdimY)
        #
        #
        # psi = [psiindeg,]*len(Markers)
        #
        # nangles = len(projIndices)
        # nmarkers = len(Markers)
        #
        # data = np.zeros((nmarkers*nangles),dtype=ALIGNMENT_ERRORS)
        #
        # errors = np.zeros(nmarkers)
        # excluded = [0,]*nmarkers
        # for imark in range(nmarkers):
        #     dx, dy = diffX[:,imark]-shiftX[:], diffY[:,imark]-shiftY[:]
        #
        #
        #
        #     dist = np.sqrt(dx**2 + dy**2)
        #
        #     dist[Markers[imark].xProj < -0.001] = -1
        #     dist[Markers[imark].yProj < -0.001] = -1
        #
        #     errors[imark] = dist[dist>-0.5].mean()
        #
        #
        #     print(f'Mean residual marker_{imark:03d}: {dist[dist>-0.5].mean()}')
        #
        #     data['AlignmentError'][imark*nangles: (imark+1)*nangles] = dist
        #
        #     data['TiltAngle'][     imark*nangles: (imark+1)*nangles] = tilt_angles[projIndices]
        #     data['MarkerIndex'][   imark*nangles: (imark+1)*nangles] = imark
        #
        # if fname: np.savetxt(fname, data, fmt=fmt, header=headerAlignmentErrors)


        return errors, shiftX, shiftY, diffX, diffY, x, y, z, psi

    def find_coordinates(self, markers, centers, tiltangles, shifts, ireftilt, psiindeg, coordinates, foundFid, rawData, excluded):

        for imark, marker in enumerate(markers):
            for itilt, tiltangleindeg in enumerate(tiltangles[(1-excluded) > 0.5]):
                x,y = marker.get_xProj(itilt), marker.get_yProj(itilt)
                if x < 0 or y< 0:
                    v0 = projectMarkerToFrame(centers[imark], tiltangleindeg, psiindeg, shifts[ireftilt], shifts[itilt])
                    for v1 in foundFid[itilt,:,:]:
                        if self.distance2(v0, v1, psiindeg, 4) < 4 and self.diff(v0,v1) < 10:
                            markers[imark].set_xProj(v0[0])
                            markers[imark].set_yProj(v0[1])
                            self.coordinates[itilt,imark,:] = v0
                            # self.assignedFiducials[itilt][imark] = 1

        return markers

    def remove_outliers(self, Markers, centers, tilt_angles, shifts, psiindeg, excluded, ireftilt, max_shift=20):

        for imark, center in enumerate(centers):
            for itilt, tiltangleindeg in enumerate(tilt_angles):
                cx, cy = projectMarkerToFrame(center, tiltangleindeg, psiindeg, shifts[ireftilt], shifts[itilt])
                mx, my = Markers[imark].get_xProj(), Markers[imark].get_yProj()
                dist = (mx-cy)**2 + (my-cy)**2
                if dist > max_shift**2:
                    Markers[imark].set_xProj(itilt, -1)
                    Markers[imark].set_yProj(itilt, -1)

    def redoCenters(self, coordinates, tilt_angles, imdimX, imdimY, max_angle=10):
        from pytom.reconstruction.tiltAlignmentFunctions import alignmentFixMagRot
        from pytom.reconstruction.TiltAlignmentStructures import Marker
        from pytom.gui.guiFunctions import ALIGNMENT_ERRORS, fmtAE as fmt, headerAlignmentErrors

        projIndices = np.arange(0, len(coordinates), 1).astype(int)[abs(tilt_angles) <= max_angle + 0.1]
        Markers = []

        for markerID in range(coordinates.shape[1]):
            coords = coordinates[:,markerID,:]
            marker = Marker(projIndices)
            num_hits = 0
            for n, (cy,cx) in enumerate(coords[projIndices]):

                if cx > 0.005 and cy > 0.005:
                    marker.set_xProj(n, cx)
                    marker.set_yProj(n, cy)
                    num_hits += 1
                else:
                    marker.set_xProj(n, -1)
                    marker.set_yProj(n, -1)

            if num_hits: Markers.append(marker)



        cTilt = np.cos(tilt_angles[projIndices]/ 180. * np.pi)
        sTilt = np.sin(tilt_angles[projIndices]/ 180. * np.pi)
        ireftilt = np.abs(tilt_angles[projIndices]).argmin()


        psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY, shiftVarX, shiftVarY =\
                alignmentFixMagRot(Markers, cTilt, sTilt, ireftilt, irefmark=1,imdimX=imdimX, imdimY=imdimY)
        centers = list(zip(x, y, z))

        return centers

    def new_markerset(self, markers, errors, threshold, excluded, irefmark):
        new_markers = []
        for n, ex in enumerate(excluded):
            if ex or errors[n] > threshold:
                excluded[n] = 1
            else:
                new_markers.append(markers[n])

        irefmark = 0 if excluded[irefmark] else irefmark

        return new_markers, excluded, irefmark

    def calcscore(self, Markers, cTilt, sTilt, ireftilt, imdimX, imdimY, projIndices, irefmark):
        from pytom.reconstruction.tiltAlignmentFunctions import alignmentFixMagRot

        irefmark = irefmark if irefmark < len(Markers) else 0

        psiindeg, shiftX, shiftY, x, y, z, distLine, diffX, diffY, shiftVarX, shiftVarY =\
                alignmentFixMagRot(Markers, cTilt, sTilt, ireftilt, irefmark=irefmark,imdimX=imdimX, imdimY=imdimY)

        nangles = len(projIndices)
        nmarkers = len(Markers)

        #data = np.zeros((nmarkers*nangles),dtype=ALIGNMENT_ERRORS)

        errors = np.zeros(nmarkers)
        #excluded = [0,]*nmarkers
        for imark in range(nmarkers):
            dx, dy = diffX[:,imark]-shiftX[:], diffY[:,imark]-shiftY[:]



            dist = np.sqrt(dx**2 + dy**2)

            dist[Markers[imark].xProj < -0.001] = -1
            dist[Markers[imark].yProj < -0.001] = -1

            errors[imark] = dist[dist>-0.5].mean()

        return [psiindeg, errors, shiftX, shiftY, diffX, diffY, x,y,z]

    def determine_markerdata(self, coordinates, errors, excluded, add_marker):
        incl = 1 - np.array(excluded, dtype=int)

        for markerID in range(coordinates.shape[1]):
            coords = coordinates[:,markerID,:]
            num_hits = 0
            for n, (cx,cy) in enumerate(coords):
                if incl[n] and cx > 0.005 and cy > 0.005:
                    num_hits += 1

            if num_hits:
                add_marker([f'Marker_{markerID:03d}', f'{num_hits:02d}/{incl.sum():02d}', f'{errors[markerID]:.2f}'])


def projectMarkerToFrame(center, tiltangleindeg, psiindeg, shift0, shift1):

    cpsi, spsi = np.cos(psiindeg * np.pi / 180), np.sin(psiindeg * np.pi / 180)
    ctlt, stlt = np.cos(tiltangleindeg * np.pi / 180), np.sin(tiltangleindeg * np.pi / 180)
    cx, cy, cz = center

    cx -= shift0[0]
    cy -= shift0[1]

    cy1 = cy * cpsi - cx * spsi
    cx1 = cy * spsi + cx * cpsi

    cy2 = cy1 * ctlt  - cz * stlt

    cx = cx1 * cpsi - cy2 * spsi
    cy = cx1 * spsi + cy2 * cpsi

    cx += shift1[0]
    cy += shift1[1]

    return np.array((cx, cy))
