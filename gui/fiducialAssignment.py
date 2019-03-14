import sys
import os
import glob
import numpy as np
import copy

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtWidgets, QtCore, QtGui

from pytom.gui.guiStyleSheets import *
from pytom.gui.mrcOperations import *

import pyqtgraph as pg
from pyqtgraph.GraphicsScene.mouseEvents import MouseClickEvent
from pyqtgraph import ImageItem

from scipy.ndimage.filters import gaussian_filter
from pytom.gui.fiducialPicking import *

global mf_write

try:
    from pytom.gui.additional.markerPositionRefinement import refineMarkerPositions
    from pytom.gui.additional.tiltAlignmentFunctions import alignmentFixMagRot
    from pytom.gui.additional.TiltAlignmentStructures import Marker
    from pytom.basic.files import write_em
    from pytom_volume import vol, read, transform
    from pytom_numpy import vol2npy
    mf_write=1
except:
    print ('marker file refinement is not possible')
    mf_write = 0

import sys
import os
import glob
import numpy
import time
import shutil


from scipy.ndimage import sobel, gaussian_filter, laplace, label, median_filter, maximum_filter, convolve
from scipy.ndimage import binary_closing, binary_fill_holes, binary_erosion, center_of_mass
from scipy.ndimage.filters import minimum_filter
from scipy.ndimage.morphology import generate_binary_structure
from scipy.signal import wiener, argrelextrema

from skimage.morphology import remove_small_objects
from skimage.feature import canny
from skimage.morphology import watershed
from skimage.feature import peak_local_max


from multiprocessing import Process, Event, Manager, Pool, cpu_count

from pytom.gui.guiStructures import *

def sort_str( obj, nrcol ):
    obj.sort(key=lambda i: str(i[nrcol]))

class TiltImage():
    def __init__(self, parent, imnr, excluded=0):

        self.frame_full = parent.frames_full[imnr]
        self.frame = parent.frames[imnr]
        self.fiducials = []
        self.indexed_fiducials = []
        self.labels = []
        self.imnr = imnr

        self.tiltangle = parent.tiltangles[imnr]
        self.excluded = excluded
        self.main_image = parent.main_image
        self.sel_image = parent.bottom_image
        self.sel_circles = parent.bottom_circles
        self.main_circles = parent.main_circles
        self.parent = parent

    def clear(self):
        self.fiducials = []
        self.indexed_fiducials = []
        self.labels = []

    def add_fiducial(self, x, y, FX, FY, cutoff=5, check=True, draw=True, label=''):
        added = False
        new = True
        radius = self.parent.radius

        if check:
            for n, (fx, fy, fi) in enumerate(self.fiducials):
                if fx >0 and fy >0 and abs(x-fx) < cutoff and abs(y-fy)< cutoff:
                    new = False
                    return new


        self.fiducials.append([FX, FY, self.imnr])
        self.indexed_fiducials.append( Qt.red )
        self.labels.append( label )

        if draw:
            self.sel_circles.append(circle(QPoint(x - radius, y - radius), size=radius * 2, color=Qt.blue))
            self.sel_image.addItem(self.sel_circles[-1])

            self.main_circles.append(circle(QPoint(FX - radius, FY - radius),size=radius * 2,color=self.indexed_fiducials[-1]))
            self.main_image.addItem(self.main_circles[-1])

        return new

    def remove_fiducial(self, x, y, cutoff=5):
        removed = 0
        for n, (fx, fy, fi) in enumerate(self.fiducials):

            if abs(fx - x) < cutoff and abs(fy - y) < cutoff:
                self.fiducials = self.fiducials[:n-removed] + self.fiducials[n + 1 - removed:]
                self.indexed_fiducials = self.indexed_fiducials[:n - removed] + self.indexed_fiducials[n + 1 - removed:]
                self.labels = self.labels[:n - removed] + self.labels[n + 1 - removed:]
                removed += 1

    def update_indexing(self, points, cutoff=0.2):
        self.indexed_fiducials = [Qt.red,]*len(self.fiducials)
        self.labels = ['',]*len(self.fiducials)
        for m, (py,px) in enumerate(points):
            for n, (fx,fy,imnr) in enumerate(self.fiducials):
                if abs(fx-px) < cutoff and abs(fy-py) < cutoff:
                    self.indexed_fiducials[n] = Qt.green
                    self.labels[n] = "{:03d}".format(m)
                    break


class FiducialAssignment(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, fname=''):
        super(FiducialAssignment, self).__init__(parent)
        self.size_policies()
        if 1:
            self.pytompath = self.parent().pytompath
            self.projectname = self.parent().projectname
        #except:
            #self.projectname = '/Users/gijs/Documents/PostDocUtrecht/Data/Juliette'
            #if not os.path.exists(self.projectname): self.projectname = '/home/gijsvds/testcase'
        #    pass

        self.tomofolder = os.path.join(self.projectname, '03_Tomographic_Reconstruction')
        self.tomogram_names = sorted( [line.split()[0] for line in os.listdir(self.tomofolder) if line.startswith('tomogram_')] )
        self.loaded_data = False
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 900, 600)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.settings = SettingsFiducialAssignment(self)
        self.selectMarkers = SelectAndSaveMarkers(self)
        self.manual_adjust_marker = ManuallyAdjustMarkers(self)

        self.dim = 0
        self.radius = 8
        self.sizeCut = 200
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.xmin, self.ymin = 0,0
        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.main_canvas = w0 = pg.GraphicsWindow(size=(600, 600), border=True)
        self.main_image = w0.addViewBox(row=0, col=0, lockAspect=True)
        self.main_image.setMouseEnabled(False)
        self.main_image.invertY(True)


        #self.top_canvas = w1 = KeyPressGraphicsWindow(size=(300, 300), border=True)
        #self.top_image = w1.addViewBox(row=0, col=0, lockAspect=True)
        #self.top_image.setMenuEnabled(False)
        self.actionwidget = QWidget(self,width=300,height=300)
        self.al = parent = QGridLayout(self)
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.widgets = {}
        self.insert_label(parent, rstep=1)
        self.insert_pushbutton(parent,text='Find Fiducials',rstep=1,tooltip='Automatically detect fiducials.',
                               action=self.find_fid,params=0, wname='findButton', state=False)

        self.insert_pushbutton(parent,text='Detect Frame Shifts',rstep=1, wname='detectButton',state=False,
                               tooltip='Detect global x-y shifts between tilt images.',
                               action=self.detect_frameshift,params=0)
        self.insert_pushbutton(parent,text='Index Fiducials',rstep=1,tooltip='Group Fiducials into marker sets.',
                               action=self.index_fid,params=0, wname='indexButton',state=False)
        self.insert_label(parent, rstep=1)
        self.insert_pushbutton(parent,text='Manually Adjust Markers',rstep=1,tooltip='Manually adjust marker sets.',
                               action=self.raise_window, params=self.manual_adjust_marker)
        self.insert_pushbutton(parent, text='Create Markerfile', rstep=1, tooltip='Select and save marker sets.',
                               action=self.raise_window, params=self.selectMarkers)
        self.insert_label(parent, rstep=1)
        self.insert_pushbutton(parent,text='Settings',rstep=1,tooltip='Adjust Settings.',
                               action=self.raise_window,params=self.settings)
        self.insert_label(parent,rstep=1,sizepolicy=self.sizePolicyA)
        self.insert_label(parent, text='MANUALLY PICK FIDUCIALS', rstep=1,alignment=Qt.AlignHCenter)


        self.actionwidget.setLayout(self.al)

        self.bottom_canvas = w2 = pg.GraphicsWindow(size=(300, 300), border=True)
        self.bottom_image = w2.addViewBox(row=0, col=0, lockAspect=True)
        #self.bottom_image.setMouseEnabled(False, False)
        self.bottom_image.setMenuEnabled(False)
        self.bottom_image.invertY(True)

        self.main_canvas.wheelEvent = self.wheelEvent


        self.image_list = [self.main_image, self.bottom_image]
        self.main_image.scene().sigMouseClicked.connect(self.selectPartTiltImage)
        self.bottom_image.scene().sigMouseClicked.connect(self.manuallyAdjustFiducials)

        self.img1a = pg.ImageItem(zeros((10,10)))
        #self.img1b = pg.ImageItem(zeros((10,10)))
        self.img1c = pg.ImageItem(zeros((10,10)))

        self.main_image.addItem(self.img1a)
        #self.top_image.addItem(self.img1b)
        self.bottom_image.addItem(self.img1c)

        self.layout.addWidget(w0, 0, 0, 2, 1)
        self.layout.addWidget(self.actionwidget, 0, 1)
        self.layout.addWidget(w2, 1, 1)
        #self.layout.addWidget(self.operationbox, 1, 0)
        self.title = 'Assign fiducials' #parent.widgets['v03_manualpp_tomogramFname'].text()
        self.setWindowTitle(self.title)


        self.bottom_circles = []
        self.main_circles = []
        self.raise_window(self.settings)
        pg.QtGui.QApplication.processEvents()
        self.loaded_data = False

    def raise_window(self,window):
        window.close()
        window.show()

    def selectPartTiltImage(self,event):

        pos = self.main_image.mapSceneToView( event.scenePos() )

        if pos.x() > 0 and pos.y() >0 and pos.x() < self.dim and pos.y() < self.dim:
            for child in self.bottom_image.allChildren():
                if type(child) == MyCircleOverlay:
                    self.bottom_image.removeItem(child)

            xmin,ymin = int(round(max(0,pos.x()-self.sizeCut//2))), int(round(max(0,pos.y()-self.sizeCut//2)))

            if xmin + self.sizeCut >= self.dim: xmin = self.dim - self.sizeCut
            if ymin + self.sizeCut >= self.dim: ymin = self.dim - self.sizeCut

            self.xmin, self.ymin = xmin, ymin

            self.img1c.setImage(image=self.frames_full[self.imnr][xmin:xmin+self.sizeCut, ymin:ymin+self.sizeCut])
            self.replot()
        #self.top_image.scen().mouseHasMoved()

        pg.QtGui.QApplication.processEvents()

    def manuallyAdjustFiducials(self,event):
        if not self.dim: return
        pos = self.bottom_image.mapSceneToView(event.scenePos())
        x, y = int(round(pos.x())), int(round(pos.y()))

        if x>0 and y>0 and y<self.sizeCut and x < self.sizeCut:
            fx, fy = x + self.xmin, y + self.ymin

            if event.button() == 1:
                self.tiltimages[self.imnr].add_fiducial(x,y,fx,fy,cutoff=self.radius)

            elif event.button() == 2:
                self.tiltimages[self.imnr].remove_fiducial(fx, fy, cutoff=self.radius)

                #THE CHILDREN DO NOT HAVE REPRODUCIBLE COORDINATES
                self.remove_all_circles()
                self.add_circles_new_frame()

    def wheelEvent(self, event):

        step = event.angleDelta().y() // 120

        if self.imnr + step < len(self.fnames) and self.imnr + step > -1:
            self.imnr += step
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_F == evt.key():
            print ('Find Fiducials')
            self.find_fid()

        elif Qt.Key_D == evt.key():
            print ('Detect Frame Offsets')
            self.detect_frameshift()

        elif Qt.Key_I == evt.key():
            print ('Index Fiducials')
            self.index_fid()

        elif Qt.Key_E == evt.key():
            self.exclude_status_change()

        elif Qt.Key_R == evt.key():
            print ('Refine Marker Positions')

        elif Qt.Key_Right == evt.key():
            if self.imnr + 1 < len(self.fnames):
                self.imnr += 1
                self.replot2()

        elif Qt.Key_Left == evt.key():
            if self.imnr - 1 >= 0:

                self.imnr -= 1
                self.replot2()

        elif evt.key() == Qt.Key_Escape:
            self.close()

    def exclude_status_change(self):
        i = self.imnr
        f = self.fnames[i].replace('/excluded','')
        fn,dn = os.path.basename(f),os.path.dirname(f)

        self.excluded[i] = 1-self.excluded[i]
        if self.excluded[i]:
            os.system('mv {} {}/excluded/'.format(f,dn) )
            print ('Excluded Tilt Image {}'.format(fn))
        else:
            os.system('mv {}/excluded/{} {}'.format(dn,fn,dn) )
            print ('Included Tilt Image {}'.format(fn))

        self.replot()

    def replot(self):
        self.img1a.setImage(image=self.frames_full[self.imnr])
        #self.img1b.setImage(image=self.frames_full[self.imnr])
        xmax,ymax = self.xmin+self.sizeCut,self.ymin+self.sizeCut
        self.img1c.setImage(image=self.frames_full[self.imnr][self.xmin:xmax,self.ymin:ymax])
        self.remove_all_circles()
        self.add_circles_new_frame2()
        self.main_image.setRange(xRange=(0,self.dim),yRange=(0,self.dim),padding=None)
        title = '{} (tilt = {:5.1f} deg){}'.format(self.fnames[self.imnr].split('/')[-1],
                                                           self.tiltangles[self.imnr],
                                                           ' excluded'*self.excluded[self.imnr])
        self.setWindowTitle(title)

    def add_circles_new_frame(self):
        for n,(fx,fy,imnr) in enumerate(self.tiltimages[self.imnr].fiducials):
            color = self.tiltimages[self.imnr].indexed_fiducials[n]
            self.main_circles.append( circle(QPoint(fx - self.radius, fy - self.radius),
                                               size=self.radius*2, color=color))
            self.main_image.addItem(self.main_circles[-1])

            if abs(self.xmin+self.sizeCut//2 - fx) < self.sizeCut//2 and \
                abs(self.ymin+self.sizeCut//2 - fy) < self.sizeCut//2:
                self.bottom_circles.append( circle(QPoint(fx - self.xmin - self.radius,
                                                        fy - self.ymin - self.radius),
                                                 size=self.radius*2, color=Qt.blue))
                self.bottom_image.addItem(self.bottom_circles[-1])

    def replot2(self):
        self.img1a.setImage(image=self.frames_full[self.imnr])
        #self.img1b.setImage(image=self.frames_full[self.imnr])
        xmax,ymax = self.xmin+self.sizeCut,self.ymin+self.sizeCut
        self.img1c.setImage(image=self.frames_full[self.imnr][self.xmin:xmax,self.ymin:ymax])
        self.remove_all_circles()
        self.add_circles_new_frame2()
        self.main_image.setRange(xRange=(0,self.dim),yRange=(0,self.dim),padding=None)
        title = '{} (tilt = {:5.1f} deg){}'.format(self.fnames[self.imnr].split('/')[-1], self.tiltangles[self.imnr],
                                                   ' excluded'*self.excluded[self.imnr])
        self.setWindowTitle(title)

    def add_circles_new_frame2(self):

#        for n,(fx,fy) in enumerate(self.mark_frames[self.imnr]):
        for n,(fx,fy, dummy) in enumerate(self.tiltimages[self.imnr].fiducials):
            if fx < 0.1 or fy<0.1: continue

            fc = self.bin_alg/self.bin_read
            color = self.tiltimages[self.imnr].indexed_fiducials[n]
            label = self.tiltimages[self.imnr].labels[n]
            self.main_circles.append( circle(QPoint(  fx - self.radius, fy - self.radius),
                                               size=self.radius*2, color=color, label=label))
            self.main_image.addItem(self.main_circles[-1])

            if abs(self.xmin+self.sizeCut//2 - fx) < self.sizeCut//2 and \
                abs(self.ymin+self.sizeCut//2 - fy) < self.sizeCut//2:
                self.bottom_circles.append( circle(QPoint(fx - self.xmin - self.radius, fy - self.ymin - self.radius),size=self.radius*2, color=Qt.blue))
                self.bottom_image.addItem(self.bottom_circles[-1])

    def remove_point(self, n, z, from_particleList=False):
        if from_particleList:
            self.particleList = self.particleList[:n] + self.particleList[n + 1:]

        if abs(z - self.slice) <= self.radius:
            self.centimage.removeItem(self.circles_cent[n])

        self.circles_cent = self.circles_cent[:n] + self.circles_cent[n + 1:]

        self.leftimage.removeItem(self.circles_left[n])
        self.circles_left = self.circles_left[:n] + self.circles_left[n + 1:]

        self.bottomimage.removeItem(self.circles_bottom[n])
        self.circles_bottom = self.circles_bottom[:n] + self.circles_bottom[n + 1:]

    def remove_all_circles(self):
        for image in self.image_list:
            for child in image.allChildren():
                if type(child) == MyCircleOverlay:
                    image.removeItem(child)

        self.main_circles = []
        self.bottom_circles = []

    def load_images(self,folder='', bin_read=8, bin_alg=12):
        for name in ('findButton','detectButton','indexButton'):
            self.widgets[name].setEnabled(False)
        pg.QtGui.QApplication.processEvents()


        self.tomogram_name = self.settings.widgets['tomogram_name'].currentText()
        self.tilt_axis = float(self.settings.widgets['tilt_axis'].text())

        self.bin_read = int( self.settings.widgets['bin_read'].text() )
        self.bin_alg = int( self.settings.widgets['bin_alg'].text() )

        self.pixel_size = float(self.settings.widgets['pixel_size'].text())
        self.fiducial_size = float(self.settings.widgets['fiducial_size'].text() )
        self.settings.update_radius()

        self.algorithm = self.settings.widgets['algorithm'].currentText()


        folder = os.path.join(self.tomofolder,self.tomogram_name,'sorted/')

        self.xmin, self.ymin = 0, 0

        excl = os.path.join(folder,'excluded/')
        fnames = [[folder,line] for line in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted')] +\
                 [[excl,line]   for line in os.listdir( excl ) if line.endswith('.mrc') and line.startswith('sorted')]

        if len(fnames) == 0:
            print ('Directory is empty')
            return

        sort_str(fnames,1)

        for n, line in enumerate(fnames):
            fnames[n] = os.path.join(line[0],line[1])

        tiltfile = [os.path.join(folder,line) for line in os.listdir(folder) if line.endswith('.mdoc')][0]

        cmd = "cat {} | grep TiltAngle | sort -nk3 | awk '{{print $3}}'".format(tiltfile)

        self.tiltangles = numpy.array([float(line) for line in os.popen(cmd).readlines()], dtype=float)

        self.excluded = [0,]*len(fnames)
        for n, excl in enumerate(fnames):
            if 'excluded' == excl.split('/')[-2]:
                self.excluded[n] = 1
        self.read_data(fnames)

        self.loaded_data = True

        self.mark_frames = [[], ] * len(self.fnames)

        self.widgets['findButton'].setEnabled(True)

        self.replot2()

    def read_list(self, fnames, proc_id, nr_procs, frames, frames_full, frames_adj, transpose):
        print ("Start reading files process {}/{}".format(proc_id + 1, nr_procs))

        for nr, fname in enumerate(fnames):
            dataf = read_mrc('{}'.format(str(fname)), binning=[self.bin_read, self.bin_read, 1])
            if transpose: dataf = dataf.T
            w = wiener(dataf)
            hpf = dataf - gaussian_filter(dataf, 10) * 0
            data = downsample(w + hpf, int(round(self.bin_alg * 1. / self.bin_read)))
            data -= data.min()
            data /= data.max()

            dataf -= dataf.min()
            dataf /= dataf.max()
            frames[nr_procs * nr + proc_id] = (data) ** 0.75
            frames_full[nr_procs * nr + proc_id] = wiener(dataf)

    def read_data(self, fnames):

        s = time.time()
        self.fnames = fnames

        try:
            del self.frames
            del self.frames_full
        except:
            pass
        nr_procs = min(len(self.fnames), cpu_count() * 2)
        frames = ['', ] * len(fnames)
        frames_full = ['', ] * len(fnames)
        frames_full_adj = ['', ] * len(fnames)

        self.list_cx_cy_imnr = []
        self.coordinates = []
        self.user_coordinates = []
        self.mark_frames = -1 * numpy.ones((len(self.fnames), 1, 2))
        self.tiltimages = []

        manager = Manager()
        f = manager.list(frames)
        ff = manager.list(frames_full)
        ffa = manager.list(frames_full_adj)
        procs = []

        for proc_id in range(nr_procs):
            proc = Process(target=self.read_list,
                           args=(fnames[proc_id::nr_procs], proc_id, nr_procs, f, ff, ffa, True))
            procs.append(proc)
            proc.start()

        while len(procs):
            procs = [proc for proc in procs if proc.is_alive()]
            time.sleep(1)

        self.frames, self.frames_full, self.frames_adj = f, ff, ffa

        for n in range(len(self.fnames)):
            self.tiltimages.append(TiltImage(self, n, excluded=self.excluded[n]))

        self.imnr = numpy.abs(self.tiltangles).argmin()
        self.settings.widgets['ref_frame'].setValue( int(round(self.imnr)) )

        self.markermove = numpy.zeros_like(self.frames_full[self.imnr])

        self.dim = self.frames_full[0].shape[0]

        #self.update()

    def find_fid(self):
        #if not self.mf:
        #    self.find_fiducials()

        #else:
        #    if askokcancel(title="Overwrite existing fiducial list",
        #               message="Are you sure you want to overwrite the existing fiducial list by an automatd search? User changes will be lost."):

        if self.widgets['findButton'].isEnabled()==True:
            self.find_fiducials()
            self.replot2()

    def find_fiducials(self):
        if self.widgets['findButton'].isEnabled()==False: return
        print ('find potential fiducials time: ',)
        procs = []
        self.algorithm = self.settings.widgets['algorithm'].currentText()
        num_procs = min(len(self.fnames), cpu_count() * 2)

        self.list_cx_cy_imnr = []

        self.cent = numpy.zeros_like(self.frames_full[self.imnr])
        self.bin = 4
        self.mf = True

        s = time.time()

        fid_list = []
        manager = Manager()
        out = [0.] * num_procs
        for i in range(num_procs): out[i] = manager.list(fid_list)

        for proc_id in range(num_procs):
            if self.algorithm == 'sensitive':
                level = find_potential_fiducials_sensitive
            else:
                level = find_potential_fiducials

            proc = Process(target=level, args=( self.frames[proc_id::num_procs], self.frames_full[proc_id::num_procs],
                                                self.bin_alg, self.bin_read, 50, self.imnr, out[proc_id], proc_id,
                                                num_procs))
            procs.append(proc)
            proc.start()

        time.sleep(0.1)
        while len(procs):
            procs = [proc for proc in procs if proc.is_alive()]

        for i in range(num_procs):
            self.list_cx_cy_imnr += [el for el in out[i]]

        self.mark_frames = -1 * numpy.ones((len(self.fnames), 300, 2), dtype=float)

        sort(self.list_cx_cy_imnr, 1)

        cntr = numpy.zeros((200), dtype=int)


        # Clear existing fiducials in tiltimages
        for n in range(len(self.tiltimages)):
            self.tiltimages[n].clear()

        for cx, cy, imnr in self.list_cx_cy_imnr:
            CX,CY =cx*self.bin_alg*1./self.bin_read,cy*self.bin_alg*1./self.bin_read
            self.tiltimages[imnr].add_fiducial(CX-self.xmin,CY-self.ymin,CX,CY, check=False,draw=self.imnr==imnr)
            self.mark_frames[imnr][cntr[imnr]][:] = numpy.array((cy, cx))
            cntr[imnr] += 1

        self.mark_frames = self.mark_frames[:, :cntr.max(), :]

        self.coordinates = numpy.zeros_like(self.mark_frames)
        self.user_coordinates = numpy.zeros_like(self.mark_frames)

        self.bin = 0
        print (time.time()-s)
        self.widgets['detectButton'].setEnabled(True)

    def update_mark(self):

        self.mark_frames = -1 * numpy.ones((len(self.fnames), 300, 2), dtype=float)
        cntr = numpy.zeros((200), dtype=int)
        for tiltNr in range(len(self.fnames)):
            print(numpy.array(self.tiltimages[tiltNr].fiducials))
            temp_fid = numpy.array(self.tiltimages[tiltNr].fiducials)[:,:2]

            dx,dy = temp_fid.shape
            data = temp_fid.flatten()[::-1].reshape(dx,dy)[::-1]
            self.mark_frames[tiltNr][:len(data),:] = data
            cntr[tiltNr] += data.shape[0]

        self.mark_frames = self.mark_frames[:, :cntr.max(), :]/(1.*self.bin_alg/self.bin_read)

    def detect_frameshift(self):
        if self.widgets['detectButton'].isEnabled()==False: return

        self.update_mark()
        detect_shifts = detect_shifts_few
        if len(self.mark_frames[0]) > 5:
            detect_shifts = detect_shifts_many
        self.frame_shifts, self.numshifts, self.outline_detect_shifts, self.fs = detect_shifts(self.mark_frames,
                                                                                               diag=True,
                                                                                               image=self.frames[0])
        self.widgets['indexButton'].setEnabled(True)

    def index_fid(self):
        if self.widgets['indexButton'].isEnabled()==False: return

        self.deleteAllMarkers()
        self.update_mark()
        # Ensure all markers are frame shifted
        tx, ty, tz = self.mark_frames.shape
        temp = zeros((tx, ty, tz), dtype=float)
        cnt = zeros((tx, ty, 1))

        tiltnr, marknr, dummy = self.mark_frames.shape
        for itilt in range(tiltnr):
            for imark in range(marknr):
                if self.mark_frames[itilt, imark, 0] > 0.01 and self.mark_frames[itilt, imark, 1] > 0.01:
                    temp[itilt][imark] = self.fs[itilt]
                    cnt[itilt][imark] = 1
        print ('num_particles = ', cnt.sum())
        self.frame_shifts = temp
        ref_frame = int(self.settings.widgets['ref_frame'].value())
        tiltaxis = int(self.settings.widgets['tilt_axis'].value())
        self.coordinates, self.index_map, \
        self.frame_shifts_sorted, self.listdx = index_potential_fiducials(self.fnames, self.mark_frames,
                                                                          self.frame_shifts, tiltangles=self.tiltangles,
                                                                          plot=False, user_coords=self.user_coordinates,
                                                                          zero_angle=ref_frame, excluded=self.excluded,
                                                                          diag=True, add_marker=self.add_markers,cut=3.,
                                                                          tiltaxis=tiltaxis)

        for tiltNr in range(len(self.fnames)):
            self.tiltimages[tiltNr].update_indexing(self.coordinates[tiltNr]*1.*self.bin_alg/self.bin_read)


        # self.coordinates -= self.frame_shifts_sorted
        if len(self.user_coordinates) == 0:
            self.user_coordinates = numpy.zeros_like(self.coordinates)

        # try:
        #    self.recenter_fid()
        # except:
        #dd = (self.cmap + self.cmap + self.cmap + self.cmap + self.cmap + self.cmap + self.cmap)[:len(self.listdx)]
        #self.controller.cm.mclist1.clear()
        #self.controller.cm.mclist0.clear()

        #self.controller.cm.mclist0.reload(zip(self.controller.markers[1:], self.listdx, ["", ] * len(self.listdx)),
        #                                  color=dd[:len(self.listdx)])

        self.replot2()
        #self.controller.RecenterFid.config(state='active')

    def add_markers(self, data):
        self.selectMarkers.addMarker(self.selectMarkers.model, data)
        self.manual_adjust_marker.addMarker(self.manual_adjust_marker.MRMmodel, data)

    def deleteAllMarkers(self):
        self.selectMarkers.deleteAll()
        self.manual_adjust_marker.deleteAll()

    def recenter(self, markerFileName='markerfile_ref_TEMP.em', outFileName='markerfile_ref_TEMP.em',
                     tiltSeriesFormat='mrc', return_ref_coords=True, selected_markers=True, save=False):
        if not mf_write:
            print ("NO RECENTERING: loading pytom failed")
            return 0
        markerFileName = os.path.join(os.path.dirname(self.fnames[0]).replace('/excluded', ''), markerFileName)
        ret = self.save(markerFileName=markerFileName, selected_markers=selected_markers)
        if ret == 0: return ret


        # Exclude the excluded frames.
        projIndices = numpy.arange(len(self.frames))
        take = 1 - numpy.array(self.excluded, dtype=int)
        projIndices = list(projIndices[take.astype(bool)])
        prefix = '{}/{}/sorted/sorted_'.format(self.tomofolder, self.tomogram_name)

        ref_frame = int(self.settings.widgets['ref_frame'].text())
        # Calculate the refined fiducial positions.
        self.ref_coords, self.errors, tX, tY = refineMarkerPositions(prefix, markerFileName, 0,
                                                                     len(self.frames) - 1, outFileName, dimBox=64,
                                                                     projIndices=projIndices,
                                                                     tiltSeriesFormat=tiltSeriesFormat,
                                                                     ret=return_ref_coords, write=False,
                                                                     ireftilt=ref_frame)
        print(self.errors)
        self.ref_coords /= (self.bin_alg* 1.)

        self.old_coords = numpy.zeros_like(self.coordinates)



        # Update respective tables.
        num_markers = self.manual_adjust_marker.MRMmodel.rowCount()
        markIndices = []
        names = []
        num_fid = []
        if selected_markers:
            model = self.selectMarkers.model1
        else:
            model = self.selectMarkers.model

        selected_marker_IDs = range(model.rowCount())
        num_markers = len(selected_marker_IDs)
        markIndices = []
        for num, iid in enumerate(selected_marker_IDs):
            name = model.data(model.index(iid,0))
            names.append(name)
            num_fid.append(model.data(model.index(iid,1)) )
            imark = int(name.split('_')[-1])
            markIndices.append(imark)


        if len(self.errors):
            if selected_markers:
                model = self.selectMarkers.model1
            else:
                model = self.selectMarkers.model

            data = list(zip(names,num_fid,self.errors))
            model.removeRows(0,model.rowCount())

            

            self.selectMarkers.addMarker(model,data)


        self.update_mark()

        # Save the refined positions if save == True.

        if save:
            for index, imnr in enumerate(projIndices):
                for num, imark in enumerate(markIndices):

                    (xx, yy) = self.ref_coords[index][num]

                    for m, (fx, fy) in enumerate(self.mark_frames[imnr]):
                        if fx < 0.1 or fy < 0.1: continue
                        if abs(self.coordinates[imnr][imark][0] - fx) + abs(
                                self.coordinates[imnr][imark][1] - fy) < 0.5:
                            self.coordinates[imnr][imark] = [xx, yy]
                            self.mark_frames[imnr][m] = [xx, yy]
                            self.old_coords[imnr][imark] = [fx, fy]

        for name in (markerFileName, outFileName):
            if os.path.exists(name): os.system('rm {}'.format(name))


        # Update tiltimages to include refined fiducial locations.
        for tiltNr in range(len(self.fnames)):
            labels = copy.deepcopy(self.tiltimages[tiltNr].labels)
            index = copy.deepcopy(self.tiltimages[tiltNr].indexed_fiducials)
            
            self.tiltimages[tiltNr].clear()

            for n, (cx, cy) in enumerate( self.mark_frames[tiltNr] ):
                if cx < 0 and cy < 0: continue
                
                FY,FX = cx*self.bin_alg/self.bin_read, cy*self.bin_alg/self.bin_read
                self.tiltimages[tiltNr].add_fiducial(FX-self.xmin,FY-self.ymin,FX,FY,label=labels[n])
                self.tiltimages[tiltNr].indexed_fiducials[n] = index[n]

        self.replot2()

    def save(self, markerFileName='markerfile.em', ask_recent=1, selected_markers=True):

        output_type = markerFileName.split('.')[-1]
        markerFileName = os.path.join(os.path.dirname(self.fnames[0]).replace('/excluded', ''), markerFileName)
        
        tiltangles = self.tiltangles

        num_markers = self.manual_adjust_marker.MRMmodel.rowCount()
        markIndices = []

        if selected_markers:
            selected_marker_IDs = range(self.selectMarkers.model1.rowCount())
            num_markers = len(selected_marker_IDs)
            markIndices = []
            for num, iid in enumerate(selected_marker_IDs):
                name = self.selectMarkers.model1.data(self.selectMarkers.model1.index(iid,0))
                imark = int(name.split('_')[-1])
                markIndices.append(imark)
        else:
            for name in range(num_markers):
                name = self.manual_adjust_marker.MRMmodel.data(self.manual_adjust_marker.MRMmodel.index(iid, 0))
                imark = int(name.split('_')[-1])
                markIndices.append(imark)

        if num_markers == 0:
            print ('No markerfile selected', "No markers selected, no markerfile saved.")
            return num_markers

        projIndices = numpy.arange(len(self.frames))
        take = 1 - numpy.array(self.excluded, dtype=int)
        projIndices = list(projIndices[take.astype(bool)])
        locX, locY = 1, 2

        if output_type == 'mrc':
            markerFileVol = -1. * numpy.ones((num_markers, len(projIndices), 12), dtype='float64')
            for num, imark in enumerate(markIndices):
                for (itilt, iproj) in enumerate(projIndices):
                    markerFileVol[num][itilt][0] = self.tiltangles[iproj]
                    markerFileVol[num][itilt][locX] = int(self.coordinates[iproj][imark][1] * self.bin_alg)
                    markerFileVol[num][itilt][locY] = int(self.coordinates[iproj][imark][0] * self.bin_alg)
            convert_numpy_array3d_mrc(markerFileVol, markerFileName)
        elif output_type == 'em':
            
            markerFileVol = vol(12, len(projIndices), num_markers)
            markerFileVol.setAll(-1)
            for (imark, Marker) in enumerate(markIndices):
                for (itilt, TiltIndex) in enumerate(projIndices):
                    if self.coordinates[TiltIndex][Marker][1] < 1 and self.coordinates[TiltIndex][Marker][0] < 1:
                        continue

                    markerFileVol.setV(int(round(self.tiltangles[TiltIndex])), 0, itilt, imark)
                    markerFileVol.setV(int(round(self.coordinates[TiltIndex][Marker][1] * self.bin_alg)), locX,
                                       itilt, imark)
                    markerFileVol.setV(int(round(self.coordinates[TiltIndex][Marker][0] * self.bin_alg)), locY,
                                       itilt, imark)

            write_em(markerFileName, markerFileVol)

        print ('output_type: ', output_type, markerFileName)

    def load(self):
        pass


class SettingsFiducialAssignment(QMainWindow, CommonFunctions):
    def __init__(self,parent=None):
        super(SettingsFiducialAssignment, self).__init__(parent)
        self.setGeometry(900, 0, 300, 100)
        self.layout = self.grid = QGridLayout(self)
        self.setWindowTitle('Settings')
        self.settings = QWidget(self)
        self.settings.setLayout(self.layout)
        self.row, self.column = 0, 0
        self.logbook = {}
        self.widgets = {}
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label_combobox(self.grid, 'Tomogram Name', 'tomogram_name', self.parent().tomogram_names, rstep=1, cstep=-1,
                                   tooltip='Algorithm used for automatic fiducial detection.')

        self.insert_label_spinbox(self.grid, 'tilt_axis', text='Angle Tilt Axis (degrees)', rstep=1,
                               value=270, maximum=359, stepsize=5,
                               tooltip='Angle of the Tiltaxis, North direction = 0, West = 270.')

        self.insert_label_spinbox(self.grid, 'ref_frame', text='Reference Frame', rstep=1,
                                  value=19,minimum=1,maximum=91,stepsize=1,
                                  tooltip='Fiducial sets are calculated from the reference frame.')

        self.insert_label(self.grid,rstep=1)

        self.insert_label_spinbox(self.grid, 'pixel_size', text=r'Pixel Size (Å)', rstep=1,
                                  value=1.75,maximum=300, stepsize=1., wtype=QDoubleSpinBox,
                                  tooltip=r'Pixel Size in Ångstrom.')

        self.insert_label_spinbox(self.grid, 'fiducial_size', text='Fiducial Size (Å)', rstep=1,
                                  value=50, maximum=350, stepsize=10, tooltip='Size of Fiducuials in nanometer.')

        self.insert_label(self.grid,rstep=1)


        self.insert_label_spinbox(self.grid,'bin_read', text='Binning Factor Reading', rstep=1, minimum=1, maximum=16,
                                  stepsize=2,tooltip='Binning factor for reading.',value=2,wtype=QSpinBox,cstep=-1)

        self.insert_label_spinbox(self.grid,'bin_alg', text='Binning Factor Finding Fiducials',rstep=1,
                                  minimum=1,maximum=16,stepsize=2,value=8,wtype=QSpinBox,cstep=-1,
                                  tooltip='Binning factor for finding fiducials, used to improve contrast.')

        self.insert_label(self.grid,rstep=1)

        self.insert_label_combobox(self.grid, 'Accuracy level', 'algorithm', ['normal', 'sensitive'], rstep=1,cstep=-1,
                                   tooltip='Algorithm used for automatic fiducial detection.')

        self.insert_label(self.grid,rstep=1, cstep=1)

        self.insert_pushbutton(self.grid, text='Load Tilt Images', rstep=1, action=self.parent().load_images,params='',
                               tooltip='Load tilt images of tomogram set in settings.')

        self.setCentralWidget(self.settings)

        #CONNECT
        self.widgets['pixel_size'].valueChanged.connect(self.update_radius)
        self.widgets['fiducial_size'].valueChanged.connect(self.update_radius)
        self.widgets['ref_frame'].valueChanged.connect(self.update_ref_frame)
        self.update_radius()

    def update_ref_frame(self):
        self.parent().ref_frame = int(self.widgets['ref_frame'].text())

    def update_radius(self):
        w = self.widgets
        fiducial_size,pixel_size,bin_read = map(float,(w['fiducial_size'].text(),
                                                       w['pixel_size'].text(),
                                                       w['bin_read'].text()) )

        self.parent().radius = fiducial_size/(pixel_size*bin_read*2.)

        if self.parent().loaded_data: self.parent().replot2()

class SelectAndSaveMarkers(QMainWindow,CommonFunctions):
    def __init__(self,parent=None):
        super(SelectAndSaveMarkers,self).__init__(parent)
        self.setGeometry(670, 438, 600, 225)
        self.layout = self.grid = QGridLayout(self)
        self.general = QWidget(self)
        self.general.setLayout(self.layout)
        self.setWindowTitle('Select and Save Marker Sets')
        self.row, self.column = 4, 1
        self.logbook = {}
        self.widgets = {}
        rows, columns = 15,5
        self.header_names = ['Name','#Markers','Rscore']
        self.num_rows, self.num_columns = rows, len(self.header_names)
        self.items = [['', ] * columns, ] * rows

        self.dataGroupBox = QGroupBox("All Marker Sets")
        self.dataView = QTreeView()
        self.dataView.setRootIsDecorated(False)
        self.dataView.setAlternatingRowColors(True)
        self.dataView.setSortingEnabled(True)
        self.dataView.setSelectionMode(3)
        self.dataView.setEditTriggers(QAbstractItemView.NoEditTriggers)
        dataLayout = QHBoxLayout()
        dataLayout.addWidget(self.dataView)
        self.dataGroupBox.setLayout(dataLayout)
        self.model = self.createMailModel(self)
        self.dataView.setModel(self.model)

        self.selectGroupBox = QGroupBox("Selected Marker Sets")
        self.selectView = QTreeView()
        self.selectView.setRootIsDecorated(False)
        self.selectView.setAlternatingRowColors(True)
        self.selectView.setSortingEnabled(True)
        self.selectView.setSelectionMode(3)
        self.selectView.setEditTriggers(QAbstractItemView.NoEditTriggers)
        selectLayout = QHBoxLayout()
        selectLayout.addWidget(self.selectView)
        self.selectGroupBox.setLayout(selectLayout)
        self.model1 = self.createMailModel(self)
        self.selectView.setModel(self.model1)

        self.layout.addWidget(self.dataGroupBox,0,0,10,1)
        self.layout.addWidget(self.selectGroupBox, 0, 3, 10, 2)

        self.insert_pushbutton(self.grid,text='>>',rstep=1,action=self.move,params=[self.dataView,self.selectView])
        self.insert_pushbutton(self.grid,text='<<',rstep=8,cstep=2, action=self.move,
                               params=[self.selectView,self.dataView])
        self.insert_pushbutton(self.grid,text='Recenter Markers',cstep=1,action=self.parent().recenter,params=0)
        self.insert_pushbutton(self.grid,text='Save Markers',cstep=-4,action=self.parent().save,params=0)
        self.insert_pushbutton(self.grid, text='Load Markerfile', action=self.parent().load,width=150,params=0)

        self.setCentralWidget(self.general)

    def move(self,params):
        model = params[0]
        model1 = params[1]

        data = ['',]*self.num_columns
        tbr = []

        for n, index in enumerate(model.selectedIndexes()):
            data[index.column()] = index.data()
            if n%self.num_columns == self.num_columns-1:
                self.addMarker(model1.model(),data)
                tbr.append(index.row())

        for row in sorted(tbr,reverse=True):
            model.model().removeRow(row)

    def createMailModel(self,parent):
        model = QStandardItemModel(0, self.num_columns, parent)
        for i in range(len(self.header_names)):
            model.setHeaderData(i, Qt.Horizontal, self.header_names[i])
        return model

    def addMarker(self, model, data):
        model.insertRow(model.rowCount())
        for i in range(min(len(data),len(self.header_names) )):
            model.setData(model.index(model.rowCount()-1, i), data[i])
            for view in (self.dataView,self.selectView):
            #for i in range(self.num_columns):
                view.resizeColumnToContents(i)

    def deleteAll(self):
        for m in (self.model,self.model1):
            m.removeRows(0,m.rowCount())

class ManuallyAdjustMarkers(QMainWindow, CommonFunctions):
    def __init__(self, parent=None):
        super(ManuallyAdjustMarkers, self).__init__(parent)
        self.setGeometry(670, 0, 600, 200)
        self.layout = self.grid = QGridLayout(self)
        self.general = QWidget(self)
        self.general.setLayout(self.layout)
        self.setWindowTitle('Select and Save Marker Sets')
        self.header_names = ['Name','#Markers']
        self.row, self.column = 0, 1
        self.logbook = {}
        self.widgets = {}
        rows, columns = 20, 4
        self.num_rows, self.num_columns = rows, columns
        self.items = [['', ] * columns, ] * rows

        self.dataGroupBox = QGroupBox("All Marker Sets")
        self.dataView = QTreeView()
        self.dataView.setRootIsDecorated(False)
        self.dataView.setAlternatingRowColors(True)
        self.dataView.setSortingEnabled(True)
        self.dataView.setSelectionMode(3)
        self.dataView.setEditTriggers(QAbstractItemView.NoEditTriggers)
        dataLayout = QHBoxLayout()
        dataLayout.addWidget(self.dataView)
        self.dataGroupBox.setLayout(dataLayout)
        self.MRMmodel = self.createMarkResultModel(self)
        self.dataView.setModel(self.MRMmodel)

        self.layout.addWidget(self.dataGroupBox, 0, 0, 20, 1)

        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_pushbutton(self.grid, text='Add Marker', rstep=1, action=self.add_marker, cstep=0)
        self.insert_pushbutton(self.grid, text='Delete Marker', rstep=1, action=self.delete_marker, cstep=0)
        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_pushbutton(self.grid, text='Select Marker', rstep=1, action=self.select_marker, cstep=0)
        self.insert_pushbutton(self.grid, text='Deselect Marker', rstep=1, action=self.deselect_marker, cstep=0)
        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_pushbutton(self.grid, text='Next Missing Fiducial', rstep=1, cstep=0, action=self.next_missing)
        self.insert_pushbutton(self.grid, text='Prev Missing Fiduical', rstep=1, cstep=0, action=self.prev_missing)
        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_label(self.grid, rstep=1, cstep=0)

        self.setCentralWidget(self.general)

    def createMarkResultModel(self, parent):
        model = QStandardItemModel(0, len(self.header_names), parent)
        for i in range(len(self.header_names)):
            model.setHeaderData(i, Qt.Horizontal, self.header_names[i])
        return model

    def addMarker(self, model, data):
        self.MRMmodel.insertRow(self.MRMmodel.rowCount())
        for i in range(min(len(data),self.num_columns)):
            self.MRMmodel.setData(model.index(self.MRMmodel.rowCount() - 1, i), data[i])
            self.dataView.resizeColumnToContents(i)

    def deleteAll(self):
        self.MRMmodel.removeRows(0,self.MRMmodel.rowCount())

    def add_marker(self, params=None):
        self.addMarker(self.model,['Marker_{:03d}'.format(self.model.rowCount()), ''])
        pass

    def delete_marker(self, params=None):
        pass

    def select_marker(self, params=None):
        try:
            print (self.dataView.selectedIndexes()[0].data()    , ' selected')
            self.parent().selected_marker = int( self.dataView.selectedIndexes()[0].data().split("_")[-1] )
        except: pass

    def deselect_marker(self, params=None):
        if self.parent().selected_marker > -1:
            print ('Marker_{:03d} deselected.'.format(self.parent().selected_marker))
            self.parent().selected_marker = -1
        else:
            pass

    def next_missing(self, params=None):
        pass

    def prev_missing(self,params=None):
        pass

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setWindowIcon(QIcon('/Users/gijs/Documents/PostDocUtrecht/GUI/pp.jpg'))
    ex = FiducialAssignment()
    ex.show()
    sys.exit(app.exec_())

if __name__=='__main__':
    main()
