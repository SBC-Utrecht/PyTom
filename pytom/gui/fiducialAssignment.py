import copy
import atexit
import sys
import os
import time
import mrcfile

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

# pytom gui imports
from pytom.gui.fiducialPicking import PickingFunctions, projectMarkerToFrame
from pytom.gui.guiStructures import circle, CommonFunctions, MyCircleOverlay, Worker
from pytom.gui import guiFunctions
from pytom.gui.mrcOperations import downsample, convert_numpy_array3d_mrc

import numpy as np
from scipy.signal import wiener
from multiprocessing import Process, Manager, cpu_count
from pytom.agnostic.io import read_size

# dont know why this is all needed...
global mf_write
try:
    from pytom.reconstruction.markerPositionRefinement import refineMarkerPositions
    from pytom.reconstruction.tiltAlignmentFunctions import alignmentFixMagRot
    from pytom.reconstruction.TiltAlignmentStructures import Marker
    from pytom.basic.files import write_em
    from pytom_volume import vol, read, transform
    from pytom_numpy import vol2npy
    mf_write=1
except Exception as e:
    print(e)
    print('marker file refinement is not possible')
    mf_write = 0


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

    def add_fiducial(self, x, y, FX, FY, cutoff=5, check=True, draw=True, label='', color=Qt.red):
        added = False
        new = True
        radius = self.parent.radius

        if check:
            for n, (fx, fy, fi) in enumerate(self.fiducials):
                if fx > 0 and fy > 0 and abs(x-fx) < cutoff and abs(y-fy) < cutoff:
                    new = False
                    return new


        self.fiducials.append([FX, FY, self.imnr])
        self.indexed_fiducials.append( color )
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


class FiducialAssignment(QMainWindow, CommonFunctions, PickingFunctions ):
    def __init__(self, parent=None, fname=''):
        super(FiducialAssignment, self).__init__(parent)
        self.size_policies()

        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname

        self.threadPool = self.parent().threadPool
        self.threadPool.setMaxThreadCount(1)
        self.activeProcesses = {}
        self.ID = 0

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
        self.metafile = ''
        self.settings = SettingsFiducialAssignment(self)
        self.selectMarkers = SelectAndSaveMarkers(self)
        self.manual_adjust_marker = ManuallyAdjustMarkers(self)
        self.error_window = ErrorWindow(self,'alignmentErrors.txt')
        self.markerAdjust = False
        self.markerHotKeys = {Qt.Key_4: 4, Qt.Key_5: 5, Qt.Key_6: 6, Qt.Key_7: 7, Qt.Key_8: 8, Qt.Key_9: 9, Qt.Key_0: 0}

        self.idReferenceImage = -1
        self.dim = 0
        self.radius = 8
        self.sizeCut = 200
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.xmin, self.ymin = 0,0
        self.selected_marker = -1
        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.main_canvas = w0 = pg.GraphicsWindow(size=(600, 600), border=True)
        self.main_image = w0.addViewBox(row=0, col=0, lockAspect=True)
        self.main_image.setMouseEnabled(False)
        self.main_image.invertY(True)

        self.actionwidget = QWidget(self,width=300,height=300)
        self.al = parent = QGridLayout(self)
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.widgets = {}
        self.insert_label(parent, rstep=1)
        self.insert_pushbutton(parent,text='Find Fiducials',rstep=1,tooltip='Automatically detect fiducials.',
                               action=self.find_fid,params=0, wname='findButton', state=False)

        self.insert_pushbutton(parent,text='Detect Frame Shifts',rstep=1, wname='detectButton',
                               tooltip='Detect global x-y shifts between tilt images.',
                               action=self.detect_frameshift,params=0, state=False)
        self.insert_pushbutton(parent,text='Index Fiducials',rstep=1,tooltip='Group Fiducials into marker sets.',
                               action=self.index_fid,params=0, wname='indexButton',state=False)
        self.insert_pushbutton(parent, text='Check Align Errors', rstep=1, tooltip='Check Alignment Errors.',
                               action=self.raise_window, params=self.error_window, wname='errorButton', state=False)
        self.insert_label(parent, rstep=1)
        self.insert_pushbutton(parent,text='Manually Adjust Markers',rstep=1,tooltip='Manually adjust marker sets.',
                               action=self.raise_window, params=self.manual_adjust_marker, wname='adjustManually',
                               state=False)
        self.insert_pushbutton(parent, text='Create Markerfile', rstep=1, tooltip='Select and save marker sets.',
                               action=self.raise_window, params=self.selectMarkers,wname='createMarkerfile',state=False)
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

        self.img1a = pg.ImageItem(np.zeros((10,10)))
        #self.img1b = pg.ImageItem(np.zeros((10,10)))
        self.img1c = pg.ImageItem(np.zeros((10,10)))

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
        self.settings.show()
        pg.QtGui.QApplication.processEvents()
        self.loaded_data = False


    def raise_window(self,window):
        try:
            window.close()
            window.show()
        except Exception as e:
            print(e)
        pass

    def selectPartTiltImage(self,event):
        if self.idReferenceImage == -1: return

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

    def updatePartTiltImageGuess(self, x, y):

        self.bottom_circles.append(
            circle(QPoint(x - self.xmin - self.radius, y - self.ymin - self.radius), size=self.radius * 2,
                   color=Qt.yellow))
        self.bottom_image.addItem(self.bottom_circles[-1])

    def manuallyAdjustFiducials(self,event):
        if self.idReferenceImage == -1:
            return

        if not self.dim: return
        pos = self.bottom_image.mapSceneToView(event.scenePos())
        x, y = int(round(pos.x())), int(round(pos.y()))

        if x>0 and y>0 and y<self.sizeCut and x < self.sizeCut:
            fx, fy = x + self.xmin, y + self.ymin

            if event.button() == 1:
                if self.selected_marker > -1:
                    px,py = self.coordinates[self.imnr,self.selected_marker]*self.bin_alg / self.bin_read
                    if px > 0.1 or py > 0.1:
                        self.tiltimages[self.imnr].remove_fiducial(py, px, cutoff=self.radius)
                    self.coordinates[self.imnr,self.selected_marker] = [fy*self.bin_read/self.bin_alg,fx*self.bin_read/self.bin_alg]
                    self.tiltimages[self.imnr].add_fiducial(x, y, fx, fy, cutoff=self.radius,color=Qt.green,label="{:03d}".format(self.selected_marker))
                    self.replot2()
                else:
                    self.tiltimages[self.imnr].add_fiducial(x,y,fx,fy,cutoff=self.radius)


            elif event.button() == 2:
                self.tiltimages[self.imnr].remove_fiducial(fx, fy, cutoff=self.radius)

                #THE CHILDREN DO NOT HAVE REPRODUCIBLE COORDINATES
                self.remove_all_circles()
                self.add_circles_new_frame()

    def wheelEvent(self, event):
        if self.idReferenceImage == -1: return
        step = event.angleDelta().y() // 120

        if self.imnr + step < len(self.fnames) and self.imnr + step > -1:
            self.imnr += step
            self.replot()

    def markerHotKeyPressed(self, id):
        if self.selected_marker != id:
            self.selected_marker = id
            print(f'Marker_{id:03d} selected')
        else:
            print(f'Marker_{id:03d} deselected')
            self.selected_marker = -1

    def keyPressEvent(self, evt):
        if Qt.Key_F == evt.key():
            if self.idReferenceImage == -1: return

            print ('Find Fiducials')
            self.find_fid()

        elif Qt.Key_D == evt.key():
            if self.idReferenceImage == -1: return

            print ('Detect Frame Offsets')
            self.detect_frameshift()

        elif Qt.Key_I == evt.key():
            if self.idReferenceImage == -1: return

            print ('Index Fiducials')
            self.index_fid()

        elif Qt.Key_1 == evt.key():
            if not self.markerAdjust:
                self.imnr = 0
                self.replot2()
            else:
                self.markerHotKeyPressed(1)

        elif Qt.Key_2 == evt.key():
            if not self.markerAdjust:
                self.imnr = self.idReferenceImage
                self.replot2()
            else:
                self.markerHotKeyPressed(2)

        elif Qt.Key_3 == evt.key():
            if self.idReferenceImage == -1: return
            if not self.markerAdjust:
                self.imnr = len(self.fnames) -1
                self.replot2()
            else:
                self.markerHotKeyPressed(3)

        elif evt.key() in self.markerHotKeys.keys():
            if self.idReferenceImage == -1: return

            if self.markerAdjust:
                self.markerHotKeyPressed(self.markerHotKeys[evt.key()])

        elif Qt.Key_L == evt.key():
            if self.idReferenceImage == -1: return

            for i in range(len(self.fnames)):
                self.imnr = i
                self.replot2()
                pg.QtGui.QApplication.processEvents()
                time.sleep(0.04)

            for i in range(len(self.fnames)):
                self.imnr -= 1
                self.replot2()
                pg.QtGui.QApplication.processEvents()
                time.sleep(0.04)

        elif Qt.Key_P == evt.key():
            if self.idReferenceImage == -1: return

            save_index = self.imnr
            for i in range(len(self.fnames)):
                self.imnr = i
                self.replot2()
                from pyqtgraph.exporters import ImageExporter
                exporter = ImageExporter( self.main_image )
                exporter.export('assignedMarkers_{:02d}.png'.format(i))

            self.imnr = save_index
            self.replot2()

        elif Qt.Key_E == evt.key():
            if self.idReferenceImage == -1: return

            self.exclude_status_change()

        elif Qt.Key_R == evt.key():
            print ('Refine Marker Positions')

        elif Qt.Key_Right == evt.key():
            if self.idReferenceImage == -1: return

            if self.imnr + 1 < len(self.fnames):
                self.imnr += 1
                self.replot2()

        elif Qt.Key_Left == evt.key():
            if self.idReferenceImage == -1: return

            if self.imnr - 1 >= 0:

                self.imnr -= 1
                self.replot2()

        elif evt.key() == Qt.Key_Escape:
            self.close()

        elif Qt.Key_Comma == evt.key():
            if self.idReferenceImage == -1: return

            self.manual_adjust_marker.prev_missing()

        elif Qt.Key_Period == evt.key():
            if self.idReferenceImage == -1: return
            try:
                self.manual_adjust_marker.next_missing()
            except:
                pass

        elif Qt.Key_M == evt.key():
            if self.idReferenceImage == -1: return

            if self.widgets['adjustManually'].isEnabled():
                self.markerAdjust = (self.markerAdjust != True)
                mah = 'activated' if self.markerAdjust else 'deactivated'
                print(f'Marker Adjustment Hotkeys are {mah}.')
                if mah == 'deactivated': self.disableMarkerAdjustment()

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
        if self.idReferenceImage == -1: return
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
        self.disableMarkerAdjustment()
        self.idReferenceImage = -1

        for name in ('findButton','indexButton', 'detectButton', 'adjustManually', 'errorButton'):
            self.widgets[name].setEnabled(False)
        pg.QtGui.QApplication.processEvents()

        self.tomogram_name = self.settings.widgets['tomogram_name'].currentText()
        folder = os.path.join(self.tomofolder, self.tomogram_name, 'sorted/')
        excl = os.path.join(folder,'excluded/')
        fnames = [[folder,line] for line in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted')] +\
                 [[excl,line]   for line in os.listdir( excl ) if line.endswith('.mrc') and line.startswith('sorted')]

        self.metafile = [os.path.join(folder, line) for line in os.listdir(folder) if line.endswith('.meta')][0]
        try:
            self.metadata = guiFunctions.loadstar(self.metafile,dtype=guiFunctions.datatype)
        except:
            metadata_old = guiFunctions.loadstar(self.metafile,dtype=guiFunctions.datatype0)
            self.metadata = np.rec.array([(0.,)*len(guiFunctions.datatype),]*len(fnames), dtype=guiFunctions.datatype)

            for key, value in guiFunctions.datatype0:
                print(key, value)
                self.metadata[key] = metadata_old[key]

        self.tilt_angles = self.metadata['TiltAngle']
        ps,fs = float(self.metadata['PixelSpacing'][0]), int(self.metadata['MarkerDiameter'][0])
        self.settings.widgets['tilt_axis'].setValue(int(self.metadata['InPlaneRotation'][0]))
        self.settings.widgets['pixel_size'].setValue(ps)
        self.settings.widgets['fiducial_size'].setValue(fs)


        self.tilt_axis = float(self.settings.widgets['tilt_axis'].text())
        self.bin_read = int( self.settings.widgets['bin_read'].text() )
        self.bin_alg = int( self.settings.widgets['bin_alg'].text() )
        self.pixel_size = float(self.settings.widgets['pixel_size'].text())
        self.fiducial_size = float(self.settings.widgets['fiducial_size'].text() )

        self.settings.update_radius()

        self.algorithm = self.settings.widgets['algorithm'].currentText()

        self.xmin, self.ymin = 0, 0



        if len(fnames) == 0:
            print ('Directory is empty')
            return

        sort_str(fnames,1)

        for n, line in enumerate(fnames):
            fnames[n] = os.path.join(line[0],line[1])


        self.tiltangles = self.metadata['TiltAngle']
        #cmd = "cat {} | grep TiltAngle | sort -nk3 | awk '{{print $3}}'".format(tiltfile)
        #self.tiltangles = np.array([float(line) for line in os.popen(cmd).readlines()], dtype=float)

        self.excluded = [0,]*len(fnames)
        for n, excl in enumerate(fnames):
            if 'excluded' == excl.split('/')[-2]:
                self.excluded[n] = 1

        self.read_data(fnames)

    def read_list(self, fnames, proc_id, nr_procs, frames, frames_full, dataRaw, dummy):

        print ("Start reading files process {}/{}".format(proc_id + 1, nr_procs))
        from copy import deepcopy
        from pytom.agnostic.io import read as readNPY
        for nr, fname in enumerate(fnames):

            #m = mrcfile.open(fname, permissive=True)
            #dataf = copy.deepcopy(m.data)
            #m.close()
            #dataf = downsample(dataf,self.bin_read)
            #dataf = read_mrc('{}'.format(str(fname)), binning=[1, 1, 1])
            #frames_adj[nr_procs * nr + proc_id] = mrcfile.open(fname,permissive=True)

            i = nr*nr_procs + proc_id
            datafile = readNPY(self.fnames[i]).squeeze()

            dataRaw[i] = deepcopy(datafile)
            #datafile = read_mrc('{}'.format(self.fnames[i]), binning=[1, 1, 1])
            fa = deepcopy(datafile)
            fa[fa > fa.mean() + 5 * fa.std()] = np.median(fa)
            self.frames_adj[i, :, :] = fa
            dataf = self.frames_adj[nr*nr_procs + proc_id, :,:].copy()

            dataf = downsample(dataf, self.bin_read)

            w = wiener(dataf)
            hpf = dataf #- gaussian_filter(dataf, 10) * 0
            data = downsample(w + hpf, int(round(self.bin_alg * 1. / self.bin_read)))
            data -= data.min()
            data /= data.max()

            dataf -= dataf.min()
            dataf /= dataf.max()
            frames[nr_procs * nr + proc_id] = (data) ** 0.75
            frames_full[nr_procs * nr + proc_id] = wiener(dataf)

    def read_data(self, fnames):
        self.selected_marker = -1
        from copy import deepcopy

        self.fnames = fnames

        try:
            del self.frames
            del self.frames_full
            del self.frames_adj
        except:
            pass
        nr_procs = min(len(self.fnames), cpu_count() )

        frames = ['', ] * len(fnames)
        frames_full = ['', ] * len(fnames)
        frames_full_adj = ['', ] * len(fnames)

        self.list_cx_cy_imnr = []
        self.coordinates = []
        self.user_coordinates = []
        self.mark_frames = -1 * np.ones((len(self.fnames), 1, 2))
        self.tiltimages = []

        manager = Manager()
        f = manager.list(frames)
        ff = manager.list(frames_full)
        raw = manager.list(frames_full_adj)


        temp = mrcfile.open(fnames[0], permissive=True)
        dimx, dimy = temp.data.shape
        temp.close()
        self.frames_adj = np.zeros((len(self.fnames), dimy, dimx))
        self.dataRaw = np.zeros_like(self.frames_adj)

        procs = []

        self.error_window.logfile = os.path.join(os.path.dirname(fnames[0]), 'alignmentErrors.txt')

        for proc_id in range(nr_procs):
            proc = Process(target=self.read_list,
                           args=(fnames[proc_id::nr_procs], proc_id, nr_procs, f, ff, raw, True))
            procs.append(proc)
            proc.start()
            atexit.register(guiFunctions.kill_proc, proc)

        while len(procs):
            procs = [proc for proc in procs if proc.is_alive()]
            time.sleep(1)



        self.frames = np.zeros((len(fnames), f[0].shape[0], f[0].shape[1]) )
        self.frames_full = np.zeros((len(fnames),ff[0].shape[0],ff[0].shape[1]) )
        for i in range(len(fnames)):
            self.frames[i, :, :] = f[i]
            self.frames_full[i, :, :] = ff[i]
            self.dataRaw[i,:,:] = raw[i]

        for n in range(len(self.fnames)):
            self.tiltimages.append(TiltImage(self, n, excluded=self.excluded[n]))

        self.imnr = np.abs(self.tiltangles).argmin()
        self.settings.widgets['ref_frame'].setValue( int(round(self.imnr)) )

        self.markermove = np.zeros_like(self.frames_full[self.imnr])

        self.dim = self.frames_full[0].shape[0]
        #self.update()


        # finsihg
        self.loaded_data = True

        self.mark_frames = [[], ] * len(self.fnames)

        self.idReferenceImage = np.abs(self.metadata['TiltAngle']).argmin()

        self.widgets['findButton'].setEnabled(True)
        self.widgets['createMarkerfile'].setEnabled(True)
        self.replot2()

    def submitJob(self, job, args=[]):
        self.activeProcesses[self.ID] = Worker(fn=job, args=args, sig=False)
        self.threadPool.start(self.activeProcesses[self.ID])
        self.ID += 1

    def disableMarkerAdjustment(self):
        if self.markerAdjust:
            print(f'Marker Adjustment and Hotkeys are deactivated.')
        self.selected_marker = -1
        self.markerAdjust = False

    def find_fid(self):
        #if not self.mf:
        #    self.find_fiducials()

        #else:
        #    if askokcancel(title="Overwrite existing fiducial list",
        #               message="Are you sure you want to overwrite the existing fiducial list by an automatd search? User changes will be lost."):
        self.disableMarkerAdjustment()
        if self.widgets['findButton'].isEnabled()==True:
            # self.activeProcesses[self.ID] = Worker(fn=self.find_fiducials, args=[], sig=False)
            # self.threadPool.start(self.activeProcesses[self.ID])
            # self.ID += 1
            self.find_fiducials()

    def create_average_markers(self, start_frame, end_frame, display=False, markersize=128):
        r = markersize // 2
        dimx, dimy = self.frames_adj[0,:,:].shape
        image = np.zeros((r, r))
        total = 0

        for nn in range(start_frame, end_frame):

            for y, x in self.mark_frames[nn, :, :]*self.bin_alg:

                if x < 0 or y<0: continue
                y = int(round(y))
                x = int(round(x))

                if x > r // 2 and y > r // 2 and x < dimx - r // 2 and y < dimy - r // 2:
                    total += 1
                    print(x,y)
                    imp = self.frames_adj[nn, x - r // 2:x + r // 2, y - r // 2:y + r // 2]
                    print(imp.sum())
                    if image.sum():
                        ar = np.fft.fftshift(np.fft.fftn(imp))
                        br = np.conj(np.fft.fftshift(np.fft.fftn(image)))
                        dd = np.fft.fftshift(abs(np.fft.ifftn(np.fft.fftshift(ar * br))))

                        ddd = dd.flatten()
                        cx, cy = ddd.argmax() // r - r // 2, ddd.argmax() % r - r // 2

                        imp = self.frames_adj[nn, max(0, x + cx - r // 2):x + cx + r // 2, max(0, y + cy - r // 2):y + cy + r // 2]



                    if imp.shape[0] == image.shape[0] and imp.shape[1] == image.shape[1]:
                        image += imp

        image /= total
        print(total, image.sum())
        imgx, imgy = np.array(image.shape) // 2
        marker = np.ones_like(self.frames_adj[0,:,:]) * np.median(image)
        mx, my = np.array(marker.shape) // 2
        marker[mx - imgx:mx + imgx, my - imgy:my + imgy] = image
        from pylab import imshow, show

        imshow(marker)
        show()
        return marker

    def find_fiducials(self):
        print('Active Threads: ', self.threadPool.activeThreadCount())

        self.idReferenceImage = -1
        self.selected_marker = -1
        if self.widgets['findButton'].isEnabled()==False: return

        self.disableMarkerAdjustment()

        print ('find potential fiducials ',)
        procs = []
        self.algorithm = self.settings.widgets['algorithm'].currentText()
        num_procs = min(len(self.fnames), cpu_count() * 2)


        if self.settings.widgets['applyToAll'].isChecked():
            self.list_cx_cy_imnr = []

        self.cent = np.zeros_like(self.frames_full[self.imnr])
        self.bin = 4
        self.mf = True


        fid_list = []
        manager = Manager()
        out = [0.] * num_procs
        for i in range(num_procs): out[i] = manager.list(fid_list)
        ref_frame = int(self.settings.widgets['ref_frame'].value())
        threshold = float(self.settings.widgets['threshold'].value())
        radius = float(self.settings.widgets['radius_watershed'].value())
        average_marker = None
        self.mark_frames = -1 * np.ones((len(self.fnames), 3000, 2), dtype=float)
        # self.assignedFiducials = -1 * np.ones((len(self.fnames), 3000), dtype=int)

        if len(self.mark_frames[ref_frame]) > 2 and self.algorithm == 'cross_correlation':
            average_marker = self.create_average_markers(ref_frame - 5, ref_frame + 6)
        for proc_id in range(num_procs):

            if self.algorithm == 'cross_correlation':
                if len(self.mark_frames[ref_frame]) > 2:
                    level = self.find_potential_fiducials_crosscorr
                else:
                    level = self.find_potential_fiducials

            elif self.algorithm == 'sensitive':
                level = self.find_potential_fiducials_sensitive

            elif self.algorithm == 'normal':
                level = self.find_potential_fiducials

            elif self.algorithm == 'watershed':
                level = self.find_potential_fiducial_stdv

            elif self.algorithm == 'LoG':
                level = self.find_potential_fiducial_log

            else: return


            if not self.settings.widgets['applyToAll'].isChecked():
                frame = self.frames[self.imnr:self.imnr+1, :,:]
                frame_full = self.frames_full[self.imnr:self.imnr+1,:,:]
                raw = self.dataRaw[self.imnr:self.imnr+1,:,:]
                level( frame, frame_full, raw, self.bin_alg, self.bin_read, 50, self.imnr, out[0], self.imnr, 1, average_marker, threshold, radius )
                additional = [[el[0], el[1],self.imnr] for el in out[0]]
                print((additional))
                self.tiltimages[self.imnr].clear()
                for cx, cy, imnr in self.list_cx_cy_imnr:
                    if imnr != self.imnr:
                        additional.append([cx,cy,imnr])
                self.list_cx_cy_imnr = additional
                break

            else:
                proc = Process(target=level, args=( self.frames[proc_id::num_procs], self.frames_full[proc_id::num_procs],
                                                    self.dataRaw[proc_id::num_procs],
                                                    self.bin_alg, self.bin_read, 50, self.imnr, out[proc_id], proc_id,
                                                    num_procs, average_marker, threshold, radius ))
                procs.append(proc)
                proc.start()
                atexit.register(guiFunctions.kill_proc, proc)

        #time.sleep(0.1)

        if self.settings.widgets['applyToAll'].isChecked():
            while len(procs):
                procs = [proc for proc in procs if proc.is_alive()]

            for i in range(num_procs):
                self.list_cx_cy_imnr += [el for el in out[i]]


        guiFunctions.sort(self.list_cx_cy_imnr, 1)

        cntr = np.zeros((200), dtype=int)

        # Clear existing fiducials in tiltimages
        for n in range(len(self.tiltimages)):
            self.tiltimages[n].clear()

        for cx, cy, imnr in self.list_cx_cy_imnr:
            CX,CY =cx*self.bin_alg*1./self.bin_read,cy*self.bin_alg*1./self.bin_read
            self.tiltimages[imnr].add_fiducial(CX-self.xmin,CY-self.ymin,CX,CY, check=False,draw=self.imnr==imnr)
            self.mark_frames[imnr][cntr[imnr]][:] = np.array((cy, cx))
            # self.assignedFiducials[imnr][cntr[imnr]] = 0
            cntr[imnr] += 1

        self.mark_frames = self.mark_frames[:, :cntr.max(), :]
        # self.assignedFiducials = self.assignedFiducials[:, :cntr.max()]

        self.coordinates = np.zeros_like(self.mark_frames)
        self.user_coordinates = np.zeros_like(self.mark_frames)

        self.bin = 0
        self.widgets['detectButton'].setEnabled(True)
        self.widgets['createMarkerfile'].setEnabled(True)
        print('fid finding jobs finished')
        self.idReferenceImage = np.abs(self.metadata['TiltAngle']).argmin()

        self.replot2()

    def update_mark(self):

        self.mark_frames = -1 * np.ones((len(self.fnames), 3000, 2), dtype=float)
        cntr = np.zeros((200), dtype=int)
        for tiltNr in range(len(self.fnames)):
            #print(np.array(self.tiltimages[tiltNr].fiducials))
            try: temp_fid = np.array(self.tiltimages[tiltNr].fiducials)[:,:2]
            except: continue
            dx,dy = temp_fid.shape
            data = temp_fid.flatten()[::-1].reshape(dx,dy)[::-1]
            self.mark_frames[tiltNr][:len(data),:] = data
            cntr[tiltNr] += data.shape[0]

        self.mark_frames = self.mark_frames[:, :cntr.max(), :]/(1.*self.bin_alg/self.bin_read)

    def detect_frameshift(self):
        print('Activce Threads: ', self.threadPool.activeThreadCount())
        if self.widgets['detectButton'].isEnabled()==False: return
        self.disableMarkerAdjustment()
        self.update_mark()
        detect_shifts = self.detect_shifts_few
        if len(self.mark_frames[0]) > 5:
            detect_shifts = self.detect_shifts_many
        self.frame_shifts, self.numshifts, self.outline_detect_shifts, self.fs = detect_shifts(self.mark_frames,
                                                                                               diag=True,
                                                                                               image=self.frames[0])
        self.widgets['indexButton'].setEnabled(True)

    def index_fid(self):
        self.selected_marker = -1
        if self.widgets['indexButton'].isEnabled()==False: return
        self.disableMarkerAdjustment()

        self.deleteAllMarkers()
        self.update_mark()
        # Ensure all markers are frame shifted
        tx, ty, tz = self.mark_frames.shape
        temp = np.zeros((tx, ty, tz), dtype=float)
        cnt = np.zeros((tx, ty, 1))

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
        max_shift = 40#self.radius*2.*self.bin_read/self.bin_alg
        print(max_shift)

        reference_marker = int(self.settings.widgets['ref_marker'].text())

        self.coordinates, self.index_map, \
        self.frame_shifts_sorted, self.listdx = self.index_potential_fiducials(self.fnames, self.mark_frames, self.frame_shifts, tiltangles=self.tiltangles,
                                                                          plot=False, user_coords=self.user_coordinates,
                                                                          zero_angle=ref_frame, excluded=self.excluded,
                                                                          diag=True, add_marker=self.add_markers,
                                                                          cut=max_shift, tiltaxis=tiltaxis)
        #np.save(os.path.join(self.projectname,'mark_frames.npy'), self.mark_frames)

        projIndices = np.array(list(range(self.coordinates.shape[0])))[(1-np.array(self.excluded)) > 0.5]

        imdimX, imdimY, z = read_size(self.fnames[0])

        cc = self.coordinates * self.bin_alg
        cc[cc < 0] = -1

        a = self.calculate_error(self.coordinates*self.bin_alg, self.tiltangles, projIndices, imdimX, imdimY,
                             fname=os.path.join(os.path.dirname(self.metafile), 'alignmentErrors.txt'), reference_marker=reference_marker)

        self.errorsModel, self.shiftXModel, self.shiftYModel, self.diffXModel, self.diffYModel, x,y,z, psi = a


        self.frame_shifts = np.array(list(zip(self.shiftXModel, self.shiftYModel)))

        self.centersModel = list(zip(x,y,z))
        self.psiindeg = psi

        try:
            self.deleteAllMarkers()
            self.update_mark()
            self.determine_markerdata(self.coordinates, self.errorsModel, self.excluded, self.add_markers)
        except Exception as e:
            print(e)
            pass

        for tiltNr in range(len(self.fnames)):
            self.tiltimages[tiltNr].update_indexing(self.coordinates[tiltNr]*1.*self.bin_alg/self.bin_read)


        # self.coordinates -= self.frame_shifts_sorted
        if len(self.user_coordinates) == 0:
            self.user_coordinates = np.zeros_like(self.coordinates)

        # try:
        #    self.recenter_fid()
        # except:
        #dd = (self.cmap + self.cmap + self.cmap + self.cmap + self.cmap + self.cmap + self.cmap)[:len(self.listdx)]
        #self.controller.cm.mclist1.clear()
        #self.controller.cm.mclist0.clear()

        #self.controller.cm.mclist0.reload(zip(self.controller.markers[1:], self.listdx, ["", ] * len(self.listdx)),
        #                                  color=dd[:len(self.listdx)])

        self.replot2()

        #self.error_window.dataView.selectionModel().setCurrentIndex(self.error_window.EWmodel.index(self.error_window.EWmodel.rowCount() - 1, 0))


        #self.controller.RecenterFid.config(state='active')
        self.widgets['adjustManually'].setEnabled(True)
        self.widgets['errorButton'].setEnabled(True)

    def add_markers(self, data):
        self.selectMarkers.addMarker(self.selectMarkers.model, data)
        self.manual_adjust_marker.addMarker(self.manual_adjust_marker.MRMmodel, data)
        self.error_window.addMarker(self.error_window.EWmodel, data)

    def deleteAllMarkers(self):
        self.selectMarkers.deleteAll()
        self.manual_adjust_marker.deleteAll()
        self.error_window.deleteAll()

    def deleteSelectedMarkers(self,ids, names=[]):
        self.selectMarkers.deleteMarkers(ids, names)
        self.manual_adjust_marker.deleteMarkers(ids)

    def recenter(self, markerFileName='markerfile_ref_TEMP.em', outFileName='markerfile_ref_TEMP.em',
                     tiltSeriesFormat='mrc', return_ref_coords=True, selected_markers=True, save=True):
        if not mf_write:
            print ("NO RECENTERING: loading pytom failed")
            return 0
        markerFileName = os.path.join(os.path.dirname(self.fnames[0]).replace('/excluded', ''), markerFileName)
        ret = self.save(markerFileName=markerFileName, selected_markers=selected_markers)
        if ret == 0: return ret


        # Exclude the excluded frames.
        projIndices = np.arange(len(self.frames))
        take = 1 - np.array(self.excluded, dtype=int)
        projIndices = list(projIndices[take.astype(bool)])
        prefix = '{}/{}/sorted/sorted'.format(self.tomofolder, self.tomogram_name)

        dimBox = int(round((self.fiducial_size / self.pixel_size) *1.5))

        ref_frame = int(self.settings.widgets['ref_frame'].text())
        # Calculate the refined fiducial positions.
        self.ref_coords, self.errors, tX, tY = refineMarkerPositions(prefix, markerFileName, 0,
                                                                     len(self.frames) - 1, outFileName, dimBox=dimBox,
                                                                     projIndices=projIndices,
                                                                     size=self.frames_full[0].shape[0]*self.bin_read,
                                                                     tiltSeriesFormat=tiltSeriesFormat,
                                                                     ret=return_ref_coords, write=False,
                                                                     ireftilt=ref_frame)

        self.ref_coords /= (self.bin_alg* 1.)

        self.old_coords = np.zeros_like(self.coordinates)



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

        model.removeRows(0,model.rowCount())


        for n in range(len(self.errors)):
            error = '{:5.2f}'.format(self.errors[n])
            self.selectMarkers.addMarker(model,[names[n],num_fid[n],error])



        self.update_mark()

        # Save the refined positions if save == True.

        if save:
            for index, imnr in enumerate(projIndices):
                for num, imark in enumerate(markIndices):

                    (xx, yy) = self.ref_coords[index][num]
                    if xx < 0.001 or yy < 0.001: continue
                    for m, (fx, fy) in enumerate(self.mark_frames[imnr]):
                        if fx < 0.1 or fy < 0.1: continue
                        if abs(self.coordinates[imnr][imark][0] - fx) + abs(self.coordinates[imnr][imark][1] - fy) < 0.5:
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

            tiltImageIndex = 0
            for n, (cx, cy) in enumerate( self.mark_frames[tiltNr] ):
                if cx < 0 and cy < 0:
                    continue

                FY,FX = cx*self.bin_alg/self.bin_read, cy*self.bin_alg/self.bin_read
                self.tiltimages[tiltNr].add_fiducial(FX-self.xmin,FY-self.ymin,FX,FY,label=labels[tiltImageIndex])
                
                try: self.tiltimages[tiltNr].indexed_fiducials[tiltImageIndex] = index[tiltImageIndex]
                except: pass
                tiltImageIndex += 1
        self.replot2()

    def save(self, markerFileName='markerfile.txt', ask_recent=1, selected_markers=True):

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
            for iid in range(num_markers):
                name = self.manual_adjust_marker.MRMmodel.data(self.manual_adjust_marker.MRMmodel.index(iid, 0))
                imark = int(name.split('_')[-1])
                markIndices.append(imark)

        if num_markers == 0:
            print ('No markerfile selected', "No markers selected, no markerfile saved.")
            return num_markers

        projIndices = np.arange(len(self.frames))
        take = 1 - np.array(self.excluded, dtype=int)
        projIndices = list(projIndices[take.astype(bool)])
        locX, locY = 1, 2

        if output_type == 'mrc':
            markerFileVol = -1. * np.ones((num_markers, len(projIndices), 12), dtype='float64')
            for num, imark in enumerate(markIndices):
                for (itilt, iproj) in enumerate(projIndices):
                    markerFileVol[num][itilt][0] = self.tiltangles[iproj]
                    markerFileVol[num][itilt][locX] = int(round(self.coordinates[iproj][imark][1] * self.bin_alg))
                    markerFileVol[num][itilt][locY] = int(round(self.coordinates[iproj][imark][0] * self.bin_alg))
            convert_numpy_array3d_mrc(markerFileVol, markerFileName)
        elif output_type == 'em':
            
            markerFileVol = vol(12, len(projIndices), num_markers)
            markerFileVol.setAll(-1)

            for (imark, Marker) in enumerate(markIndices):
                for (itilt, TiltIndex) in enumerate(projIndices):
                    if self.coordinates[TiltIndex][Marker][1] < 1 and self.coordinates[TiltIndex][Marker][0] < 1:
                        markerFileVol.setV(int(round(self.tiltangles[TiltIndex])), 0, itilt, imark)
                        continue

                    markerFileVol.setV(int(round(self.tiltangles[TiltIndex])), 0, itilt, imark)
                    markerFileVol.setV(int(round(self.coordinates[TiltIndex][Marker][1] * self.bin_alg)), locX,
                                       itilt, imark)
                    markerFileVol.setV(int(round(self.coordinates[TiltIndex][Marker][0] * self.bin_alg)), locY,
                                       itilt, imark)

            write_em(markerFileName, markerFileVol)

        elif output_type == 'txt':
            from pytom.basic.datatypes import HEADER_MARKERFILE, FMT_MARKERFILE as fmtMarkerfile
            markerFile = np.ones((num_markers,len(projIndices),4))*-1
            for iMark, Marker in enumerate(markIndices):
                markerFile[iMark,:,0] = iMark

                for (itilt, TiltIndex) in enumerate(projIndices):
                    markerFile[iMark, itilt, 1] = self.tiltangles[TiltIndex]
                    if self.coordinates[TiltIndex][Marker][1] < 1 and self.coordinates[TiltIndex][Marker][0] < 1:
                        continue
                    markerFile[iMark, itilt, 2:] = self.coordinates[TiltIndex][Marker][::-1] * self.bin_alg

            with open(markerFileName, 'w') as outfile:
                np.savetxt(outfile,[],header=HEADER_MARKERFILE)

                for data_slice in markerFile:
                    np.savetxt(outfile, data_slice, fmt=fmtMarkerfile)

        print ('output_type: ', output_type, markerFileName)

    def load(self):
        try:
            markfilename = QFileDialog.getOpenFileName(self, 'Open file', self.projectname, "Marker files (*.*)")
            markfilename = markfilename[0]

            if not markfilename or not os.path.exists(markfilename): return

            mark_frames = guiFunctions.read_markerfile(markfilename,self.tiltangles)

            if markfilename.endswith('.wimp'):
                imodShiftFile = QFileDialog.getOpenFileName(self, 'Open file', self.projectname, "Imod transf file (*.xf)")[0]
                if not imodShiftFile: return
                shifts = guiFunctions.parseImodShiftFile(imodShiftFile)
                mark_frames, shift = guiFunctions.addShiftToMarkFrames(mark_frames, shifts, self.metadata,
                                                                       self.excluded)

            self.deleteAllMarkers()

            if type(mark_frames) in (int, float):
                print('loading file failed')
                return

            self.mark_frames = mark_frames / float(self.bin_read)
            self.coordinates = copy.deepcopy(self.mark_frames)
            self.fs = np.zeros( (len(self.fnames),2) )
            self.frame_shifts = np.zeros( (len(self.fnames),2) )

            for n in range(len(self.tiltimages)):
                self.tiltimages[n].clear()

            itilt, ifid, dummy = self.mark_frames.shape

            for imnr in range(itilt):
                for index in range(ifid):
                    CX,CY = self.mark_frames[imnr][index]

                    if CX < 0 or CY < 0: continue

                    self.tiltimages[imnr].add_fiducial(CX - self.xmin, CY - self.ymin, CX, CY, check=False, draw=False)

            self.widgets['detectButton'].setEnabled(True)
            #self.detect_frameshift()
            if not markfilename.endswith('.npy'):
                self.detect_frameshift()
                self.index_fid()
        except Exception as e:
            print(e)


class SettingsFiducialAssignment(QMainWindow, CommonFunctions):
    def __init__(self,parent=None):
        super(SettingsFiducialAssignment, self).__init__(parent)
        self.setGeometry(900, 0, 300, 100)
        self.layout = self.grid = QGridLayout(self)
        self.setWindowTitle('Settings')
        self.settingsWindow = QWidget(self)
        self.settingsWindow.setLayout(self.layout)
        self.row, self.column = 0, 0
        self.logbook = {}
        self.widgets = {}
        rows, columns = 20, 20
        self.refscores = { 'LoG': 7, 'normal': 1.75, 'sensitive':1.75, 'watershed': 2}



        self.items = [['', ] * columns, ] * rows

        self.insert_label_combobox(self.grid, 'Tomogram Name', 'tomogram_name', self.parent().tomogram_names, rstep=1, cstep=-1,
                                   tooltip='Algorithm used for automatic fiducial detection.')

        self.insert_label_spinbox(self.grid,'bin_read', text='Binning Factor Reading', rstep=1, minimum=1, maximum=16,
                                  stepsize=2,tooltip='Binning factor for reading.',value=4,wtype=QSpinBox,cstep=-1)

        self.insert_label_spinbox(self.grid,'bin_alg', text='Binning Factor Finding Fiducials',rstep=1,
                                  minimum=1,maximum=16,stepsize=2,value=12,wtype=QSpinBox,cstep=0,
                                  tooltip='Binning factor for finding fiducials, used to improve contrast.\n'
                                          'Must be a multiple of the binning factor for reading.')

        self.insert_pushbutton(self.grid, text='Load Tilt Images', rstep=1, action=self.parent().load_images,params='',
                               tooltip='Load tilt images of tomogram set in settings.', cstep=-1)

        self.insert_label(self.grid,rstep=1)

        self.insert_label_spinbox(self.grid, 'tilt_axis', text='Angle Tilt Axis (degrees)', rstep=1,
                               value=270, maximum=359, stepsize=5,
                               tooltip='Specifies how much the tilt axis deviates from vertical (Y axis), clockwise.')

        self.insert_label_spinbox(self.grid, 'ref_frame', text='Reference Frame', rstep=1,
                                  value=19,minimum=1,maximum=91,stepsize=1,
                                  tooltip='Which marker is the reference marker.')
        self.insert_label_spinbox(self.grid, 'ref_marker', text='Reference Marker', rstep=1,
                                  value=0, minimum=0, maximum=91, stepsize=1,
                                  tooltip='Indexing and saving the markerfile uses the reference marker.')

        self.insert_label(self.grid,rstep=1)

        self.insert_label_spinbox(self.grid, 'pixel_size', text=r'Pixel Size (Å)', rstep=1,
                                  value=1.75,maximum=300, stepsize=1., wtype=QDoubleSpinBox,
                                  tooltip=r'Pixel Size in Ångstrom.')

        self.insert_label_spinbox(self.grid, 'fiducial_size', text='Fiducial Size (Å)', rstep=1,
                                  value=50, maximum=350, stepsize=10, minimum=20,
                                  tooltip='Size of Fiducuials in nanometer.')

        self.insert_label(self.grid,rstep=1)

        self.insert_label_checkbox(self.grid, 'applyToAll', 'Search All Frames',
                                   tooltip='Set algorithm to find fiducials in all frames.\nUnchecking means only the '
                                           'current frame is updated.')
        self.insert_label_combobox(self.grid, 'Algorithm', 'algorithm', self.refscores.keys(), rstep=1,cstep=-1,
                                   tooltip='Algorithm used for automatic fiducial detection.')

        self.insert_label_spinbox(self.grid, 'threshold', text='Threshold', rstep=1,
                                  value=7,maximum=10, stepsize=0.1, wtype=QDoubleSpinBox,
                                  tooltip='Threshold detecting a peak in  correlation map.')
        self.insert_label_spinbox(self.grid, 'radius_watershed', text='Radius', rstep=1,
                                  value=3, minimum=1, maximum=100, stepsize=0.5, wtype=QDoubleSpinBox,
                                  tooltip='Radius of particle')



        self.insert_label(self.grid,rstep=1, cstep=1)


        self.widgets['applyToAll'].setChecked(True)
        self.setCentralWidget(self.settingsWindow)

        #CONNECT
        self.widgets['tilt_axis'].valueChanged.connect(self.update_tilt_axis)
        self.widgets['pixel_size'].valueChanged.connect(self.update_radius)
        self.widgets['fiducial_size'].valueChanged.connect(self.update_radius)
        self.widgets['ref_frame'].valueChanged.connect(self.update_ref_frame)
        self.widgets['threshold'].valueChanged.connect(self.update_refscores)
        self.widgets['algorithm'].currentIndexChanged.connect(self.update_baseVValue)
        self.update_radius()

        self.widgets['bin_read'].valueChanged.connect(self.keepMultiple)
        self.widgets['bin_alg'].valueChanged.connect(self.keepMultiple)

    def keepMultiple(self):
        v = self.widgets['bin_read'].value()
        w = self.widgets['bin_alg'].value()

        if v % w != 0:
            w = np.around(((w+v-1)//v)*v,0)

        self.widgets['bin_read'].setValue(v)
        self.widgets['bin_alg'].setValue(w)

    def update_baseVValue(self):
        self.widgets['threshold'].setValue(self.refscores[self.widgets['algorithm'].currentText()])

    def update_refscores(self):
        self.refscores[self.widgets['algorithm'].currentText()] = self.widgets['threshold'].value()

    def update_ref_frame(self):
        self.parent().ref_frame = int(self.widgets['ref_frame'].text())

    def update_tilt_axis(self):
        tilt_axis = float(self.widgets['tilt_axis'].text())
        metafile = self.parent().metafile
        if metafile:
            metadata = self.parent().metadata
            metadata['InPlaneRotation'] = tilt_axis
            guiFunctions.savestar(metafile, metadata, fmt=guiFunctions.fmt, header=guiFunctions.headerText)

    def update_radius(self):
        w = self.widgets
        fiducial_size,pixel_size,bin_read = map(float,(w['fiducial_size'].text(),
                                                       w['pixel_size'].text(),
                                                       w['bin_read'].text()))

        self.parent().radius = fiducial_size/(pixel_size*bin_read*2.)
        w['radius_watershed'].setValue(fiducial_size/pixel_size/w['bin_alg'].value())
        metafile = self.parent().metafile
        if metafile:
            metadata = self.parent().metadata
            metadata['PixelSpacing'] = pixel_size
            metadata['MarkerDiameter'] = fiducial_size
            guiFunctions.savestar(metafile, metadata, fmt=guiFunctions.fmt, header=guiFunctions.headerText)

        if self.parent().loaded_data: self.parent().replot2()


class SelectAndSaveMarkers(QMainWindow,CommonFunctions):
    def __init__(self,parent=None):
        super(SelectAndSaveMarkers,self).__init__(parent)
        self.setGeometry(570, 438, 700, 225)
        self.layout = self.grid = QGridLayout(self)
        self.general = QWidget(self)
        self.general.setLayout(self.layout)
        self.setWindowTitle('Select and Save Marker Sets')
        self.row, self.column = 4, 1
        self.logbook = {}
        self.widgets = {}
        rows, columns = 15,5
        self.header_names = ['Name','#Markers','Rscore']
        self.dtypes = [str, str, float]
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

            model.setData(model.index(model.rowCount()-1, i), self.dtypes[i](data[i]))
            for view in (self.dataView,self.selectView):
            #for i in range(self.num_columns):
                view.resizeColumnToContents(i)

    def deleteMarkers(self,ids, names):
        for m in (self.model,self.model1):
            for row in range(m.rowCount()):
                if self.dataView.model().item(row, 0).text() in names:
                    m.removeRows(row,row+2)

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
        self.setWindowTitle('Manually Adjust Marker Sets')
        self.header_names = ['Name', '#Markers', 'RScore']
        self.dtypes = [str, str, float]
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
        self.MRMmodel = self.model = self.createMarkResultModel(self)
        self.dataView.setModel(self.MRMmodel)

        self.layout.addWidget(self.dataGroupBox, 0, 0, 20, 1)

        self.insert_label(self.grid, rstep=1, cstep=0)
        self.insert_label(self.grid, rstep=1, cstep=0)
        #self.insert_pushbutton(self.grid, text='Add Marker', rstep=1, action=self.add_marker, cstep=0)
        #self.insert_pushbutton(self.grid, text='Delete Marker', rstep=1, action=self.delete_marker, cstep=0)
        #self.insert_label(self.grid, rstep=1, cstep=0)
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

    def deleteMarkers(self,ids):
        for id in reversed(ids):
            print(id)
            self.MRMmodel.removeRows(id,id+1)

    def deleteAll(self):
        self.MRMmodel.removeRows(0,self.MRMmodel.rowCount())

    def add_marker(self, params=None):

        self.parent().add_markers(['Marker_{:03d}'.format(self.MRMmodel.rowCount()), ''])

    def delete_marker(self, params=None):
        #self.MRMmodel.removeRows()
        #for d in dir(self.dataView.selectionModel()): print(d)
        pass
        '''
        idsx = self.dataView.selectionModel().selectedIndexes()
        ids = [id.row() for id in idsx[::2]]

        markernames = [self.dataView.model().item(id.row(),id.column()).text() for id in idsx[::2]]

        if len(ids):
            print(markernames)
            self.parent().deleteSelectedMarkers(ids, markernames)
        '''

    def select_marker(self, params=None):
        try:
            print(self.dataView.selectedIndexes()[0].data()    , ' selected')
            self.parent().selected_marker = int( self.dataView.selectedIndexes()[0].data().split("_")[-1] )
        except: pass

    def deselect_marker(self, params=None):
        if self.parent().selected_marker > -1:
            print ('Marker_{:03d} deselected.'.format(self.parent().selected_marker))
            self.parent().selected_marker = -1
        else:
            pass

    def next_missing(self, params=None):
        try: sel_marker = self.parent().selected_marker
        except: return
        sizeCut = self.parent().sizeCut
        bin_alg = self.parent().bin_alg
        bin_read = self.parent().bin_read
        dim = self.parent().dim
        fs = self.parent().fs
        ref_frame = int(self.parent().settings.widgets['ref_frame'].value())

        if sel_marker > -1:
            imnr = self.parent().imnr
            for i in range(imnr, len(self.parent().fnames)):
                tx,ty = self.parent().coordinates[i, sel_marker]
                if tx < 0.01 or ty < 0.01:

                    found =False
                    for prev in range(1,i+1):
                        px, py = self.parent().coordinates[i-prev, sel_marker]
                        if px > 0 and py > 0:
                            found = True
                            break
                    if not found: break

                    self.parent().imnr = i
                    fx,fy = fs[i-1][0]-fs[i][0], fs[i-1][1]-fs[i][1]
                    print(px*bin_alg/bin_read - fx, py*bin_alg/bin_read-fy, fx, fy)
                    xmin = int(max(0, py*bin_alg/bin_read + fy - sizeCut // 2))
                    ymin = int(max(0, px*bin_alg/bin_read + fx - sizeCut // 2))
                    if xmin + sizeCut >= dim: xmin = dim - sizeCut
                    if ymin + sizeCut >= dim: ymin = dim - sizeCut


                    self.parent().ymin = ymin #int(max(0, px*bin_alg/bin_read + fx - sizeCut // 2))
                    self.parent().xmin = xmin #int(max(0, px*bin_alg/bin_read + fy - sizeCut // 2))
                    self.parent().replot2()
                    break

    def prev_missing(self,params=None):

        try: sel_marker = self.parent().selected_marker
        except: return
        sizeCut = self.parent().sizeCut
        bin_alg = self.parent().bin_alg
        bin_read = self.parent().bin_read
        dim = self.parent().dim
        fs = self.parent().fs
        numim = len(self.parent().fnames)
        if sel_marker > -1:
            imnr = self.parent().imnr
            for i in np.arange(imnr, -1, -1):
                tx, ty = self.parent().coordinates[i, sel_marker]
                if tx < 0.01 or ty < 0.01:
                    found = False
                    for next in range(1, numim-i):
                        px, py = self.parent().coordinates[i + 1, sel_marker]
                        if px > 0 and py >0:
                            found = True
                            break
                    if not found: break

                    self.parent().imnr = i
                    fx, fy = fs[i][0] - fs[i+1][0], fs[i][1] - fs[i+1][1]
                    print(px * bin_alg / bin_read - fx, py * bin_alg / bin_read - fy, fx, fy)
                    xmin = int(max(0, py * bin_alg / bin_read - fy - sizeCut // 2))
                    ymin = int(max(0, px * bin_alg / bin_read - fx - sizeCut // 2))

                    self.parent().ymin = ymin  # int(max(0, px*bin_alg/bin_read + fx - sizeCut // 2))
                    self.parent().xmin = xmin  # int(max(0, px*bin_alg/bin_read + fy - sizeCut // 2))
                    self.parent().replot2()
                    break



        pass



def addLegendGUI(self, offset=(30,30), **kwargs):
    """
    Create a new :class:`~pyqtgraph.LegendItem` and anchor it over the
    internal ViewBox. Plots will be automatically displayed in the legend
    if they are created with the 'name' argument.

    If a LegendItem has already been created using this method, that
    item will be returned rather than creating a new one.

    Accepts the same arguments as :meth:`~pyqtgraph.LegendItem`.
    """

    if self.legend is None:
        self.legend = LegendItemGUI(offset=offset, **kwargs)
        self.legend.setParentItem(self.vb)
    return self.legend

import pyqtgraph as pg

from pyqtgraph.graphicsItems.GraphicsWidget import GraphicsWidget
from pyqtgraph.graphicsItems.LabelItem import LabelItem
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.functions as fn
from pyqtgraph.Point import Point
from pyqtgraph.graphicsItems.GraphicsWidgetAnchor import GraphicsWidgetAnchor



class LegendItemGUI(GraphicsWidget, GraphicsWidgetAnchor):
    """
    Displays a legend used for describing the contents of a plot.
    LegendItems are most commonly created by calling PlotItem.addLegend().

    Note that this item should not be added directly to a PlotItem. Instead,
    Make it a direct descendant of the PlotItem::

        legend.setParentItem(plotItem)

    """

    def __init__(self, size=None, offset=None):
        """
        ==========  ===============================================================
        Arguments
        size        Specifies the fixed size (width, height) of the legend. If
                    this argument is omitted, the legend will autimatically resize
                    to fit its contents.
        offset      Specifies the offset position relative to the legend's parent.
                    Positive values offset from the left or top; negative values
                    offset from the right or bottom. If offset is None, the
                    legend must be anchored manually by calling anchor() or
                    positioned by calling setPos().
        ==========  ===============================================================

        """

        GraphicsWidget.__init__(self)
        GraphicsWidgetAnchor.__init__(self)
        self.setFlag(self.ItemIgnoresTransformations)
        self.layout = QtGui.QGraphicsGridLayout()
        self.setLayout(self.layout)
        self.items = []
        self.size = size
        self.offset = offset
        if size is not None:
            self.setGeometry(QtCore.QRectF(0, 0, self.size[0], self.size[1]))

    def setParentItem(self, p):
        ret = GraphicsWidget.setParentItem(self, p)
        if self.offset is not None:
            offset = Point(self.offset)
            anchorx = 1 if offset[0] <= 0 else 0
            anchory = 1 if offset[1] <= 0 else 0
            anchor = (anchorx, anchory)
            self.anchor(itemPos=anchor, parentPos=anchor, offset=offset)
        return ret

    def addItem(self, item, name):
        """
        Add a new entry to the legend.

        =========== ========================================================
        Arguments
        item        A PlotDataItem from which the line and point style
                    of the item will be determined or an instance of
                    ItemSample (or a subclass), allowing the item display
                    to be customized.
        title       The title to display for this item. Simple HTML allowed.
        =========== ========================================================
        """
        label = LabelItem(name)
        if isinstance(item, ItemSample):
            sample = item
        else:
            sample = ItemSample(item)
        row = len(self.items)
        self.items.append((sample, label))
        self.layout.addItem(sample, row, 0)
        self.layout.addItem(label, row, 1)
        self.updateSize()

    def removeItem(self, name):
        """
        Removes one item from the legend.

        =========== ========================================================
        Arguments
        title       The title displayed for this item.
        =========== ========================================================
        """
        # Thanks, Ulrich!
        # cycle for a match
        for sample, label in self.items:
            if label.text == name:  # hit
                self.items.remove((sample, label))  # remove from itemlist
                self.layout.removeItem(sample)  # remove from layout
                sample.close()  # remove from drawing
                self.layout.removeItem(label)
                label.close()
                self.updateSize()  # redraq box

    def updateSize(self):
        if self.size is not None:
            return

        height = 0
        width = 0
        # print("-------")
        for sample, label in self.items:
            height += max(sample.height(), label.height()) + 3
            width = max(width, sample.width() + label.width())
            # print(width, height)
        # print width, height
        self.setGeometry(0, 0, width + 25, height)

    def boundingRect(self):
        return QtCore.QRectF(0, 0, self.width(), self.height())

    def paint(self, p, *args):
        p.setPen(fn.mkPen(255, 255, 255, 100))
        p.setBrush(fn.mkBrush(100, 100, 100, 50))
        p.drawRect(self.boundingRect())

    def hoverEvent(self, ev):
        ev.acceptDrags(QtCore.Qt.LeftButton)

    def mouseDragEvent(self, ev):
        if ev.button() == QtCore.Qt.LeftButton:
            dpos = ev.pos() - ev.lastPos()
            self.autoAnchor(self.pos() + dpos)


class ItemSample(GraphicsWidget):
    """ Class responsible for drawing a single item in a LegendItem (sans label).

    This may be subclassed to draw custom graphics in a Legend.
    """

    def __init__(self, item):
        GraphicsWidget.__init__(self)
        self.item = item

    def boundingRect(self):
        return QtCore.QRectF(0, 0, 20, 20)

    def paint(self, p, *args):
        # p.setRenderHint(p.Antialiasing)  # only if the data is antialiased.
        opts = self.item.opts

        if opts.get('fillLevel', None) is not None and opts.get('fillBrush', None) is not None:
            p.setBrush(fn.mkBrush(opts['fillBrush']))
            p.setPen(fn.mkPen(None))
            p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(2, 18), QtCore.QPointF(18, 2), QtCore.QPointF(18, 18)]))

        if not isinstance(self.item, pg.ScatterPlotItem):
            p.setPen(fn.mkPen(opts['pen']))
            p.drawLine(2, 18, 18, 2)

        symbol = opts.get('symbol', None)
        if symbol is not None:
            if isinstance(self.item, pg.PlotDataItem):
                opts = self.item.scatter.opts

            pen = pg.mkPen(opts['pen'])
            brush = pg.mkBrush(opts['brush'])
            size = opts['size']

            p.translate(10, 10)
            path = pg.graphicsItems.ScatterPlotItem.drawSymbol(p, symbol, size, pen, brush)


class ErrorWindow(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, logfile=''):
        super(ErrorWindow, self).__init__(parent)
        self.setGeometry(670, 0, 1000, 300)
        self.layout = self.grid = QGridLayout(self)
        self.general = QWidget(self)
        self.general.setLayout(self.layout)
        self.setWindowTitle('Show Alignment Errors')
        self.header_names = ['Name','#Markers', 'RScore']
        self.dtypes = [str, str,float]
        self.logfile = logfile
        self.row, self.column = 0, 1
        self.logbook = {}
        self.widgets = {}
        rows, columns = 20, 20
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
        self.EWmodel = self.model = self.createErrorWindowModel(self)
        self.dataView.setModel(self.EWmodel)

        self.layout.addWidget(self.dataGroupBox, 0, 0, 20, 1)
        self.setCentralWidget(self.general)

        self.view = pg.GraphicsLayoutWidget()

        self.scatterPlot = self.view.addPlot()

        self.scatterPlot.addLegend = addLegendGUI
        self.scatterPlot.addLegend(self.scatterPlot)


        n = 61
        self.s1 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(color='y'), brush=pg.mkBrush(255, 255, 255, 90), identical=True, name='   assigned')
        self.s1.sigClicked.connect(self.mouseHasMoved)
        self.s2 = pg.ScatterPlotItem(size=10, pen=pg.mkPen(color='r'), brush=pg.mkBrush(255, 255, 255, 90), identical=True, name='   unassigned')
        self.s2.sigClicked.connect(self.mouseHasMoved)
        import numpy as np
        pos = np.zeros((2,n),dtype=np.float32)
        pos[0,:] = range(-60,61,2)
        spots = [{'pos': pos[:, i], 'data': 1} for i in range(n)]
        self.s1.addPoints(spots)
        self.scatterPlot.addItem(self.s1)
        #self.s2.addPoints(spots)
        self.scatterPlot.addItem(self.s2)
        self.layout.addWidget(self.view, 0, 1, 20, 15)

        self.scatterPlot.legend.addItem(self.s1, '   assigned')
        self.scatterPlot.legend.addItem(self.s2, '   unassigned')

        self.scatterPlot.setLabels(
            bottom='Tilt Angle (deg)',
            left='Error (pixel)')
        self.scatterPlot.setTitle('Error between expected location fiducial vs assigned location, per tilt image.')
        self.replot()
        self.dataView.selectionModel().selectionChanged.connect(self.replot)

    def mouseHasMoved(self, plot, points):

        try:
            d = points[0].viewPos()
            x, y = d.x(), d.y()

            incl = 1 - np.array(self.parent().excluded)

            itiltFull = np.abs(self.parent().tilt_angles - float(x)).argmin()
            id = np.abs(self.parent().tilt_angles[incl > 0.5] - float(x)).argmin()
            ireftilt = np.abs(self.parent().tilt_angles[incl > 0.5]).argmin()
            itilt = int(np.around(id))
            smallest_difference_angle = np.abs(self.parent().tilt_angles[incl > 0.5] - float(x)).min()
            imark = int(self.dataView.selectedIndexes()[0].data().split("_")[-1])

            bin_read, bin_alg = self.parent().bin_read, self.parent().bin_alg
            # cx, cy = self.parent().coordinates[itiltFull][imark]

            if 1:
                psi, shape = self.parent().psiindeg[imark], read_size(self.parent().fnames[itiltFull])
                shift0 = [self.parent().shiftXModel[ireftilt], self.parent().shiftYModel[ireftilt]]
                shift1 = [self.parent().shiftXModel[itilt], self.parent().shiftYModel[itilt]]

                ccx, ccy = projectMarkerToFrame(self.parent().centersModel[imark], x, psi, shift0, shift1)
                ccx = (ccx + shape[0]//2) / bin_read
                ccy = (ccy + shape[1]//2) / bin_read

            # diffX, diffY = self.parent().diffXModel[itilt][imark], self.parent().diffYModel[itilt][imark]
            # shiftX, shiftY = self.parent().shiftXModel[itilt], self.parent().shiftYModel[itilt]
            sizeCut, (sizeX, sizeY) = self.parent().sizeCut, self.parent().frames_full[0,:,:].shape
            # gx = int(np.around(cx*bin_alg/bin_read + (diffX-shiftX)*bin_read/bin_alg))
            # gy = int(np.around(cy*bin_alg/bin_read + (diffY-shiftY)*bin_read/bin_alg))



            if smallest_difference_angle < 0.2:
                self.parent().imnr = itiltFull
                self.parent().xmin = min(max(0, int(round(ccx)) - sizeCut//2),sizeY-sizeCut//2)
                self.parent().ymin = min(max(0, int(round(ccy)) - sizeCut//2),sizeX-sizeCut//2)
                self.parent().replot2()
                self.parent().updatePartTiltImageGuess(ccx, ccy)

        except Exception as e:
            print(e)
            pass

    def replot(self):
        fname = self.logfile
        id = -1
        try:
            print (self.dataView.selectedIndexes()[0].data(), ' selected for error')
            id = int( self.dataView.selectedIndexes()[0].data().split("_")[-1] )
        except Exception as e:
            print(e)
            return
        if not os.path.exists(fname) or id < 0: return

        data = guiFunctions.loadstar(fname,dtype=guiFunctions.ALIGNMENT_ERRORS)

        markID = data['MarkerIndex']
        angles = data['TiltAngle']
        errors = data['AlignmentError']

        pos = np.array(list(zip(list(angles[markID == id]), list(errors[markID==id]))))
        self.s1.clear()
        spots = [{'pos': pos[i, :], 'data': 1} for i in range(len(pos)) if pos[i,1] > -0.001]
        self.s1.addPoints(spots)
        self.s2.clear()
        spots = [{'pos': [pos[i, 0], 0], 'data': 1} for i in range(len(pos)) if pos[i, 1] <= -0.001]
        self.s2.addPoints(spots)

    def createErrorWindowModel(self, parent):
        model = QStandardItemModel(0, len(self.header_names), parent)
        for i in range(len(self.header_names)):
            model.setHeaderData(i, Qt.Horizontal, self.header_names[i])
        return model

    def addMarker(self, model, data):
        self.EWmodel.insertRow(self.EWmodel.rowCount())
        for i in range(min(len(data), self.num_columns)):
            self.EWmodel.setData(model.index(self.EWmodel.rowCount() - 1, i), self.dtypes[i](data[i]))
            self.dataView.resizeColumnToContents(i)

    def deleteMarkers(self, ids):
        for id in reversed(ids):
            print(id)
            self.EWmodel.removeRows(id, id + 1)

    def deleteAll(self):
        self.EWmodel.removeRows(0,self.EWmodel.rowCount())

def main():
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    app.setWindowIcon(QIcon('/Users/gijs/Documents/PostDocUtrecht/GUI/pp.jpg'))
    ex = FiducialAssignment()
    ex.show()
    sys.exit(app.exec_())

if __name__=='__main__':
    main()
