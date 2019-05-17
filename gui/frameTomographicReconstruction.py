import sys
import os
import random
import glob
import numpy
import time
import shutil
import copy

from os.path import dirname, basename
from ftplib import FTP_TLS, FTP
from multiprocessing import Manager, Event, Process

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets

from pytom.gui.guiStyleSheets import *
from pytom.gui.fiducialAssignment import FiducialAssignment
from pytom.gui.guiStructures import *
from pytom.gui.guiFunctions import avail_gpu, sort
from pytom.gui.guiSupportCommands import *
from pytom.basic.files import read
from pytom_numpy import vol2npy
from pytom.gui.mrcOperations import read_mrc
import pytom.gui.guiFunctions as guiFunctions
from pytom.gui.mrcOperations import square_mrc

class TomographReconstruct(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(TomographReconstruct, self).__init__(parent)

        self.stage = 'v02_'
        self.projectname = self.parent().projectname
        self.tomogram_folder = self.parent().tomogram_folder
        self.rawnanographs_folder = self.parent().rawnanographs_folder
        self.motioncor_folder = self.parent().motioncor_folder
        self.widgets = {}

        self.widgets['pytomPath'] = QLineEdit()
        self.widgets['pytomPath'].setText(self.parent().pytompath)

        self.pytompath = self.parent().pytompath
        headers = ["Select Tomograms", "Create Markerfile", "Alignment", 'CTF Determination', 'Reconstruction']
        subheaders = [[], [], ['Individual Alignment', 'Batch Alignment'], [], ['INFR', 'WBP', 'Batch Reconstruction']]
        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}

        self.subprocesses = 10

        self.tabs = {'tab1': self.tab1,
                     'tab2':  self.tab2,
                     'tab31': self.tab31, 'tab32': self.tab32,
                     'tab4': self.tab4,
                     'tab51': self.tab51, 'tab52': self.tab52, 'tab53': self.tab53}

        self.tab_actions = {'tab1':  self.tab1UI,
                            'tab2':  self.tab2UI,
                            'tab31': self.tab31UI, 'tab32': self.tab32UI,
                            'tab4': self.tab4UI,
                            'tab51': self.tab51UI, 'tab52': self.tab52UI, 'tab53': self.tab53UI}

        for i in range(len(headers)):
            t = 'tab{}'.format(i+1)
            empty = 1*(len(subheaders[i]) == 0)
            for j in range(len(subheaders[i])+empty):
                tt = t+str(j+1)*(1-empty)
                if tt in ('tab2', 'tab31', 'tab32', 'tab4', 'tab51', 'tab52'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()
                button = QPushButton('Refresh Tab')
                button.setSizePolicy(self.sizePolicyC)
                button.clicked.connect(self.tab_actions[tt])

                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)
                self.checkbox[tt] = QCheckBox('sbatch')

                if tt in ('tab1','tab32','tab53'):
                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])

                if tt in ('tab2','tab31','tab4', 'tab51','tab52'):
                    self.tab_actions[tt]()


                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])
                #self.tab_actions[tt]()

    def startFidAssignment(self,parent=None):
        self.fidass = FiducialAssignment(self)
        self.fidass.show()
        pass

    def fill_tab(self, id, headers, types, values, sizes, tooltip=[],connect=0):
        try:
            self.tables[id].setParent(None)
            self.pbs[id].setParent(None)
            self.ends[id].setParent(None)
            self.checkbox[id].setParent(None)
        except:
            pass

        self.tables[id] = SimpleTable(headers, types, values, sizes, tooltip=tooltip, connect=connect)
        self.widgets['v02_batch_aligntable{}'.format(id)] = self.tables[id]
        self.pbs[id] = QPushButton('Run')
        self.pbs[id].setSizePolicy(self.sizePolicyC)
        self.ends[id] = QWidget()
        self.ends[id].setSizePolicy(self.sizePolicyA)

        for a in (self.tables[id], self.checkbox[id], self.pbs[id], self.ends[id]):
            self.table_layouts[id].addWidget(a)

    def tab1UI(self):
        id = 'tab1'

        self.filepath_tomodata = {}
        self.filepath_tomodata['Motion Corrected'] = self.motioncor_folder
        self.filepath_tomodata['Raw Nanographs']  = self.rawnanographs_folder

        headers = ["meta file", "create", 'pututive name',"input files", 'Name tomogram', 'redo', 'delete']
        types = ['txt', 'checkbox', 'txt', 'combobox', 'txt', 'checkbox', 'checkbox']
        sizes = [0, 80, 500, 0, 200, 80, 80]

        tooltip = ['Name of meta files. This file indicates which images belong together, as well as the tiltangles.',
                   'Create new tomogram folders.',
                   'Putative name of tomogram',
                   'Names of existing tomogram folder.',
                   'Redo the creation of a tomogram file.',
                   'Select items to delete tomogram folders.']

        processed = sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder)))
        processed_fn = [basename(line) for line in processed]
        unprocessed = numpy.array(sorted(glob.glob('{}/*.meta'.format(self.rawnanographs_folder))))
        unprocessed = unprocessed[[not (basename(u_item) in processed_fn) for u_item in unprocessed]]

        values = []

        for t in processed:
            values.append([t, False, '', ['Raw Nanographs', 'Motion Corrected'], t.split('/')[-3], True, True])
        for t in unprocessed:
            values.append([t, True, '', ['Raw Nanographs','Motion Corrected'], '', False, False])

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, connect=self.update_create_tomoname)
        #for row in range(self.tables[id].table.rowCount()):

        #   if self.tables[id].table.item(row,2):
        #        for d in dir(self.tables[id].table.item(row,2) ):
        #            for d in dir(self.tables[id].table): print d
        #            print self.tables[id].table.item(row,2).property('isChecked')
        #            break
        #self.tables[id].table.cellChanged.connect(lambda d, ID=id: self.update_create_tomoname(ID))

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.create_tomogram_folders(pid, v))

    def update_create_tomoname(self, id):
        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables['tab1'].table
        widgets = self.tables['tab1'].widgets
        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row,1)
            if wname in widgets.keys() and widgets[wname]:
                if widgets[wname].isChecked():
                    widgets['widget_{}_{}'.format(row, 2)].setText('tomogram_{:03d}'.format(n))
                    n+=1
                else:
                    widgets['widget_{}_{}'.format(row, 2)].setText('')

        table.resizeColumnsToContents()

    def create_tomogram_folders(self, id, values):
        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables['tab1'].table
        widgets = self.tables['tab1'].widgets
        procs = []

        jobs = []

        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row, 1)

            # Create
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = widgets['widget_{}_{}'.format(row, 2)].text()
                metafile = values[row][0]
                folder = self.filepath_tomodata[widgets['widget_{}_{}'.format(row, 3)].currentText()]

                jobs.append([self.create_tomodir_instance, (tomofoldername,metafile,folder)])

                #self.create_tomodir_instance(tomofoldername,metafile,folder)

            # Redo
            wname = 'widget_{}_{}'.format(row, 5)
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = widgets['widget_{}_{}'.format(row, 4)].text()
                metafile = os.path.basename(values[row][0])
                folder = self.filepath_tomodata[widgets['widget_{}_{}'.format(row, 3)].currentText()]
                metafile = os.path.join( self.filepath_tomodata[folder], metafile )
                shutil.rmtree( os.path.join(self.tomogram_folder, tomofoldername) )
                jobs.append( [self.create_tomodir_instance, (tomofoldername,metafile,folder)])

            # Delete
            wname = 'widget_{}_{}'.format(row, 6)
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = widgets['widget_{}_{}'.format(row, 4)].text()
                shutil.rmtree( os.path.join(self.tomogram_folder, tomofoldername) )
                # self.create_tomodir_instance(tomofoldername,mdocfile,folder)

        events = []
        self.jobs_nr_create_tomofolders = len(jobs)

        self.run_jobs(jobs)

    def run_jobs(self, jobs):
        procs = []
        for target, args in jobs:
            procs = [proc for proc in procs if not proc.is_alive()]
            while len(procs) >= self.subprocesses:
                time.sleep(0.5)
            proc = Process(target=target, args=args)
            procs.append(proc)
            proc.start()

        self.update_create_tomoname(id)

    def create_tomodir_instance(self, tomofoldername, metafile, folder):
        src = os.path.join(self.tomogram_folder, '.tomoname')
        num_subprocesses = 10
        dst = os.path.join(self.tomogram_folder, tomofoldername)
        if os.path.exists(dst):
            print('{} exists: QUIT'.format(dst))
            return
        shutil.copytree(src, dst)
        meta_dst = os.path.join(dst, 'sorted', os.path.basename(metafile))
        os.system('cp {} {}'.format(metafile, meta_dst) )

        metafile = numpy.loadtxt(metafile, dtype=guiFunctions.datatype)
        tif_files  = metafile['FileName']
        tiltangles = metafile['TiltAngle']

        if len(tif_files) == 0: return

        for n in range(len(tif_files)):
            tif_files[n] = [line.replace('.tif', '.mrc') for line in repr(tif_files[n]).split('\\') if '.' in line][-1]
        name_angle = list(zip(tif_files, tiltangles))

        sort(name_angle, 1)

        procs = []
        num_copied = 0
        for n, (tif, angle) in enumerate(name_angle):
            while len(procs) >= num_subprocesses:
                time.sleep(0.5)
                procs = [proc for proc in procs if proc.is_alive()]

            src_mcor = os.path.join(folder, tif[1:-1])

            dst_mcor = os.path.join(os.path.dirname(meta_dst), 'sorted_{:02d}.mrc'.format(n))

            if os.path.exists(src_mcor):
                # print('test', src_mcor)
                num_copied += 1
                os.system('cp {} {}'.format(src_mcor, dst_mcor))
                #proc = Process(target=square_mrc, args=([dst_mcor]))
                #procs.append(proc)
                #proc.start()
                out = square_mrc(dst_mcor)
        if num_copied < 5:
            shutil.rmtree(dst)

    def tab2UI(self):
        id = 'tab2'
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = self.table_layouts[id]

        last, reftilt = 10, 5

        #self.insert_label(parent, text='Start Assignment', cstep=1, rstep=1, sizepolicy=self.sizePolicyB, width=200)
        self.insert_pushbutton(parent,cstep=1,text='Create Markerfile',action=self.startFidAssignment, params=[parent],
                               wname='startFidAss',rstep=1)
        self.widgets['startFidAss'].setSizePolicy(self.sizePolicyC)
        self.insert_label(parent,sizepolicy=self.sizePolicyA)
        pass

    def tab31UI(self):
        id = 'tab31'
        self.row, self.column = 0,0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        h = mode = 'v02_Align_'

        last, reftilt = 10, 5
        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted Tilt Images', 'v02_Align_FolderSorted',
                                    'Select the folder where the sorted tiltimages are located.\n')
        self.insert_label_spinbox(parent, mode+'FirstIndex', text='First Index', tooltip='Index of first image.',
                                  value=1,minimum=1,stepsize=1)
        self.insert_label_spinbox(parent, mode +'LastIndex', text='Last Index', tooltip='Index of last image.',
                                  value=last, minimum=2,stepsize=1)
        self.insert_label_spinbox(parent, mode +'RefTiltIndex', text='Reference Tilt Image', value=reftilt,minimum=1,
                                  tooltip='Index of reference tilt image. Typically zeros degree image.')
        self.insert_label_spinbox(parent, mode +'RefMarkerIndex', text='Reference Marker', value=1,
                                  tooltip='Index of reference marker. See previous step.')
        self.insert_label_spinbox(parent, mode + 'RotationTiltAxis', text='Angle Tilt Axis (degrees)',
                                  value=0, minimum=0, maximum=359,
                                  tooltip='Angle of the tilt axis (degrees). 0 degrees is facing norther, '+
                                          '90 degrees is facing east.')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor',value=8, cstep=0,
                                  tooltip='Binning factor used for reconstruction')
        #self.insert_label(parent,'Orientation tilt axis', rstep=1, cstep=1,
        #                  tooltip='Orientation of the tiltaxis with respect to the orientation of the camera.'+
        #                 'The point of the arrow indicates the rotation direction, following the left hand rule.')
        self.widgets[mode + 'tomofolder'] = QLineEdit()
        self.widgets[mode + 'tomogramNR'] = QLineEdit()
        self.widgets[mode + 'FolderSorted'].textChanged.connect(lambda dummy, m=mode: self.updateTomoFolder(m))
        self.updateTomoFolder(mode)

        execfilename = [mode + 'tomofolder', 'alignment/alignment.sh']


        paramsSbatch = guiFunctions.createGenericDict(fname='Alignment',folder=execfilename)
        paramsCmd    = [mode + 'tomofolder', self.parent().pytompath, mode + 'FirstIndex', mode + 'LastIndex',
                        mode + 'RefTiltIndex', mode + 'RefMarkerIndex', mode + 'BinningFactor', mode+'RotationTiltAxis',
                        templateAlignment]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=execfilename, paramsSbatch = paramsSbatch,
                                 paramsCmd=paramsCmd, action=self.convert_em, paramsAction=[mode,'alignment','sorted'])
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)

    def tab32UI(self):
        id='tab32'
        headers = ["name tomogram", "align", 'First Index',"Last Index", 'Ref. Image', 'Ref. Marker', 'Exp. Rot. Angle']
        types = ['txt', 'checkbox', 'lineedit', 'lineedit', 'lineedit', 'combobox','lineedit']
        sizes = [0, 80, 0, 0, 0, 0, 0, 0]

        tooltip = ['Names of existing tomogram folders.',
                   'Do alignment.',
                   'First index of tiltimages',
                   'Last index of tiltimages',
                   'Reference image number.',
                   'Redo the creation of a tomogram file.',
                   'Select items to delete tomogram folders.',
                   'Expected Rotation Angle']

        tomodir = 'Juliette/03_Tomographic_Reconstruction'

        markerfiles = sorted(glob.glob('{}/tomogram_*/sorted/markerfile.em'.format(self.tomogram_folder)))
        print('{}/tomogram_*/sorted/markerfile.em'.format(self.tomogram_folder))


        values = []

        for markerfile in markerfiles:
            qmarkerfile = markerfile.replace('markerfile.em','*.meta')
            qsortedfiles = markerfile.replace('markerfile.em','sorted_*.mrc')
            metafile = glob.glob( qmarkerfile )[0]
            #tangs = numpy.abs(numpy.array(sorted([float(line.split()[-1]) for line in
            #                            os.popen('cat {} | grep TiltAngle '.format(metafile)).readlines()])) )
            metadata = numpy.loadtxt(metafile,dtype=guiFunctions.datatype)
            tangs = metadata['TiltAngle']
            sortedfiles = sorted(glob.glob(qsortedfiles))
            last_frame = len(sortedfiles)
            index_zero_angle = 0
            mm = 9999
            for n, sortedfile in enumerate(sortedfiles):
                index_s = int(sortedfile.split('_')[-1].split('.')[0])
                if abs(tangs[index_s]) < mm:
                    mm = abs(tangs[index_s])
                    index_zero_angle = n+1

            d = read(markerfile)
            data = copy.deepcopy( vol2npy(d) )
            if len(data.shape) < 3: continue
            options_reference = list(map(str, range( data.shape[2] ))) + ['all']
            expect = int(float(metadata['InPlaneRotation'][0]))
            values.append( [markerfile.split('/')[-3], True, 1, last_frame, index_zero_angle, options_reference, expect] )

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip)

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.run_multi_align(pid, v))

        '''
        self.W= QWidget()
        checkButton = QCheckBox('sbatch')
        layout=QVBoxLayout(self.W)
        layout.addWidget(checkButton)
        layout.setAlignment(Qt.AlignLeft)

        #self.checkButton.setA.AlignLeft)
        #self.checkButton.setStyleSheet('QChecbox::indicator{background-color: red;}')
        #self.checkButton.setPalette(QtGui.QPalette(QtGui.QColor(255, 0, 0)))
        self.table_layouts[id].addWidget(self.W,3,1)
        '''

    def run_multi_align(self,id,values):
        print('multi_align', id)
        num_procs = 20
        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables[id].table
        widgets = self.tables[id].widgets
        file_tomoname = os.path.join(self.tomogram_folder, '.multi_alignment.txt')
        tomofolder_file = open(file_tomoname, 'w')
        number_tomonames = 0
        num_procs_per_proc = 0
        expected = values[0][6]
        print(expected)
        for row in range(table.rowCount()):
            wname = 'widget_{}_{}'.format(row, 1)

            # Align
            if wname in widgets.keys() and widgets[wname].isChecked():
                tomofoldername = values[row][0]
                refindex = values[row][4]
                markindex =  widgets['widget_{}_{}'.format(row,5)].currentText()
                tomofolder_file.write('{} {} {}\n'.format(tomofoldername,refindex, markindex))
                num_procs_per_proc = max(num_procs_per_proc, len(values[row][5]) - 1)
                number_tomonames += 1
                folder = os.path.join(self.tomogram_folder, tomofoldername)
                os.system('cp {}/sorted/markerfile.em {}/alignment'.format(folder,folder))

        tomofolder_file.close()

        guiFunctions.batch_tilt_alignment( number_tomonames, fnames_tomograms=file_tomoname, num_procs=20, deploy=True,
                                           projectfolder=self.tomogram_folder, num_procs_per_proc=num_procs_per_proc,
                                           tiltseriesname='sorted/sorted', markerfile='alignment/markerfile.em',
                                           targets='alignment', weightingtype=0, queue=self.checkbox[id].isChecked(),
                                           expectedRotationAngle=expected)

    def tab51UI(self):
        id = 'tab51'
        alg = 'INFR'
        mode = 'v02_{}_'.format(alg)
        self.row, self.column = 0,0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]


        last,reftilt = 10, 5
        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted Tilt Images', mode+'FolderSorted',
                                    'Select the folder where the sorted tiltimages are located.\n')
        self.insert_label_spinbox(parent, mode+'FirstIndex', text='First Index', tooltip='Index of first image.',
                                  value=1,minimum=1,stepsize=1)
        self.insert_label_spinbox(parent, mode +'LastIndex', text='Last Index', tooltip='Index of last image.',
                                  value=last, minimum=2,stepsize=1)
        self.insert_label_spinbox(parent, mode +'RefTiltIndex', text='Reference Tilt Image', value=reftilt,minimum=1,
                                  tooltip='Index of reference tilt image. Typically zeros degree image.')
        self.insert_label_spinbox(parent, mode +'RefMarkerIndex', text='Reference Marker', value=1,
                                  tooltip='Index of reference marker. See previous step.')
        self.insert_label_spinbox(parent, mode + 'RotationTiltAxis', text='Angle Tilt Axis (degrees)',
                                  value=0, minimum=0, maximum=359,
                                  tooltip='Angle of the tilt axis (degrees). 0 degrees is facing norther, ' +
                                          '90 degrees is facing east.')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor',value=8, cstep=0,
                                  tooltip='Binning factor used for reconstruction')

        #self.insert_label(parent,'Orientation tilt axis', rstep=1, cstep=1,
        #                  tooltip='Orientation of the tiltaxis with respect to the orientation of the camera.'+
        #                 'The point of the arrow indicates the rotation direction, following the left hand rule.')

        self.widgets[mode + 'tomofolder'] = QLineEdit()
        self.widgets[mode + 'tomogramNR'] = QLineEdit()
        self.widgets[mode + 'FolderSorted'].textChanged.connect(lambda dummy, m=mode: self.updateTomoFolder(m))
        self.updateTomoFolder(mode)

        execfilename = [mode+'tomofolder','reconstruction/INFR/INFR_reconstruction.sh']
        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'ReconstructionINFR'
        paramsSbatch['folder'] = execfilename
        paramsSbatch['partition'] = 'fastq'
        paramsSbatch['time'] = 1
        paramsSbatch['num_jobs_per_node'] = 1

        self.insert_gen_text_exe(parent,mode,jobfield=False,action=self.convert_em,paramsAction=[mode,'reconstruction/INFR','sorted'],
                                 exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=[mode+'tomofolder',self.parent().pytompath,mode+'FirstIndex',mode+'LastIndex',
                                            mode + 'RefTiltIndex',mode + 'RefMarkerIndex',mode+'BinningFactor',
                                            self.parent().pytompath, mode + 'tomogramNR',  mode+'RotationTiltAxis',
                                            templateINFR])

        #self.insert_label_action_label(parent,'Generate command',cstep=-2, sizepolicy=self.sizePolicyB)
        #self.insert_textfield(parent, h+'CommandText', columnspan=3, rstep=1, cstep=2)
        #self.insert_label_action_label(parent,'Execute command', rstep=1)
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)

    def updateTomoFolder(self, mode):

        folderSorted = self.widgets[mode+'FolderSorted'].text()
        print(folderSorted)
        if not folderSorted: return
        t = folderSorted.replace('/sorted','')
        self.widgets[mode+'tomofolder'].setText(t)
        self.widgets[mode+ 'tomogramNR'].setText( os.path.basename(t) )

        files = [line for line in os.listdir(folderSorted) if line.startswith('sorted') and line.endswith('.mrc')]
        lastIndex = len(files)
        self.widgets[mode+'LastIndex'].setValue(lastIndex)
        metafiles = [line for line in os.listdir(folderSorted) if line.endswith('.meta')]

        if len(metafiles) ==1:
            metafile = metafiles[0]

            try:
                metadata = numpy.loadtxt(os.path.join(folderSorted, metafile), dtype=guiFunctions.datatype)
            except:
                metadata = numpy.loadtxt(os.path.join(folderSorted,metafile),dtype=guiFunctions.datatype0)

            angles = metadata['TiltAngle']
            self.widgets['RotationTiltAxis'] = metadata
            for i in range(len(files)):
                if not 'sorted_{:02d}.mrc'.format(i) in files:
                    angles[i]+=10000

            refIndex = 1 + abs(angles).argmin()
            self.widgets[mode+'RefTiltIndex'].setValue(refIndex)

        if 'WBP' in mode: self.updateVoldims(mode)

    def tab52UI(self):
        id = 'tab52'
        alg = 'WBP'
        h =  mode = 'v02_{}_'.format(alg)
        self.row, self.column = 0,0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]


        last,reftilt = 10, 5

        self.widgets[h + 'Voldims'] = QLineEdit()

        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted Tilt Images', h+'FolderSorted',
                                    'Select the folder where the sorted tiltimages are located.\n')
        self.insert_label_spinbox(parent, mode+'FirstIndex', text='First Index', tooltip='Index of first image.',
                                  value=1,minimum=1,stepsize=1)
        self.insert_label_spinbox(parent, mode +'LastIndex', text='Last Index', tooltip='Index of last image.',
                                  value=last, minimum=2,stepsize=1)
        self.insert_label_spinbox(parent, mode +'RefTiltIndex', text='Reference Tilt Image', value=reftilt,minimum=1,
                                  tooltip='Index of reference tilt image. Typically zeros degree image.')
        self.insert_label_spinbox(parent, mode +'RefMarkerIndex', text='Reference Marker', value=1,
                                  tooltip='Index of reference marker. See previous step.')
        self.insert_label_spinbox(parent, mode + 'RotationTiltAxis', text='Angle Tilt Axis (degrees)',
                                  value=0, minimum=0, maximum=359,
                                  tooltip='Angle of the tilt axis (degrees). 0 degrees is facing norther, '+
                                          '90 degrees is facing east.')

        self.insert_label_spinbox(parent, mode + 'WeightingType', text='Weighting Type',
                                  value=1, minimum=-1, maximum=3000, stepsize=1,
                                  tooltip='Select weighting type:\n\t 0: no weighting\n\t-1: analytical weighting'+
                                          '\n\t 1: "exact" weighting')
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor',value=8, minimum=1, cstep=0,
                                  tooltip='Binning factor used for reconstruction')

        self.widgets[h + 'tomofolder'] = QLineEdit()
        self.widgets[h + 'tomogramNR'] = QLineEdit()
        self.widgets[h + 'FolderSorted'].textChanged.connect(lambda dummy, m=mode: self.updateTomoFolder(m))
        self.updateTomoFolder(mode)

        self.widgets[h + 'BinningFactor'].valueChanged.connect(lambda dummy, m=mode: self.updateVoldims(m))
        self.updateVoldims(mode)
        execfilename = [mode + 'tomofolder', 'reconstruction/WBP/WBP_reconstruction.sh']

        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'ReconstructionWBP'
        paramsSbatch[ 'folder' ] = execfilename
        paramsSbatch['partition'] = 'fastq'
        paramsSbatch['time'] = 1
        paramsSbatch['num_jobs_per_node'] = 1

        paramsCmd = [mode + 'tomofolder', self.parent().pytompath, mode + 'FirstIndex', mode + 'LastIndex',
                     mode + 'RefTiltIndex', mode + 'RefMarkerIndex', mode + 'BinningFactor', mode + 'tomogramNR', 'mrc',
                     mode + 'Voldims', mode + 'WeightingType', mode+'RotationTiltAxis', templateWBP]

        self.insert_gen_text_exe(parent,mode, jobfield=False, action=self.convert_em, exefilename=execfilename,
                                 paramsAction=[mode,'reconstruction/WBP','sorted'],paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.table_layouts[id].addWidget(label)

    def updateVoldims(self,mode):
        folderSorted = self.widgets[mode+'FolderSorted'].text()
        if not folderSorted: return
        files = [line for line in os.listdir(folderSorted) if line.startswith('sorted') and line.endswith('.mrc')]
        if not files: return
        files = files[0]
        imdim = read_mrc(os.path.join(folderSorted, files)).shape[0]
        self.widgets[mode+'Voldims'].setText(str(int(float(imdim)/float(self.widgets[mode+'BinningFactor'].text())+.5)))

    def convert_em(self,params):
        mode = params[0]
        directory = self.widgets[mode+'tomofolder'].text()
        output_folder = '{}/{}'.format(directory, params[1])
        prefix = params[2]
        sorted_folder = self.widgets[mode+'FolderSorted'].text()



        if not os.path.exists(f'{output_folder}/temp_files_unweighted'): os.mkdir(f'{output_folder}/temp_files_unweighted')


        if os.path.basename(params[1]) == 'INFR':
            if len([line for line in os.listdir(output_folder) if line.startswith(prefix.split('/')[-1])]):
                os.system('rm {}/sorted*.em'.format(output_folder))
            guiFunctions.conv_mrc2em(self.widgets[mode+'FolderSorted'].text(), output_folder)
            guiFunctions.renumber_gui2pytom(output_folder, prefix.split('/')[-1])
        #print (self.pytompath)
        #print ('{}/bin/pytom rename_renumber.py {} {} {}'.format(self.pytompath, sorted_folder, output_folder, prefix))


        # print '{}/reconstruction/{}/reconstruction.sh'.format(tomofolderpath,p)
        if 1:
            if os.path.exists('{}/{}_reconstruction.sh'.format(output_folder,os.path.basename(output_folder))):
                for i in range(1000):
                    if not os.path.exists('{}/backup'.format(output_folder)): os.mkdir(
                        '{}/backup'.format(output_folder))
                    fname = "{}/backup/reconstruction_{:03d}".format(output_folder, i)
                    if os.path.exists(fname): continue

                    os.mkdir(fname)
                    for f in ('temp_files_binned', 'temp_files_unweighted', 'temp_files_weighted', 'tomogram.em',
                              '{}_reconstruction.sh'.format(os.path.basename(output_folder)), 'markerfile.em'):
                        if os.path.exists('{}/{}'.format(output_folder, f)):
                            os.system('mv {}/{} {}/'.format(output_folder, f, fname))
                        if not '.' in f and not os.path.exists("{}/{}".format(output_folder, f)):
                            os.mkdir("{}/{}".format(output_folder, f))
                    break
        else:
            pass

        try: os.system(f'cp {directory}/sorted/markerfile.em {output_folder}/markerfile.em')
        except: pass

    def tab53UI(self):
        id='tab53'
        headers = ["name tomogram", "INFR", 'WBP', 'First Index', "Last Index", 'Ref. Image', 'Ref. Marker', 'Exp. Rot. Angle',
                   'Bin Factor']
        types = ['txt', 'checkbox', 'checkbox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit']
        sizes = [0, 80, 80, 0, 0, 0, 0, 0, 0]

        tooltip = ['Names of existing tomogram folders.',
                   'Do INFR Reconstruction.',
                   'Do WBP Reconstruction.',
                   'First index of tiltimages',
                   'Last index of tiltimages',
                   'Reference Image number.',
                   'Reference Marker number.',
                   'Expected in-plane Rotation Angle',
                   'Binning factor applied to images.']

        markerfiles = sorted(glob.glob('{}/tomogram_*/sorted/markerfile.em'.format(self.tomogram_folder)))

        values = []

        for markerfile in markerfiles:
            qmarkerfile = markerfile.replace('markerfile.em','*.meta')
            qsortedfiles = markerfile.replace('markerfile.em','sorted_*.mrc')
            metafile = glob.glob( qmarkerfile )[0]

            metadata=numpy.loadtxt(metafile, dtype=guiFunctions.datatype)
            tilt_angles = metadata['TiltAngle']
            index_zero_angle = numpy.argmin(numpy.abs(tilt_angles)) + 1
            files =  os.listdir( os.path.dirname(markerfile) )
            fnames = sorted([fname for fname in files if fname.startswith('sorted') and fname.endswith('mrc') ])

            print(tilt_angles, fnames)
            rot = int(float(metadata['InPlaneRotation'][0]))

            values.append( [markerfile.split('/')[-3], True, True, 1, len(fnames), index_zero_angle, 1, rot, 8] )

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip)
        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.run_multi_reconstruction(pid, v))

    def run_multi_reconstruction(self, id, values):
        print('multi_reconstructions', id)

        n = len(sorted(glob.glob('{}/tomogram_*/sorted/*.meta'.format(self.tomogram_folder))))
        table = self.tables[id].table
        widgets = self.tables[id].widgets
        mode = 'batch_recon_'
        dd = {1:'reconstruction/INFR',2:'reconstruction/WBP'}

        for row in range(table.rowCount()):
            tomofolder = os.path.join(self.tomogram_folder, values[row][0])
            metafile = glob.glob(os.path.join(tomofolder,'sorted/*.meta'))
            expectedRotation = int(float(widgets['widget_{}_{}'.format(row,7)].text()))
            try:
                metadata = numpy.loadtxt(metafile[-1], dtype=guiFunctions.datatype)
            except:
                continue
            sortedFolder = os.path.join(tomofolder, 'sorted')
            self.widgets[mode + 'tomofolder'] = QLineEdit(text=tomofolder)
            self.widgets[mode + 'FolderSorted'] = QLineEdit(text=sortedFolder)

            for i in (1,2):
                widget = 'widget_{}_{}'.format(row, i)
                if widgets[widget].isChecked():

                    params = [mode,dd[i],'sorted']
                    print(params)
                    self.convert_em(params)

                    execfilename = os.path.join(tomofolder, '{}/{}_Reconstruction.sh'.format(dd[i], dd[i].split('/')[-1]))
                    paramsSbatch = guiFunctions.createGenericDict()
                    paramsSbatch['folder'] = os.path.dirname(execfilename)

                    if i == 1:
                        paramsCmd = [tomofolder, self.pytompath, values[row][3], values[row][4],
                                     values[row][5], values[row][6], values[row][8],
                                     self.pytompath, os.path.basename(tomofolder),
                                     expectedRotation]
                        commandText = templateINFR.format(d=paramsCmd)
                    elif i==2:
                        paramsCmd = [tomofolder, self.pytompath, values[row][3], values[row][4], values[row][5],
                                     values[row][6], values[row][8], os.path.basename(tomofolder), 'mrc', '464', '1',
                                     expectedRotation]
                        commandText= templateWBP.format(d=paramsCmd)
                    else:
                        print( 'No Batch Submission' )

                    if self.checkbox[id].isChecked():
                        header = guiFunctions.gen_queue_header(name=paramsSbatch['fname'], folder=paramsSbatch['folder'],
                                                               modules=paramsSbatch['modules'], time=1, num_jobs_per_node=1,
                                                               partition='fastq')
                        commandText = header + commandText

                    params = [execfilename, commandText]
                    self.submit_multi_recon_job(params)

    def submit_multi_recon_job(self, params):
        print(params[1])
        print(params[0])
        try:

            exefile = open(params[0], 'w')
            exefile.write(params[1])
            exefile.close()

            if len(params[1].split('SBATCH')) > 2:
                os.system('sbatch {}'.format(params[0]))
            else:
                os.system('sh {}'.format(params[0]))
        except:
            print ('Please check your input parameters. They might be incomplete.')

    def tab4UI(self):
        key = 'tab4'
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'CtfDetermination_', self.stage + 'CtfCorrection_'

        items += list(self.create_expandable_group(self.ctfDetermination, self.sizePolicyB, 'CTF Determination',
                                                   mode=t0))
        items[-1].setVisible(False)
        items += list(self.create_expandable_group(self.ctfCorrection, self.sizePolicyB, 'CTF Correction',
                                                   mode=t1))
        items[-1].setVisible(False)
        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def ctfDetermination(self, mode):

        title = "CTF Determination"
        tooltip = 'Run IMOD ctfplotter.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent,cstep=1,sizepolicy=self.sizePolicyB,width=400 )
        self.insert_label_line_push(parent,'Folder Sorted & Aligned Tilt Images', mode + 'FolderSortedAligned',
                                    'Select the folder where the sorted tiltimages are located.\n',
                                    initdir=self.tomogram_folder, mode='folder')
        self.insert_label_line_push(parent, 'Angle File', mode + 'AngleFile',mode='file',filetype='tlt',
                                    tooltip='File with a list of tilt angles. (.tlt extension)',
                                    initdir=self.tomogram_folder)
        self.insert_label_line_push(parent, 'Defocus File', mode + 'DefocusFile',mode='file',filetype='defocus',
                                    tooltip='Filename to store results (.defocus extension).',
                                    initdir=self.tomogram_folder)
        self.insert_label_line_push(parent, 'Config File (OPTIONAL)', mode + 'ConfigFile', mode='file',filetype='cfg',
                                    tooltip='File with a list of noise files used to estimate the noise (.cfg extension).',
                                    initdir=self.tomogram_folder)

        self.insert_label_spinbox(parent, mode +'AxisAngle', text='Axis Angle (deg)',
                                  tooltip='Specifies how much the tilt axis deviates from vertical (Y axis).',
                                  value=270, minimum=0, maximum=359, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'PixelSpacing', text='Pixel Spacing (nm)',
                                  tooltip='Size of pixels in nm.', wtype=QDoubleSpinBox, decimals=3,
                                  value=0.175, minimum=.1, maximum=10, stepsize=.1)
        self.insert_label_spinbox(parent, mode + 'ExpectedDefocus', text='Expected Defocus (nm)',
                                  value=6000, minimum=-60000, maximum=60000, wtype=QDoubleSpinBox,
                                  stepsize=100, tooltip='Expected defocus at the tilt axis in nanometers.')
        self.insert_label_spinbox(parent, mode + 'AngleRangeMin', text='Angle Range(min)',
                                  value=-60, minimum=-180, maximum=-2, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'AngleRangeMax', text='Angle Range (max)',
                                  tooltip='Maximal angle of tilt images.',
                                  value=60, minimum=2, maximum=180, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'Voltage', text='Voltage (kV)',
                                  tooltip='Specifies how much the tilt axis deviates from vertical (Y axis).',
                                  value=200, minimum=2, maximum=1300, stepsize=1)
        self.insert_label_spinbox(parent, mode + 'SphericalAberration', text='Spherical Aberration',
                                  value=2.7, stepsize=0.1, minimum=0, maximum=10., wtype=QDoubleSpinBox,
                                  tooltip='Expected defocus at the tilt axis in nanometers.')
        self.insert_label_spinbox(parent, mode + 'AmplitudeContrast', text='Amplitude Contrast',
                                  value=0.8, stepsize=0.1, minimum=0, maximum=1, wtype=QDoubleSpinBox)
        self.insert_label_checkbox(parent, mode+'FindAstigmatism', 'Find Astigmatism')
        self.insert_label_checkbox(parent, mode + 'FindPhaseShift', 'Find PhaseShift')
        self.insert_label_checkbox(parent, mode + 'FindCutOn', 'Find Cut On Frequency',rstep=1, cstep=0)

        exefilename = [mode+'tomofolder', 'ctf/CTFPlotter.sh']
        xmlfilename = [mode+'tomofolder', 'ctf/ctfplotter.com']

        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'CTF Plotting'
        paramsSbatch[ 'folder' ] = exefilename

        paramsXML = [mode + 'StackFile', mode + 'AngleFile', mode + 'DefocusFile', mode + 'ConfigFile',
                     mode + 'AxisAngle', mode + 'PixelSpacing', mode + 'ExpectedDefocus', mode + 'AngleRangeMin',
                     mode + 'AngleRangeMax', mode + 'Voltage', mode + 'SphericalAberration', mode + 'AmplitudeContrast',
                     mode + 'FindAstigmatism', mode + 'FindPhaseShift', mode + 'FindCutOn', mode + 'ConfigHeader',
                     ParamsFileCTFPlotter]

        paramsCmd = [mode + 'tomofolder', 'ctf/ctfplotter.com', templateCTFPlotter]

        self.insert_gen_text_exe(parent, mode, jobfield=True, exefilename=exefilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, paramsXML=paramsXML, xmlfilename=xmlfilename)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        parent.addWidget(label)

        self.widgets[mode+'tomofolder'] = QLineEdit()
        self.widgets[mode+'StackFile'] = QLineEdit()
        self.widgets[mode+'ConfigHeader'] = QLineEdit()
        self.widgets[mode + 'FolderSortedAligned'].textChanged.connect(lambda dummy, m=mode: self.updateCTFPlotter(m))
        self.widgets[mode + 'ConfigFile'].textChanged.connect(lambda  dummy, m=mode: self.updateConfigHeader(m))
        self.updateCTFPlotter(mode,newstack=False)
        self.updateConfigHeader(mode)

        setattr(self, mode + 'gb_CD', groupbox)
        return groupbox

    def updateConfigHeader(self, mode):
        if self.widgets[mode + 'ConfigFile'].text():
            self.widgets[mode + 'ConfigHeader'].setText('ConfigFile')
        else:
            self.widgets[mode + 'ConfigHeader'].setText('          ')

    def updateCTFPlotter(self, mode, newstack=True):
        folder = self.widgets[mode + 'FolderSortedAligned'].text()
        if not folder:
            print('No Folder Selected.')
            return
        if 1:
            tomogramID = folder.split('tomogram_')[-1][:3]
            sortedFolder = '{}/tomogram_{}/sorted/'.format(self.tomogram_folder, tomogramID)
            tomoname = '{}/tomogram_{}'.format(self.tomogram_folder, tomogramID)

            outstack = '{}/tomogram_{}_{}.st'.format(folder,tomogramID, os.path.basename(folder))
        else:
            print('update CTF Plotter failed.')
            return
        files = [line for line  in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted_')]
        if not files:
            print('No files wuth MRC-format found in: {}'.format(folder))
            return
        dd = []
        for file in files:
            id = int(file.split('_')[-1].split('.')[0])
            dd.append([file,id])
        guiFunctions.sort(dd,1)
        files, ids = zip(*dd)
        if newstack:
            infile, outfile = open(os.path.join(folder, 'filein.txt'), 'w'), open(os.path.join(folder, 'fileout.txt'),'w')
            cmd = 'cd {}; newstack -filei filein.txt -fileo fileout.txt '.format(folder)
            outfile.write('{}\n{}\n{}\n'.format(1, outstack, len(files)))
            infile.write('{}\n'.format(len(files)))
            for fname in files:
                infile.write('{}\n0\n'.format(fname))
            infile.close()
            outfile.close()
            os.system(cmd)

        metafile = [os.path.join(sortedFolder, line) for line in os.listdir(sortedFolder) if line.endswith('.meta')][0]
        metadata = numpy.loadtxt(metafile,dtype=guiFunctions.datatype)
        outAngle   = outstack.replace('.st','.tlt')
        out = open(outAngle,'w')
        for id in ids: out.write('{}\n'.format(metadata['TiltAngle'][id]))
        out.close()
        outDefocus = os.path.join(tomoname, 'ctf', os.path.basename(outstack).replace('.st','.defocus'))

        self.widgets[mode + 'tomofolder'].setText(tomoname)
        self.widgets[mode + 'StackFile'].setText(outstack)
        self.widgets[mode + 'AngleFile'].setText(outAngle)
        self.widgets[mode + 'DefocusFile'].setText(outDefocus)

        self.widgets[mode + 'AxisAngle'].setValue(metadata['InPlaneRotation'][0])
        self.widgets[mode + 'PixelSpacing'].setValue(metadata['PixelSpacing'][0]/10.)
        self.widgets[mode + 'ExpectedDefocus'].setValue(1000.*metadata['DefocusU'][abs(metadata['TiltAngle']).argmin()])
        self.widgets[mode + 'AngleRangeMin'].setValue(numpy.floor(metadata['TiltAngle'][ids[0]]))
        self.widgets[mode + 'AngleRangeMax'].setValue(numpy.ceil(metadata['TiltAngle'][ids[-1]]))
        self.widgets[mode + 'Voltage'].setValue(metadata['Voltage'][0])
        self.widgets[mode + 'SphericalAberration'].setValue(metadata['SphericalAberration'][0])
        self.widgets[mode + 'AmplitudeContrast'].setValue(metadata['AmplitudeContrast'][0])

        print(self.widgets[mode + 'tomofolder'].text())

    def prep_value(self, params):

        mode = params[0]
        outtext = '''# command file to run ctfplotter\n#\n'''

        folder = self.widgets[mode +'FolderSortedAligned'].text()
        metafile = [line for line in os.listdir('{}/../../sorted/'.format(folder)) if line.endswith('.meta')][0]
        metadata = numpy.loadtxt(metafile,dtype=guiFunctions.datatype)

        files = [os.path.join(folder,line) for line in os.listdir(folder) if line.endswith('.mrc') and line.startswith('sorted_aligned')]

        files = sorted(files)

    def ctfCorrection(self, mode):
        title = "CTF Correction"
        tooltip = 'CTF Correction for aligned tilt images.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB, width=400)

        self.insert_label_line_push(parent, 'Folder Sorted & Aligned Tilt Images', mode + 'FolderSortedAligned',
                                    'Select the folder where the sorted tiltimages are located.\n',
                                    initdir=self.tomogram_folder, mode='folder')
        self.insert_label_line_push(parent, 'Defocus File', mode + 'DefocusFile', mode='file', filetype='defocus',
                                    tooltip='Filename to store results (.defocus extension).',
                                    initdir=self.tomogram_folder)
        self.insert_label_spinbox(parent, mode + 'GridSpacing', text='Grid Spacing',
                                  tooltip= 'grid spacing [2,4,6,...] is the size of the area which is corrected with a'+
                                           'constant ctf. The ctf value is taken from the center of this area. This '+
                                           'area builds up the corrected projection.',
                                  value=1, stepsize=1, minimum=1, maximum=1000)
        self.insert_label_spinbox(parent, mode + 'FieldSize', text='FieldSize',
                                  tooltip='fieldsize [2,4,6,...] & (fs>=gs) is the size of the area which is extracted'+
                                          ' from the projection and corrected with a constant ctf. Fieldsize is also '+
                                          'the size of the modelled ctf.',
                                  value=1, stepsize=1, minimum=1, maximum=1000)
        self.insert_label_spinbox(parent, mode + 'BinningFactor', text='Binning Factor',
                                  value=1, stepsize=1, minimum=1, maximum=16, rstep=1, cstep=0)

        exefilename = [mode + 'tomofolder', 'ctf/CTFCorrection.sh']


        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'CTF_Correction'
        paramsSbatch['folder'] = exefilename


        paramsCmd = [mode + 'tomofolder', self.pytompath, mode+'uPrefix', mode + 'cPrefix', mode + 'MetaFile',
                     mode + 'GridSpacing', mode + 'FieldSize', mode + 'BinningFactor', templateCTFCorrection]

        self.insert_gen_text_exe(parent, mode, exefilename=exefilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, action=self.updateMetaFile, paramsAction=[mode])

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        parent.addWidget(label)

        self.widgets[mode + 'tomofolder'] = QLineEdit()
        self.widgets[mode + 'uPrefix'] = QLineEdit()
        self.widgets[mode + 'cPrefix'] = QLineEdit()
        self.widgets[mode + 'MetaFile'] = QLineEdit()

        self.widgets[mode + 'FolderSortedAligned'].textChanged.connect(lambda dummy, m=mode: self.updateCTFCorrection(m))
        self.widgets[mode + 'GridSpacing'].valueChanged.connect(lambda dummy, m=mode: self.updateGridAndFieldSize(m))

        self.updateCTFCorrection(mode)
        self.updateGridAndFieldSize(mode)

        setattr(self, mode + 'gb_CD', groupbox)
        return groupbox

    def updateCTFCorrection(self, mode):
        folder = self.widgets[mode + 'FolderSortedAligned'].text()
        if not folder:
            print('CTF Correction: empty folder sorted aligned')
            return
        try:
            tomogramID = folder.split('tomogram_')[-1][:3]
            tomoname = '{}/tomogram_{}'.format(self.tomogram_folder, tomogramID)
            cOut = '{}/ctf/{}'.format(tomoname, os.path.basename(folder))
            if not os.path.exists(cOut): os.mkdir(cOut)
            cPrefix = os.path.join(cOut, 'sorted_aligned_ctf_')
            uPrefix = os.path.join(folder, 'sorted_aligned_')
            sortedFolder = os.path.join(tomoname, 'sorted')
            print(sortedFolder)
            metafile = glob.glob(sortedFolder + '/*.meta')[0]
            #defocusfile = glob.glob( os.path.join( os.path.dirname(cOut), '*.defocus'))
        except:
            self.popup_messagebox('Error', 'Setting values for CTF Correction Failed',
                                  'The setting of the parameters for CTF correction has failed. \n' +
                                  'Running ctf correction might result in errors. ' +
                                  'Try to select Folder Sorted Aligned again.')
            return

        self.widgets[mode +'tomofolder'].setText(tomoname)
        self.widgets[mode + 'uPrefix'].setText( uPrefix )
        self.widgets[mode + 'cPrefix'].setText( cPrefix )
        self.widgets[mode + 'MetaFile'].setText( metafile )
        
    def updateGridAndFieldSize(self, mode):
        gridSpacing = float(self.widgets[mode + 'GridSpacing'].value())
        fieldSize   = float(self.widgets[mode + 'FieldSize'].value())
        if gridSpacing > fieldSize:
            self.widgets[mode + 'GridSpacing'].setValue(fieldSize)
        self.widgets[mode+'GridSpacing'].setMinimum(fieldSize)

    def updateMetaFile(self,params):
        mode = params[0]
        metafile = self.widgets[mode + 'MetaFile'].text()
        defocusfile = self.widgets[mode + 'DefocusFile'].text()
        try:
            guiFunctions.update_metadata_from_defocusfile(metafile,defocusfile)
        except:
            self.popup_messagebox('Error','Update MetaData Failed', 'Update metadata has failed. Your job is not using the paramaters from the selected defocus file.')
