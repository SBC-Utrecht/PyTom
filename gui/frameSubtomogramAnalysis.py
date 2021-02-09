import sys
import os
import random
import glob
import numpy
import time

from multiprocessing import Manager, Event, Process
from ftplib import FTP_TLS, FTP
from os.path import dirname, basename

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets

from pytom.basic.files import read
from pytom.gui.guiStyleSheets import *
from pytom.gui.guiSupportCommands import *
from pytom.gui.guiStructures import *
from pytom.gui.fiducialAssignment import FiducialAssignment
from pytom.gui.guiFunctions import avail_gpu
import pytom.gui.guiFunctions as guiFunctions
from pytom.bin.extractTomoNameFromXML import *
from pytom.gui.guiFunctions import readMarkerfile
from pytom.tompy.io import read_size

class SubtomoAnalysis(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(SubtomoAnalysis, self).__init__(parent)
        self.stage          = 'v04_'
        self.pytompath      = self.parent().pytompath
        self.projectname    = self.parent().projectname
        self.logfolder      = self.parent().logfolder
        self.subtomodir     = self.parent().subtomo_folder
        self.tomoanalysis   = self.parent().tomogram_folder
        self.fscdir         = os.path.join(self.subtomodir, 'Validation')
        self.polishfolder   = os.path.join(self.subtomodir, 'ParticlePolishing')
        self.frmdir         = os.path.join(self.subtomodir,'Alignment/FRM')
        self.glocaldir      = os.path.join(self.subtomodir, 'GLocal/FRM')
        self.cpcadir        = os.path.join(self.subtomodir, 'Classification/CPCA')
        self.acdir          = os.path.join(self.subtomodir, 'Classification/AutoFocus')
        self.pickpartdir    = self.parent().particlepick_folder+'/Picked_Particles'
        self.tomogramfolder = os.path.join(self.parent().particlepick_folder, 'Tomograms')
        self.tomogram_folder = self.parent().tomogram_folder
        self.acpath = os.path.join(self.parent().subtomo_folder, 'Classification/AutoFocus')
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand
        self.progressBarCounters = {}
        self.progressBars = {}
        self.queueEvents = self.parent().qEvents
        self.localqID = {}
        self.activeProcesses = {}
        self.threadPool = self.parent().threadPool
        self.workerID = 0
        self.qparams = {}
        self.tabs_dict = {}
        self.tab_actions = {}


        # CHANGE THE NEXT FOUR VARIABLES WHEN ADDING OR REMOVING TABS
        # ONE NEEDS 1) UI FUNCTION TO SET UP THE GENERAL TAB (NUMBER OF DROPDOWN FIELDS VS DYNAMIC TABLE)
        #           2) FUNCTION TO FILL DROP DOWN MENU / FILL THE TABLE
        #           3) HELPER FUNCTIONS TO AUTOFILL OR UPDATE FIELDS BASED ON USER INPUT

        headers = ["Reconstruct Subtomograms", "Particle Polishing", "Align Subtomograms", "Classify Subtomograms", "Validation"]
        subheaders = [['Single Reconstruction','Batch Reconstruction'],['Single', 'Batch'], ['FRM Alignment','GLocal'],['CPCA','Auto Focus'], [] ]
        tabUIs = [[self.SubtomoReconstrSingleUI, self.SubtomoReconstrBatchUI],
                  [self.PolishSingleUI, self.PolishBatchUI],
                  [self.FRMUI,self.GLocalUI],
                  [self.CPCAUI,self.AC3DUI],
                  self.FSCUI]
        static_tabs = [[True, False], [True, False], [True, True], [True, True], [True, True]]




        self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders,tabUIs=tabUIs,tabs=self.tabs_dict, tab_actions=self.tab_actions)

        self.widgets = {}
        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.num_nodes = {}
        self.checkbox = {}

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if static_tabs[i][j]:# tt in ('tab11', 'tab31', 'tab32', 'tab41', 'tab42'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()
                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)
                self.checkbox[tt] = QCheckBox('queue')

                if not static_tabs[i][j]:#tt in ('tab12'):
                    button = QPushButton('Refresh Tab')
                    button.setSizePolicy(self.sizePolicyC)
                    button.clicked.connect(lambda d, k=tt, a=self.tab_actions[tt]: a(k))

                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])
                else:#if not tt in ('tab12'):
                    self.tab_actions[tt](tt)

                tab = self.tabs_dict[tt]
                tab.setLayout(self.table_layouts[tt])


    # General UI functions

    def dd(self, key):
        pass

    def SubtomoReconstrSingleUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleReconstruction_'

        items += list(self.create_expandable_group(self.addSubtomogramReconstructionFields, self.sizePolicyB,
                                                   'Single Reconstruction', mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def SubtomoReconstrBatchUI(self, key=''):
        try: self.extractLists.text()
        except: self.extractLists = QLineEdit()

        self.mass_extract = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.xml'],
                                        outputline=self.extractLists, id=key,
                                        run_upon_complete=self.populateSubtomoReconBatchTable,
                                        title='Select particlLists')

    def PolishSingleUI(self, key=''):

        return

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'SingleParticlePolish_'

        items += list(self.create_expandable_group(self.addParticlePolishFields, self.sizePolicyB, 'Single Reconstruction',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def PolishBatchUI(self, key=''):

        return

        try: self.polishLists.text()
        except: self.polishLists = QLineEdit()

        self.selectPLWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.xml'],
                                          outputline=self.polishLists, id=key,
                                          run_upon_complete=self.selectTemplates,
                                          title='Select particlLists')

    def FRMUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1, t2 =  self.stage + 'inputFiles_', self.stage + 'frmSetttings_',self.stage + 'sampleInformation_'

        items += list(self.create_expandable_group(self.addFRMAlignmentFields, self.sizePolicyB, 'FRM Alignment',
                                                   mode=t0))
        items[-1].setVisible(False)


        for n, item in enumerate(items):
            grid.addWidget(item, n, 0,1,4)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def GLocalUI(self, key=''):


        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        items += list(self.create_expandable_group(self.addGLocalAlignmentFields, self.sizePolicyB, 'GLocal Alignment',
                                                   mode=self.stage + 'gLocal_'))
        items[-1].setVisible(False)
        # items += list( self.create_expandable_group(self.createSubtomograms, self.sizePolicyB, 'Extract Subtomograms',
        #                                            mode=self.stage+'extract_') )
        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def CPCAUI(self, key=''):

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'CCC_', self.stage + 'CPCA_'

        items += list(self.create_expandable_group(self.addCCCFields, self.sizePolicyB, 'Pairwise Cross Correlation',
                                                   mode=t0))
        items[-1].setVisible(False)

        items += list(self.create_expandable_group(self.addCPCAFields, self.sizePolicyB, 'CPCA',
                                                   mode=t1))
        items[-1].setVisible(False)
        # items += list(self.create_expandable_group(self.sampleInformation, self.sizePolicyB, 'Sample Information',
        #                                           mode=t2))
        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def AC3DUI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'AC_', self.stage + 'CPCA_'

        items += list(self.create_expandable_group(self.addAC3DFields, self.sizePolicyB, 'Autofocussed Classification',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def FSCUI(self,key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0 = self.stage + 'FSC_'

        items += list(self.create_expandable_group(self.addFSCFields, self.sizePolicyB, 'Fourier Shell Correlation',
                                                   mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)


    # Functions to fill tabs with input fields or an input table

    def addSubtomogramReconstructionFields(self, mode=''):

        title = "Single Reconstruction"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particlelist',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Folder with aligned tilt images', mode + 'AlignedTiltDir',
                                    'Select the folder with the aligned tilt images.')
        self.insert_label_line_push(parent, 'Meta file with tilt angles', mode + 'MetaFile', mode='file',
                                    filetype='meta',
                                    tooltip='Select the corresponding metafile.')
        self.insert_label_line_push(parent, 'Result file particle polishing', mode + 'polishFile', mode='file',
                                    filetype='txt',initdir=self.polishfolder,
                                    tooltip='Select a resultfile from particle polishing.')
        self.insert_label_spinbox(parent, mode + 'BinFactorReconstruction',
                                  'Binning factor used in the reconstruction.',
                                  'Defines the binning factor used in the reconstruction of the tomogram from which' +
                                  'the particles are selected.',
                                  minimum=1, stepsize=1, value=8)

        self.insert_label_spinbox(parent, mode + 'WeightingFactor', 'Apply Weighting (0/1)',
                                  'Sets the weighting scheme applied to the tilt images.\n' +
                                  '0: no weighting.\n1: ramp filter.', minimum=-5, maximum=5, stepsize=1, value=0)

        self.insert_label_spinbox(parent, mode + 'SizeSubtomos', 'Size subtomograms.',
                                  'Sets the size of the subtomograms.',
                                  minimum=10, maximum=1000, stepsize=1, value=128)

        self.insert_label_spinbox(parent, mode + 'BinFactorSubtomos', 'Binning Factor Subtomograms.',
                                  'Sets the binning factor of the subtomograms.', rstep=1,
                                  value=1, stepsize=1, minimum=1)

        self.insert_label_spinbox(parent, mode + 'OffsetX', 'Offset in x-dimension',
                                  'Has the tomogram been cropped in the x-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the x-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetY', 'Offset in y-dimension',
                                  'Has the tomogram been cropped in the y-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the y-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetZ', 'Offset in z-dimension',
                                  'Has the tomogram been cropped in the z-dimension?\n' +
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping' +
                                  ' in the z-dimension.', cstep=0, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)


        self.widgets[mode + 'polishFlag'] = QLineEdit('')

        self.widgets[mode + 'particlelist'].textChanged.connect(lambda d, m=mode: self.updateMeta(m))
        self.widgets[mode + 'polishFile'].textChanged.connect(lambda d, m=mode: self.updatePolishFlag(m))

        execfilename = os.path.join(self.subtomodir, 'Reconstruction/reconstructSubtomograms.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='subtomoReconstr', folder=self.logfolder,
                                                      id='SingleSubtomoReconstruct')
        paramsCmd = [mode + 'particlelist', mode + 'AlignedTiltDir', mode + 'BinFactorReconstruction',
                     mode + 'SizeSubtomos', mode + 'BinFactorSubtomos', mode + 'OffsetX', mode + 'OffsetY',
                     mode + 'OffsetZ', self.subtomodir, mode + 'WeightingFactor', mode + 'MetaFile', '20',
                     mode + 'polishFlag', extractParticles]

        self.updatePolishFlag(mode)

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename, paramsSbatch=paramsSbatch)
        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def populateSubtomoReconBatchTable(self, key='tab12'):
        self.mass_extract.close()
        particleFilesStart = sorted(self.extractLists.text().split('\n'))
        particleFiles = []

        for particleFile in particleFilesStart:
            if '_tomogram_' in particleFile:
                particleFiles.append(particleFile)
            else:
                output_folder = os.path.join(self.pickpartdir, os.path.splitext(os.path.basename(particleFile))[0])
                pLs = extractParticleListsByTomoNameFromXML(particleFile, directory=output_folder)
                for pl in pLs:
                    particleFiles.append(pl)

        particleFiles = sorted(particleFiles)

        headers = ["Filename particleList", "Run", "Tilt Images", "Alignment Type", "Origin", 'Bin factor recon', 'Weighting',
                   "Size subtomos", "Bin subtomos", "Offset X", "Offset Y", "Offset Z", 'Polish Result', '']
        types = ['txt', 'checkbox', 'combobox', 'combobox', 'combobox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit',
                 'lineedit', 'lineedit', 'comboboxF', 'txt']
        a = 40
        sizes = [0, 0, 80, 80, 90, a, a, a, a, a, a, a, a, 80, a]

        tooltip = ['Names of the particleList files',
                   'Check this box to run subtomogram reconstruction.',
                   'Which folder contains the tilt-images you want to use for subtomo reconstruction?',
                   'Aligned Images',
                   'Binning factor used for the reconstruction.',
                   'Weighting Type.\n0: No Weighting,\n1: Analytical Weighting.\n-1: Ramp Weighting',
                   'Size Subtomograms', 'Binning factor for subtomograms (--projBinning)', 'Offset in X-dimension',
                   'Offset in Y-dimension', 'Offset in Z-dimension', 'Results files from particle polishing.']

        values = []
        refmarkindices = []

        for n, particleFile in enumerate(particleFiles):
            if not particleFile: continue
            base, ext = os.path.splitext(
                os.path.basename(particleFile).replace('particleList_', '').replace('coords_', '').replace('_flipped',
                                                                                                           ''))
            if '_tomogram_' in base:
                base = 'tomogram_' + base.split('_tomogram_')[1]

            for t in ('WBP', 'INFR'):
                if t in base: base = base.split(t)[0] + t

            tomoindex = base.split('tomogram_')[1].split('_')[0]

            polishfiles = glob.glob(f'{self.polishfolder}/resultsPolish*tomogram_{tomoindex}*.txt')
            polishfiles += glob.glob(f'{self.polishfolder}/*/resultsPolish*tomogram_{tomoindex}*.txt')
            polishfiles += ['']


            if base + '.mrc' in os.listdir(self.tomogramfolder) or base + '.em' in os.listdir(self.tomogramfolder):
                if os.path.exists(os.path.join(self.tomogramfolder, base + '.mrc')):
                    folder = os.popen('ls -alrt {}.mrc'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                elif os.path.exists(os.path.join(self.tomogramfolder, base + '.em')):
                    folder = os.popen('ls -alrt {}.em'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                else:
                    folder = ''
                if not folder: continue
                folder = os.path.dirname(folder.split()[-1])

                # markerfile  = os.path.join(folder, 'markerfile.txt')
                # markerdata = readMarkerfile(markerfile,61)

                al = os.path.join(os.path.dirname(os.path.dirname(folder)), 'alignment')
                ctf = os.path.join(os.path.dirname(os.path.dirname(folder)), 'ctf')
                choices = [al + '/' + f for f in os.listdir(al) if 'marker_' in f and os.path.isdir(al + '/' + f)]
                # choices = list(map(str,range(markerdata.sizeZ()))) # + ['closest']
                # a = sorted(glob.glob('{}/Reconstruction*-*.out'.format(folder)))[-1]

                try:
                    from lxml import etree
                    xmlObj = etree.parse(particleFile)
                    particles = xmlObj.xpath('Particle')
                    binning = 8 #int(particles[0].xpath('InfoTomogram')[0].get('BinningFactor'))
                    refmarkindex = int(particles[0].xpath('InfoTomogram')[0].get('RefMarkIndex'))
                except:
                    print('Default values used for {}:\n\tbin recon = 8\n\t ref mark index = 1'.format(
                        os.path.basename(particleFile)))
                    binning = 8
                    refmarkindex = 1
                # binning = os.popen('cat {} | grep "--referenceMarkerIndex" '.format(a)).read()[:-1]
                # print(binning)

                closest_options = []
                for c in choices:
                    print(c)
                    align_folder = os.path.dirname(c)
                    mm, markID, angles = os.path.basename(c).split('_')
                    closestTemp = f'{align_folder}/{mm}_CLOSEST_{angles}'
                    if not closestTemp in closest_options: closest_options.append(closestTemp)

                if choices:
                    aligntype = [f'{choices[0]}/{f}' for f in os.listdir(choices[0]) if os.path.isdir(f'{choices[0]}/{f}')]
                    if aligntype:
                        origin = [f'{aligntype[0]}/{f}' for f in os.listdir(aligntype[0]) if os.path.isdir(f'{aligntype[0]}/{f}')]
                    else:
                        origin = []

                else:
                    aligntype, origin = [],[]

                choices += closest_options

                values.append([particleFile, True, choices, aligntype, origin, binning, -1, 128, 1, 0, 0, 0, polishfiles, ''])
                refmarkindices.append(refmarkindex)

        try:
            self.num_nodes[key].setParent(None)
        except:
            pass

        if values:
            self.valuesBatchSubtomoReconstruction = values
            self.fill_tab(key, headers, types, values, sizes, tooltip=tooltip, nn=True)
            self.tab12_widgets = self.tables[key].widgets
            for n, index in enumerate(refmarkindices):
                self.tab12_widgets['widget_{}_3'.format(n)].setCurrentIndex(index)

            for row in range(len(values)):
                w = self.tab12_widgets['widget_{}_2'.format(row)]
                w.currentIndexChanged.connect(lambda d, r=row, k=key: self.updateChoices(r, k))

                ww = self.tab12_widgets['widget_{}_3'.format(row)]
                ww.currentIndexChanged.connect(lambda d, r=row, k=key: self.updateAlignmentType(r, k))


            self.particleFilesBatchExtract = particleFiles
            self.pbs[key].clicked.connect(lambda dummy, pid=key, v=values: self.massExtractParticles(pid, v))
        else:
            return

        for i in range(len(values)):
            self.updateChoices(i, key)

    def addParticlePolishFields(self, mode=''):
        title = "Particle Polish (Single)"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Folder with aligned tilt images', mode + 'AlignedTiltDir',
                                    'Select the folder with the aligned tilt images.')
        self.insert_label_line_push(parent, 'Template File', mode + 'Template',mode='file',filetype=['em', 'mrc'],
                                    tooltip='Select a template file.')
        self.insert_label_line_push(parent, 'FSC File', mode + 'FSCPath', mode='file', filetype=['dat'],
                                    tooltip='Select a FSC file.')
        self.insert_label_line_push(parent, 'Meta File', mode + 'MetaFile', mode='file', filetype=['meta'],
                                    tooltip='Select the corresponding meta file (in tomogram_???/sorted)')
        self.insert_label_line_push(parent, 'Alignment Results File', mode + 'AlignmentResultsFile', mode='file', filetype=['txt'],
                                    tooltip='Select the alignment results file (in alignment/???/)')
        self.insert_label_line_push(parent, 'Output Directory', mode + 'destination', mode='folder',
                                    tooltip='Select the destination directory.')

        self.insert_label_spinbox(parent,mode+'dimZ', 'Size tomogram in Z-dimension.',
                                  'Sets the size of z-dimension of the origin tomogram.',
                                  minimum=10,maximum=4000,stepsize=1,value=200)
        self.insert_label_spinbox(parent, mode + 'BinFactorReconstruction', 'Binning factor used in the reconstruction.',
                                  'Defines the binning factor used in the reconstruction of the tomogram from which'+
                                  'the particles are selected.',
                                  minimum=1,stepsize=1,value=8)
        self.insert_label_spinbox(parent, mode+'maxParticleShift', 'Max particle shift.',
                                  'Sets the max of the particle shift.',rstep=1,
                                  value=25, stepsize=1, minimum=1)
        self.insert_label_spinbox(parent,  mode + 'OffsetX', 'Offset in x-dimension',
                                  'Has the tomogram been cropped in the x-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the x-dimension.', cstep=-1, rstep=1,
                                  value=0, stepsize=1,minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetY', 'Offset in y-dimension',
                                  'Has the tomogram been cropped in the y-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the y-dimension.', cstep=-1,rstep=1,
                                  value=0, stepsize=1,minimum=-4000, maximum=4000)
        self.insert_label_spinbox(parent, mode + 'OffsetZ', 'Offset in z-dimension',
                                  'Has the tomogram been cropped in the z-dimension?\n'+
                                  'If so, add the cropped magnitude as an offset.\nExample: 200 for 200 px cropping'+
                                  ' in the z-dimension.', cstep=0, rstep=1,
                                  value=0, stepsize=1, minimum=-4000, maximum=4000)

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'FSCFlag'] = QLineEdit('')
        self.widgets[mode + 'MetaFileFlag'] = QLineEdit('')
        self.widgets[mode + 'AlignmentResultsFileFlag'] = QLineEdit('')


        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updateMeta(m))
        self.widgets[mode + 'FSCPath'].textChanged.connect(lambda d, m=mode: self.updateFSCFlag(m))
        self.widgets[mode + 'MetaFile'].textChanged.connect(lambda d, m=mode: self.updateMetaFileFlag(m))
        self.widgets[mode + 'AlignmentResultsFile'].textChanged.connect(lambda d, m=mode: self.updateAlignmentResultsFlag(m))

        execfilename = os.path.join(self.subtomodir, 'ParticlePolishing/reconstructSubtomograms.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='particlePolishSingle', folder=self.logfolder,
                                                      id='SingleParticlePolish')
        paramsCmd = [self.subtomodir, mode + 'numberMpiCores', self.pytompath, mode + 'particleList',
                     mode + 'AlignedTiltDir', mode + 'Template', mode + 'destination', mode + 'BinFactorReconstruction',
                     mode + 'maxParticleShift', mode + 'OffsetX', mode + 'OffsetY', mode + 'OffsetZ', mode + 'FSCFlag',
                     mode + 'MetaFileFlag', mode + 'AlignmentResultsFileFlag',polishParticles]

        [func(mode) for func in (self.updateMeta, self.updateFSCFlag, self.updateMetaFileFlag, self.updateAlignmentResultsFlag)]

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename, paramsSbatch=paramsSbatch)
        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def populatePartPolishingBatchTable(self, key=''):
        self.selectFSCFilesWindow.close()
        particleFilesStart = sorted(self.polishLists.text().split('\n'))
        templates = sorted(self.templateList.text().split('\n'))
        fscfiles = [''] + sorted(self.fscList.text().split('\n'))
        particleFiles = []

        for particleFile in particleFilesStart:
            if '_tomogram_' in particleFile:
                particleFiles.append(particleFile)
            else:
                pLs = extractParticleListsByTomoNameFromXML(particleFile, directory=os.path.dirname(particleFile))
                for pl in pLs:
                    particleFiles.append(pl)

        particleFiles = sorted(particleFiles)

        headers = ["Filename particleList", "Run", "Origin", "Tilt Images",  'Template', 'dimZ', 'Bin factor recon',
                   "Max Shift", "Offset X", "Offset Y", "Offset Z", 'Meta File', 'Align Results', 'FSC Filter', '']
        types = ['txt', 'checkbox', 'combobox', 'combobox', 'combobox', 'lineedit', 'lineedit', 'lineedit', 'lineedit',
                 'lineedit', 'lineedit', 'combobox', 'comboboxF', 'combobox', 'txt']
        a = 40
        sizes = [0, 0, 80, 80, a, a, a, a, a, a, a, a, a, a]

        tooltip = ['Names of the particleList files',
                   'Check this box to run particle polishing.',
                   'From which step do you want to use tilt-images?',
                   'Folder with aligned results.',
                   'Template File',
                   'Z-dimension of tomogram from which particles are picked.',
                   'Binning factor used for the reconstruction.',
                   'Max particle shift (pixel). Set the area used to locate maximum.'
                   'Offset in X-dimension', 'Offset in Y-dimension', 'Offset in Z-dimension',
                   'Meta file from which angles are extracted',
                   'FSC File used for filtering the cross-correlation']

        values = []
        refmarkindices = []

        for n, particleFile in enumerate(particleFiles):
            if not particleFile: continue

            path = os.path.basename(particleFile)
            for i in ('particleList_', 'coords_','_flipped'):
                path = path.replace(i,'')

            base, ext = os.path.splitext(path)

            if '_tomogram_' in base: base = 'tomogram_' + base.split('_tomogram_')[1]

            for t in ('WBP', 'INFR'):
                if t in base: base = base.split(t)[0] + t

            if base + '.mrc' in os.listdir(self.tomogramfolder) or base + '.em' in os.listdir(self.tomogramfolder):
                if os.path.exists(os.path.join(self.tomogramfolder, base + '.mrc')):
                    folder = os.popen('ls -alrt {}.mrc'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                elif os.path.exists(os.path.join(self.tomogramfolder, base + '.em')):
                    folder = os.popen('ls -alrt {}.em'.format(os.path.join(self.tomogramfolder, base))).read()[:-1]
                else:
                    folder = ''
                if not folder: continue


                dimZ = read_size(folder.split()[-1], 'z')
                folder = os.path.dirname(folder.split()[-1])

                # markerfile  = os.path.join(folder, 'markerfile.txt')
                # markerdata = readMarkerfile(markerfile,61)
                sor = os.path.join(os.path.dirname(os.path.dirname(folder)), 'sorted')
                al = os.path.join(os.path.dirname(os.path.dirname(folder)), 'alignment')
                ctf = os.path.join(os.path.dirname(os.path.dirname(folder)), 'ctf')
                choices = [al + '/' + f for f in os.listdir(al) if 'marker_' in f and os.path.isdir(al + '/' + f)]
                # choices = list(map(str,range(markerdata.sizeZ()))) # + ['closest']
                # a = sorted(glob.glob('{}/Reconstruction*-*.out'.format(folder)))[-1]

                try:
                    from lxml import etree
                    xmlObj = etree.parse(particleFile)
                    particles = xmlObj.xpath('Particle')
                    binning = int(particles[0].xpath('InfoTomogram')[0].get('BinningFactor'))
                    refmarkindex = int(particles[0].xpath('InfoTomogram')[0].get('RefMarkIndex'))
                except:
                    print('Default values used for {}:\n\tbin recon = 8\n\t ref mark index = 1'.format(
                        os.path.basename(particleFile)))
                    binning = 8
                    refmarkindex = 1
                # binning = os.popen('cat {} | grep "--referenceMarkerIndex" '.format(a)).read()[:-1]
                # print(binning)
                origin = ['alignment', 'ctf', 'sorted']
                metafiles = [os.path.join(sor, f) for f in os.listdir(sor) if f.endswith('.meta')] + ['']

                import glob
                query = os.path.join(al,'*/alignmentResults.txt')
                alignfiles= [os.path.join(sor, f) for f in glob.glob(query)] + ['']

                values.append([particleFile, True, origin, choices, templates, dimZ, binning, 25, 0, 0, 0, metafiles, alignfiles, fscfiles, ''])
                refmarkindices.append(refmarkindex)

        try:
            self.num_nodes[key].setParent(None)
        except:
            pass

        if values:
            self.valuesBatchParticlePolishing = values
            self.fill_tab(key, headers, types, values, sizes, tooltip=tooltip, nn=True)
            self.tabPolish_widgets = self.tables[key].widgets
            for n, index in enumerate(refmarkindices):
                self.tabPolish_widgets['widget_{}_3'.format(n)].setCurrentIndex(index)

            for row in range(len(values)):
                w = self.tabPolish_widgets['widget_{}_2'.format(row)]
                w.currentIndexChanged.connect(lambda d, r=row, v=values: self.updateChoicesPP(r, v))

            self.particleFilesBatchPolish = particleFiles
            self.pbs[key].clicked.connect(lambda dummy, pid=key, v=values: self.massPolishParticles(pid, v))
        else:
            return

        for i in range(len(values)):
            self.updateChoicesPP(i, values)

    def addFRMAlignmentFields(self, mode=None):
        title = "FRM Alignment"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 40, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB)
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList', initdir=self.pickpartdir,
                                    tooltip='Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask', initdir=self.frmdir,
                                    tooltip='Select the mask file.', mode='file', filetype=['em', 'mrc'], cstep=1,
                                    rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask'])
        self.insert_label_line_push(parent, 'Filename Average', mode + 'filenameAverage', initdir=self.frmdir,
                                    tooltip='Choose a filename for the average of all particles.', mode='file',
                                    filetype=['em', 'mrc'], cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Average', rstep=1, cstep=-3, action=self.gen_average,
                               params=[mode + 'particleList', mode + 'filenameAverage', mode + 'outputDir'])

        self.insert_label_line_push(parent, 'Output Directory', mode + 'outputDir',
                                    'Folder in which the output will be written.')
        self.insert_label(parent, rstep=1, cstep=0)
        self.insert_label_spinbox(parent, mode + 'bwMin', 'Min Order SH Polynomial',
                                  value=8, minimum=0, stepsize=1,
                                  tooltip='The minimal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'bwMax', 'Max Order SH Polynomial',
                                  value=64, minimum=0, stepsize=1,
                                  tooltip='The maximal order of the polynomial used for spherical harmonics alignment.')
        self.insert_label_spinbox(parent, mode + 'frequency', 'Frequency (px)',
                                  value=8, minimum=0, stepsize=1,
                                  tooltip='The minimal frequency used for reconstruction.')
        self.insert_label_spinbox(parent, mode + 'maxIterations', 'Maximum Iterations',
                                  value=8, minimum=1, stepsize=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
        self.insert_label_spinbox(parent, mode + 'peakOffset', 'Peak Offset',
                                  value=0, minimum=0, stepsize=1,
                                  tooltip='Sets the peak offset.')
        self.insert_label(parent, rstep=1, cstep=0)
        self.insert_label_spinbox(parent, mode + 'pixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox, minimum=0.1, stepsize=0.1, value=1.75)
        self.insert_label_spinbox(parent, mode + 'particleDiameter', 'Particle Diameter (A)', rstep=1, cstep=0,
                                  minimum=10, stepsize=1, value=300, maximum=10000, width=150)

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updateFRM(m))

        rscore = 'False'
        weightedAv = 'False'
        weighting = ''
        binning_mask = '1'
        sphere = 'True'
        ad_res = '0.00'
        fsc = '0.50'

        jobfilename = [mode + 'outputDir', 'job_description.xml']  # os.path.join(self.frmdir, 'job_description.xml')
        exefilename = [mode + 'outputDir', 'frmAlignment.sh']  # os.path.join(self.frmdir, 'frmAlignment.sh')

        paramsSbatch = guiFunctions.createGenericDict(fname='FRMAlign', folder=self.logfolder,
                                                      id='FRMAlignment')  # , modules=['openmpi/2.1.1', 'python/2.7', 'lib64/append', 'pytom/dev/gui'])
        paramsJob = [mode + 'bwMin', mode + 'bwMax', mode + 'frequency', mode + 'maxIterations', mode + 'peakOffset',
                     rscore, weightedAv, mode + 'filenameAverage', weighting, mode + 'filenameMask', binning_mask,
                     sphere,
                     mode + 'pixelSize', mode + 'particleDiameter', mode + 'particleList', mode + 'outputDir']
        paramsCmd = [self.subtomodir, self.pytompath, jobfilename, mode + 'numberMpiCores', templateFRMSlurm]

        self.insert_gen_text_exe(parent, self.stage, xmlfilename=jobfilename, jobfield=True, exefilename=exefilename,
                                 paramsXML=paramsJob + [templateFRMJob], paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_inputFiles', groupbox)
        return groupbox

    def addGLocalAlignmentFields(self,mode):
        title = "GLocal alignment"
        tooltip = 'Run pytom GLocal routine.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Initial reference model', mode + 'referenceModel', mode='file',
                                    filetype=['em','mrc'], enabled=True,
                                    tooltip='Reference : the initial reference - if none provided average of particle list')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask'])
        self.insert_label_line_push(parent, 'Output Directory', mode + 'destination', mode='folder',
                                    tooltip='Select the destination directory.')
        self.insert_label_spinbox(parent, mode + 'numIterations', 'Number of Iterations',
                                  minimum=1, value=4, stepsize=1,
                                  tooltip='Sets the number of iterations.')
        self.insert_label_spinbox(parent, mode + 'pixelSize', 'Pixel Size (A)',
                                  wtype=QDoubleSpinBox, minimum=0.1, value=1.75, stepsize=0.1,
                                  tooltip='Pixelsize in Angstrom ')
        self.insert_label_spinbox(parent, mode + 'particleDiameter', 'Particle Diameter (A)',
                                  minimum=10, stepsize=1, value=300, maximum=10000,
                                  rstep=1, cstep=-1, tooltip='Particle diameter in Angstrom.')
        self.insert_label_spinbox(parent, mode + 'angleShells', 'Number of angular shells',
                                  tooltip='# Angle shells used for angular refinement.',
                                  minimum=1, stepsize=1, value=3, maximum=100,
                                  rstep=1, cstep=-1)
        self.insert_label_spinbox(parent, mode + 'angleIncrement', 'Angular Increment (degrees)',
                                  minimum=1, stepsize=1, value=3, maximum=359, wtype=QDoubleSpinBox,
                                  rstep=1, cstep=-1, tooltip='Angular increment for refinement.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=-1,
                                  stepsize=1,minimum=1,value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0,
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.widgets[mode+'jobName'] = QLineEdit()
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets['referenceCommand'] = QLineEdit(self)
        self.widgets['referenceCommand'].setVisible(False)

        self.widgets[mode + 'referenceModel'].textChanged.connect(lambda dummy, m=mode: self.updateReference(m))
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.updatePixelSize(m))
        self.widgets[mode + 'destination'].textChanged.connect(lambda d, m=mode: self.updateJobname(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        self.updateJobname(mode)

        exefilename = [mode+'destination', 'GLocal_Alignment.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='GLocal', folder=self.logfolder, id='GLocalAlignment')
        paramsCmd = [self.subtomodir, self.pytompath, self.pytompath, mode+'particleList', 'referenceCommand',
                     mode+'filenameMask', mode+'numIterations', mode+'pixelSize', mode+'particleDiameter',
                     mode+'binning', mode+'jobName', mode+'destination', mode + 'angleShells',
                     mode + 'angleIncrement', mode + 'numberMpiCores', mode + 'gpuString',  templateGLocal]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        self.updateGpuString(mode)
        self.updateReference(mode)
        self.updatePixelSize(mode)

        setattr(self, mode + 'gb_GLocal', groupbox)
        return groupbox

    def addCPCAFields(self,mode=''):
        title = "Classify CPCA"
        tooltip = 'The CCC is further used for classification. This script computes \n'+\
                  'the eigenvectors of the CCC and projects the data on the first \n'+\
                  'neig eigenvectors. Subsequently, these multidimensional vectors \n'+\
                  'are clustered into nclass groups using a kmeans algorithm.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_line(parent, 'Filename Output particleList', mode + 'outputFilename',
                               tooltip='Filename for generated XML file that includes the assigned classes for each particle. No full path needed.')
        self.insert_label_line_push(parent, 'CCC File', mode + 'cccFile',
                                    'Select the constrained correlation matrix file from the previous step.', mode='file', filetype='csv')
        self.insert_label_spinbox(parent, mode + 'numEig', text='Number of Eigenvectors',
                                  value=4, minimum=1, stepsize=1, maximum=100,
                                  tooltip='Sets the number of eigenvectors (corresponding to largest eigenvectors) used for clustering.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes',
                                  value=4, minimum=1,stepsize=1, maximum=1000,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_line(parent, 'Prefix', mode + 'prefix', rstep=1, cstep=0,
                               tooltip='Root for generated averages of the corresponding classes. The files will be called "Prefix"_iclass.em.')

        # self.widgets[mode + 'particleList'].textChanged.connect(
        #     lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode, p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        exefilename = [mode + 'outFolder', 'CPCA_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CPCA', folder=self.logfolder, id='CPCA')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'outputFilename',
                     mode + 'cccFile', mode + 'numEig', mode+'numClasses', mode + 'outFolder', mode+'prefix',  templateCPCA]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        setattr(self, mode + 'gb_CPCA', groupbox)
        return groupbox

    def addCCCFields(self,mode=''):
        title = "Pairwise Constrained Cross Correlation"
        tooltip = 'Calculate the pairwise constrained cross correlation.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Filename Mask', mode + 'filenameMask',mode='file', filetype=['em', 'mrc'],
                                    tooltip='Select the mask file.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameMask'])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'lowpass', 'Lowpass filter (px)',
                                  minimum=0, maximum=1024, stepsize=1, value=20,
                                  tooltip='The lowpass filter is applied to all subtomograms after binning.')
        self.insert_label_spinbox(parent, mode + 'binning', 'Binning Factor', rstep=1, cstep=-1,
                                  minimum=1, stepsize=1, value=1,
                                  tooltip='Perform binning (downscale) of subvolumes by factor. Default=1.')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0,
                               tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')
        self.widgets[mode + 'gpuString'] = QLineEdit('')

        # self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.cpcadir: self.updateOutFolder(mode,p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        exefilename = [mode + 'outFolder', 'CCC_Classification.sh']
        paramsSbatch = guiFunctions.createGenericDict(fname='CCC_Class', folder=self.logfolder,
                                                      id='PairwiseCrossCorrelation')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'filenameMask',
                     mode + 'lowpass', mode + 'binning', mode + 'outFolder', mode + 'numberMpiCores', mode + 'gpuString',
                     templateCCC]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)
        self.updateGpuString(mode)

        setattr(self, mode + 'gb_CCC', groupbox)
        return groupbox

    def addAC3DFields(self,mode):
        title = "Autofocussed Classification"
        tooltip = 'Run autofocussed classification.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list.', mode='file', filetype='xml')
        self.insert_label_line_push(parent, 'Classification Mask', mode + 'filenameClassificationMask', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='This mask is used for constraining the calculation of the focused mask. (Optional)', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameClassificationMask'])
        self.insert_label_line_push(parent, 'Alignment Mask', mode + 'filenameAlignmentMask', mode='file', filetype=['em', 'mrc'], enabled=True,
                                    tooltip='This mask is only used for the alignment purpose. Only specify it if the particle list is not aligned.', cstep=1, rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_mask,
                               params=[mode + 'filenameAlignmentMask'])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'numClasses', 'Number of Classes', stepsize=1, value=4, minimum=2,
                                  tooltip='Number of classes used for kmeans classification.')
        self.insert_label_spinbox(parent, mode + 'maxIterations', 'Number Of Iterations',stepsize=1,value=10, minimum=1,
                                  tooltip='Sets the maximal number of iterations of alignmment.')
        self.insert_label_spinbox(parent, mode + 'binningFactor', 'Binning Factor',stepsize=1,value=1, minimum=1,
                                  tooltip='Sets the binning factor of the input data.')
        self.insert_label_spinbox(parent, mode + 'bwMax', 'Max Bandwidth Reconstruction (px)',
                                  minimum=1, maximum=1024, stepsize=1, value=20,
                                  tooltip='The maximal frequency used for reconstruction.')
        self.insert_label_spinbox(parent, mode + 'peakOffset', 'Max Peak Offset', stepsize=1,value=10, minimum=1,
                                  tooltip='Sets the peak offset.')
        self.insert_label_spinbox(parent, mode+'noisePercentage', 'Noise Percentage',
                                  wtype=QDoubleSpinBox, value=0.1, stepsize=.1,
                                  tooltip='Noise percentage (between 0 and 1). If you estimate your dataset contains certain amount of noise outliers, specify it here.')
        self.insert_label_spinbox(parent, mode + 'partDensThresh', 'Particle Density Threshold',
                                  wtype=QDoubleSpinBox, value=0., minimum=-6.0, maximum=6.0, stepsize=.1,
                                  tooltip='Particle density threshold for calculating the difference map (optional, by default 0).\n'
                                          'Two other most common choise are -2 and 2. -2 means all the values of the subtomogram \n'
                                          'below the -2 sigma will be used for calculating the difference mask (negative values\n'
                                          'count). 2 means all the values of the subtomogram above the 2 sigma will be used for\n'
                                          'calculating t he difference mask (positive values count). Finally, 0 means all the \n'
                                          'values are used for the calculation.')
        self.insert_label_spinbox(parent, mode + 'stdDiffMap', 'STD Threshold Diff Map', rstep=1, cstep=0,
                                  wtype=QDoubleSpinBox, stepsize=.1, minimum=0, maximum=1, value=0.4,
                                  tooltip='STD threshold for the difference map (optional, by default 0.4). This value should be \n'
                                          'between 0 and 1. 1 means only the place with the peak value will be set to 1 in the \n'
                                          'difference map (too much discriminative ability). 0 means all the places with the value\n'
                                          'above the average of STD will be set to 1 (not enough discriminative ability).')

        # Connected Widgets
        self.widgets[mode + 'filenameAlignmentMask'].textChanged.connect(
            lambda d, m=mode: self.updateAlignmentMaskFlag(m))
        self.widgets[mode + 'filenameClassificationMask'].textChanged.connect(
            lambda d, m=mode: self.updateClassificationMaskFlag(m))
        # self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode, p=self.acpath: self.updateOutFolder(mode, p))
        self.widgets[mode + 'outFolder'].textChanged.connect(lambda d, m=mode: self.createOutFolder(m))

        # Widgets Updated When Other Widgets Are Updated
        self.widgets[mode + 'flagAlignmentMask'] = QLineEdit('')
        self.widgets[mode + 'flagClassificationMask'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')

        # Parameters for execution
        exefilename = [mode + 'outFolder', 'AC_Classification.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='AutoFocus', folder=self.logfolder,
                                                      id='AutoFocusClassification')
        paramsCmd = [self.subtomodir, self.pytompath, mode + 'particleList', mode + 'flagAlignmentMask',
                     mode + 'flagClassificationMask', mode + 'numClasses', mode + 'bwMax', mode + 'maxIterations',
                     mode + 'peakOffset', mode + 'noisePercentage', mode + 'partDensThresh', mode + 'stdDiffMap',
                     mode + 'outFolder', mode + 'numberMpiCores', mode + 'binningFactor', templateAC]


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        # Run Update With Data From Logfile
        self.updateAlignmentMaskFlag(mode)
        self.updateClassificationMaskFlag(mode)
        self.updateOutFolder(mode, self.acpath)

        setattr(self, mode + 'gb_AC', groupbox)
        return groupbox

    def addFSCFields(self, mode):
        title = "Fourier Shell Correlation"
        tooltip = 'Run fourier shell correlation function.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        # Insert Parameter Widgets

        self.insert_label_line_push(parent, 'Volume 1', mode + 'volume1', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='The first volume path. (Optional)')
        self.insert_label_line_push(parent, 'Volume 2', mode + 'volume2', mode='file',
                                    filetype=['em', 'mrc'], enabled=True,
                                    tooltip='The second volume path. (Optional)')
        self.insert_label_line_push(parent, 'Particle List', mode + 'particleList',
                                    'Select the particle list if v1 and v2 are not available.', mode='file',
                                    filetype='xml')
        self.insert_label_line_push(parent, 'Mask', mode + 'mask',
                                    'Mask (optional, but recomended).', mode='file',
                                    filetype=['em', 'mrc'],cstep=1,rstep=0)
        self.insert_pushbutton(parent, 'Create', rstep=1, cstep=-3, action=self.gen_fsc_mask,
                               params=[mode + 'mask'])
        self.insert_label_line_push(parent, 'Output Folder', mode + 'outFolder', mode='folder',
                                    tooltip='Select/Create an output folder.')
        self.insert_label_spinbox(parent, mode + 'fsc', 'FSC cutoff', stepsize=1, value=0.143, minimum=0, maximum=1.,
                                  wtype=QDoubleSpinBox,decimals=3,
                                  tooltip='The FSC criterion. Value between 0.0 and 1.0. Standard values are 0.5, 0.3. or 0.17')
        self.insert_label_spinbox(parent, mode + 'pixelsize', 'Pixelsize (A)',stepsize=1,value=2.62, minimum=.1, decimals=3,
                                  tooltip='Pixelsize in Angstrom', wtype=QDoubleSpinBox)
        self.insert_label_spinbox(parent, mode + 'randomizePhases', 'randomizePhases', stepsize=1,value=0.8, minimum=0, maximum=1,
                                  tooltip='Check validity of FSC using phases randomization beyond spatial frequency'
                                          ' where the uncorrected FSC curve drops below set threshold. Values between '
                                          '0, 1. No phase randomization check performed when value is set to 0.', decimals=3,
                                  wtype=QDoubleSpinBox)
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', cstep=-1,
                                  tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.insert_label_checkbox(parent, mode + 'CombinedResolution', 'Combined Resolution',
                                   tooltip='Check this box to calculate the resolution of the two volumes together.')

        self.insert_label_checkbox(parent, mode + 'plot', 'Plot Results',
                                   tooltip='Check this box to plot the results.', cstep=0, rstep=1)


        # Connected Widgets
        self.widgets[mode + 'volume1'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 1))
        self.widgets[mode + 'volume2'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 2))
        self.widgets[mode + 'particleList'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 0))
        self.widgets[mode + 'mask'].textChanged.connect(
            lambda d, m=mode: self.updateFSCFlags(m, 3))
        self.widgets[mode + 'plot'].stateChanged.connect(
            lambda d, m=mode: self.updateFSCPlotFlag(m))
        self.widgets[mode + 'CombinedResolution'].stateChanged.connect(
            lambda d, m=mode: self.updateFSCCombResFlag(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        # Widgets Updated When Other Widgets Are Updated
        self.widgets[mode + 'flagVolume1'] = QLineEdit('')
        self.widgets[mode + 'flagVolume2'] = QLineEdit('')
        self.widgets[mode + 'flagParticleList'] = QLineEdit('')
        self.widgets[mode + 'flagMask'] = QLineEdit('')
        self.widgets[mode + 'flagPlot'] = QLineEdit('')
        self.widgets[mode + 'flagCombinedResolution'] = QLineEdit('')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets[mode + 'numberMpiCores'] = QLineEdit('20')

        # Parameters for execution

        exefilename = [mode + 'outFolder', 'FSC_Validation.sh'] #os.path.join(acpath, 'AC_Classification.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='FSC', folder=self.logfolder,
                                                      id='FSCValidation')
        paramsCmd = [self.fscdir, self.pytompath, mode + 'flagParticleList', mode + 'flagVolume1',
                     mode + 'flagVolume2', mode + 'flagMask', mode + 'outFolder', mode + 'fsc',
                     mode + 'pixelsize', mode + 'randomizePhases', mode + 'flagPlot', mode + 'flagCombinedResolution',
                     mode + 'gpuString', templateFSC]


        # Generation of textboxes and pushbuttons related to submission
        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=exefilename, paramsCmd=paramsCmd,
                                 paramsSbatch=paramsSbatch)

        self.widgets[mode+'queue'].stateChanged.connect(lambda d, m=mode: self.updateLog(m))

        # Run Update With Data From Logfile
        self.updateFSCFlags(mode)
        self.updateFSCPlotFlag(mode)
        self.updateLog(mode)
        self.updateFSCCombResFlag(mode)
        self.updateGpuString(mode)

        setattr(self, mode + 'gb_FSC', groupbox)
        return groupbox

    # Helper functions

    def updatePolishFlag(self, mode):
        fname = self.widgets[mode + 'polishFile'].text()
        if not fname:
            self.widgets[mode + 'polishFlag'].setText('')
            return

        self.widgets[mode + 'polishFlag'].setText(f'--particlePolishFile {fname}')

    def updateFSCFlag(self, mode):

        fname = self.widgets[mode + 'FSCPath'].text()
        if not fname:
            self.widgets[mode + 'FSCFlag'].setText('')
            return

        self.widgets[mode + 'FSCFlag'].setText(f'--FSCPath {fname}')

    def updateMetaFileFlag(self, mode):

        fname = self.widgets[mode + 'MetaFile'].text()
        if not fname:
            self.widgets[mode + 'MetaFileFlag'].setText('')
            return

        self.widgets[mode + 'MetaFileFlag'].setText(f'--metaFile {fname}')

    def updateAlignmentResultsFlag(self, mode):

        fname = self.widgets[mode + 'AlignmentResultsFile'].text()
        if not fname:
            self.widgets[mode + 'AlignmentResultsFileFlag'].setText('')
            return

        self.widgets[mode + 'AlignmentResultsFileFlag'].setText(f'--alignmentResultsFile {fname}')

    def updateMeta(self,mode):
        pl = self.widgets[mode + 'particlelist'].text()
        try: 
            tomoID = int(pl.split('_tomogram_')[-1][:3])
            tomo = os.path.join(self.tomogram_folder, 'tomogram_{:03d}/sorted/'.format(tomoID))
            a = glob.glob(tomo+'*.meta')
            if not a: print('No meta file found. auto update stopped.')
            a = a[-1]
            self.widgets[mode+'MetaFile'].setText(a)
        except:
            pass

    def updateHeaderChoices(self, rowID, key):
        print('Update Header', rowID, key)
        header, row = self.tables[key].general_widgets[3], self.tab12_widgets['widget_{}_3'.format(rowID)]
        AllItemsGeneral = [header.itemText(i) for i in range(header.count())]
        indPart = [row.itemText(i) for i in range(row.count())]
        header.clear()

        header.addItems([item for item in indPart if not item in AllItemsGeneral and 'CLOSEST' in item])
        header.addItems(AllItemsGeneral)

    def updateAlignmentType(self, rowID, table_id):
        print('Update alignment type')
        values = self.valuesBatchSubtomoReconstruction
        self.tab12_widgets['widget_{}_4'.format(rowID)].clear()
        current_index_f = self.tab12_widgets['widget_{}_2'.format(rowID)].currentIndex()
        current_index_a = self.tab12_widgets['widget_{}_3'.format(rowID)].currentIndex()
        folder = values[rowID][2][current_index_f]
        aligntype = values[rowID][3][current_index_a]

        if not '_closest_' in os.path.basename(folder.lower()):
            if aligntype:
                temp = os.path.join(folder, aligntype)
                origin = [f'{temp}/{f}' for f in os.listdir(temp) if os.path.isdir(f'{temp}/{f}')]
            else:
                origin = []
        else:
            origin = []
            angles = os.path.basename(folder).split('_')[-1]
            if aligntype:
                for ff in values[rowID][2]:
                    temp_angles = os.path.basename(ff).split('_')[-1]
                    if 'CLOSEST' in ff or temp_angles != angles:
                        continue
                    tempaligntype = [f for f in os.listdir(ff) if os.path.isdir(ff + '/' + f)]

                    if tempaligntype:
                        temporigin = [f for f in os.listdir(f'{ff}/{tempaligntype[0]}') if os.path.isdir(f'{ff}/{tempaligntype[0]}/{f}')]
                    else:
                        temporigin = []

                    for to in temporigin:
                        if not to in origin:
                            origin.append(to)

        self.valuesBatchSubtomoReconstruction[rowID][4] = origin
        for item in origin:
            self.tab12_widgets['widget_{}_4'.format(rowID)].addItem(os.path.basename(item))

    def updateChoices(self, rowID, table_id):
        print(f'Update Choices {rowID}')

        values = self.valuesBatchSubtomoReconstruction
        self.tab12_widgets['widget_{}_3'.format(rowID)].clear()
        self.tab12_widgets['widget_{}_4'.format(rowID)].clear()

        current_index = self.tab12_widgets['widget_{}_2'.format(rowID)].currentIndex()

        try:
           folder = values[rowID][2][current_index]
        except:
            return

        particleFile = values[rowID][0]
        a = 'tomogram_' + particleFile.split('_tomogram_')[1][:3]

        print('update for folder', folder)

        if not '_closest_' in os.path.basename(folder.lower()):
            aligntype = [folder + '/' + f for f in os.listdir(folder) if os.path.isdir(folder + '/' + f)]
            if aligntype:
                origin = [f'{aligntype[0]}/{f}' for f in os.listdir(aligntype[0]) if os.path.isdir(f'{aligntype[0]}/{f}')]
            else:
                origin = []
        else:
            aligntype = []
            origin = []
            angles = os.path.basename(folder).split('_')[-1]
            for ff in values[rowID][2]:
                temp_angles = os.path.basename(ff).split('_')[-1]
                if 'CLOSEST' in ff or temp_angles != angles:
                    continue
                tempaligntype = [f for f in os.listdir(ff) if os.path.isdir(ff + '/' + f)]

                if tempaligntype:
                    temporigin = [f for f in os.listdir(f'{ff}/{tempaligntype[0]}') if os.path.isdir(f'{ff}/{tempaligntype[0]}/{f}')]
                else:
                    temporigin = []

                for ta in tempaligntype:
                    if not ta in aligntype:
                        aligntype.append(ta)
                for to in temporigin:
                    if not to in origin:
                        origin.append(to)

        #
        # closest_choices = {}
        #
        # for choice in choices:
        #     f = os.path.dirname(choice)
        #     try: a = choice.split('marker_')[1]
        #     except: continue
        #     if ',' in a:
        #
        #         try:
        #             aS, aE = a.split('_')[1].split(',')
        #             angleS = '{:4.1f}'.format( float(aS) )
        #             angleE = '{:4.1f}'.format( float(aE) )
        #
        #             if 'ctf' in a: ctf = '_ctf'
        #             else: ctf = ''
        #
        #
        #             key = '{}/marker_CLOSEST_{},{}{}'.format(f, angleS, angleE, ctf)
        #             value = '{}/marker_CLOSEST_{},{}{}'.format(f, aS, aE, ctf)
        #
        #
        #             closest_choices[key] = value
        #         except:
        #             pass
        #
        #     elif 'marker_' in a:
        #         closest_choices['{}/marker_CLOSEST'.format(f)] = 1
        #     elif 'sorted' in a:
        #         closest_choices[choice] = 1
        #
        # choices  = list(closest_choices.values()) + choices
        #
        self.valuesBatchSubtomoReconstruction[rowID][3] = aligntype
        self.valuesBatchSubtomoReconstruction[rowID][4] = origin

        for item in aligntype:
            self.tab12_widgets['widget_{}_3'.format(rowID)].addItem(os.path.basename(item))
        # for item in origin:
        #     self.tab12_widgets['widget_{}_4'.format(rowID)].addItem(os.path.basename(item))
        if aligntype:
            self.updateAlignmentType(rowID, table_id)
        self.updateHeaderChoices(rowID, table_id)

    def updateChoicesPP(self, rowID, values):
        values = self.valuesBatchParticlePolishing
        self.tabPolish_widgets['widget_{}_3'.format(rowID)].clear()

        current_index = self.tabPolish_widgets['widget_{}_2'.format(rowID)].currentIndex()

        origin = values[rowID][2][current_index]
        particleFile = values[rowID][0]
        a = 'tomogram_' + particleFile.split('_tomogram_')[1][:3]

        folder = os.path.join(self.tomogram_folder, a, origin)

        choices = [folder + '/' + f for f in os.listdir(folder) if
                    'marker_' in f and os.path.isdir(folder + '/' + f)]

        self.valuesBatchParticlePolishing[rowID][3] = choices

        for item in choices:
            self.tabPolish_widgets['widget_{}_3'.format(rowID)].addItem(os.path.basename(item))

    def updateFRM(self, mode):
        item = self.widgets[mode + 'particleList'].text()
        if not item:
            return
        folder, ext = os.path.splitext( os.path.basename(item))
        outputDir = os.path.join(self.frmdir, folder.replace('particleList_', ''))
        if not os.path.exists(outputDir): os.mkdir(outputDir)
        else:
            for i in range(1000):
                if not os.path.exists(outputDir+'_{:03d}'.format(i)):
                    outputDir += '_{:03d}'.format(i)
                    os.mkdir(outputDir)
                    break


        self.widgets[mode+'outputDir'].setText(outputDir)

    def updateReference(self, mode):
        if self.widgets[mode + 'referenceModel'].text():
            self.widgets['referenceCommand'].setText( "--reference " + self.widgets[mode + 'referenceModel'].text())
        else:
            self.widgets['referenceCommand'].setText("")

    def updateGpuString(self, mode):
        id = self.widgets[mode + 'gpuID'].text()
        try:
            a = list(map(int,[el for el in id.split(',') if el != '']))
        except:
            self.widgets[mode + 'gpuID'].setText('')
            self.popup_messagebox('Warning', 'Invalid value in field', 'Impossible to parse gpu IDs, field has been cleared.')
            return

        if len(id) > 0:
            self.widgets[mode + 'gpuString'].setText(f' --gpuID {id}')
            self.widgets[mode + 'numberMpiCores'].setText(f'{len(a)+1}')
        else:
            self.widgets[mode + 'gpuString'].setText('')

            if 'gLocal'in mode and 'GLocalAlignment' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['GLocalAlignment'].values()
            elif 'CCC' in mode and 'PairwiseCrossCorrelation' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['PairwiseCrossCorrelation'].values()
            elif 'FSC' in mode and 'FSCValidation' in self.qparams.keys():
                qname, n_nodes, cores, time, modules = self.qparams['FSCValidation'].values()
            else:
                cores = 20
            self.widgets[mode + 'numberMpiCores'].setText(f'{cores}')


    def updatePixelSize(self, mode):
        from pytom.basic.structures import ParticleList
        try:
            particleList = self.widgets[mode + 'particleList'].text()
            pl = ParticleList()
            pl.fromXMLFile(particleList)
            tomid = pl[0].getPickPosition().getOriginFilename().split('/tomogram_')[1].split('_')[0]
            metaquery = os.path.join( self.tomoanalysis, f'tomogram_{tomid}/sorted/*.meta')
            a = glob.glob(metaquery)
            if len(a):
                pixelsize = guiFunctions.loadstar(a[0], dtype=guiFunctions.datatype)['PixelSpacing'][0]
                self.widgets[mode + 'pixelSize'].setValue(pixelsize)
        except:
            pass

    def updateJobname(self, mode):
        dest = self.widgets[mode+'destination'].text()
        if dest:
            self.widgets[mode + 'jobName'].setText( os.path.join(dest, 'glocal_input_params_'+os.path.basename(dest)+'.xml'))

    def updateClassificationMaskFlag(self, mode):
        try:
            mask = self.widgets[mode+'filenameClassificationMask'].text()
            if mask:
                self.widgets[mode + 'flagClassificationMask'].setText('\\\n-c {}'.format(mask) )
            else:
                self.widgets[mode + 'flagClassificationMask'].setText('')
        except:
            pass

    def updateAlignmentMaskFlag(self, mode):
        try:
            mask = self.widgets[mode+'filenameAlignmentMask'].text()
            if mask:
                self.widgets[mode + 'flagAlignmentMask'].setText('-m {}'.format(mask) )
            else:
                self.widgets[mode + 'flagAlignmentMask'].setText('-a')
        except:
            pass

    def updateOutFolder(self,mode, path, name='particleList'):
        pl = self.widgets[mode + name].text()
        if not pl or not os.path.exists(pl):
            self.widgets[mode + name].setText('')
            return

        if not self.doAllParticlesExist(pl):
            self.popup_messagebox('Error', 'Particles do not exist',
                                  f'Not all particles in {os.path.basename(pl)} have been created')
            self.widgets[mode + name].setText('')
            return

        folder = os.path.basename(pl).replace('particleList_','')[:-4]
        folder = os.path.join(path,folder)
        if os.path.exists(folder):
            folder = folder
            if not folder: return

        if not os.path.exists(folder):
            os.mkdir(folder)
            #os.system('ln -s {}/Subtomograms {}/Subtomograms'.format(self.subtomodir, folder ) )
        self.widgets[mode + 'outFolder'].setText(folder) 

    def createOutFolder(self, mode):
        folder = self.widgets[mode + 'outFolder'].text()
        if not folder: return
        if not os.path.exists(folder):
            os.mkdir(folder)

    def doAllParticlesExist(self, particleList):
        from pytom.basic.structures import ParticleList
        pl = ParticleList()
        pl.fromXMLFile(particleList)
        allParticlesExist = True
        for particle in pl:
            if not os.path.exists(os.path.join(self.subtomodir ,particle.getFilename())):
                allParticlesExist = False
                break
        return allParticlesExist

    def massExtractParticles(self, pid, values):
        num_nodes = int(self.num_nodes[pid].value())
        qIDs = []
        num_submitted_jobs = nsj = 0
        values = self.valuesBatchSubtomoReconstruction
        try:
            for row in range(self.tables[pid].table.rowCount()):
                if self.tab12_widgets['widget_{}_1'.format(row)].isChecked():
                    particleXML = values[row][0] #[row][0]
                    tomoindex = particleXML.split('_tomogram_')[-1][:3]
                    markerPath = values[row][2][self.tab12_widgets['widget_{}_{}'.format(row,2)].currentIndex()]
                    folder_align_type = values[row][3][self.tab12_widgets['widget_{}_{}'.format(row, 3)].currentIndex()]
                    folder_origin      = values[row][4][self.tab12_widgets['widget_{}_{}'.format(row, 4)].currentIndex()]

                    metafile = glob.glob('{}/03_Tomographic_Reconstruction/tomogram_{}/sorted/*.meta'.format(self.projectname,tomoindex))
                    if not metafile: continue
                    metafile = metafile[0]
                    #q = '{}/03_Tomographic_Reconstruction/tomogram_{}/alignment/unweighted_unbinned_marker_{}'
                    #folder_aligned = q.format(self.projectname,tomoindex,ref_marker)
                    bin_read = self.tab12_widgets['widget_{}_{}'.format(row,5)].text()
                    weight = self.tab12_widgets['widget_{}_{}'.format(row, 6)].text()
                    size = self.tab12_widgets['widget_{}_{}'.format(row, 7)].text()
                    bin_subtomo = self.tab12_widgets['widget_{}_{}'.format(row,8)].text()
                    offx = self.tab12_widgets['widget_{}_{}'.format(row, 9)].text()
                    offy = self.tab12_widgets['widget_{}_{}'.format(row, 10)].text()
                    offz = self.tab12_widgets['widget_{}_{}'.format(row, 11)].text()
                    ppfile = values[row][12][self.tab12_widgets['widget_{}_{}'.format(row, 12)].currentIndex()]
                    ppflag = f'--particlePolishResultFile {ppfile}' if ppfile else ''
                    refid = os.path.basename(markerPath).split('_')[1]

                    if refid.lower() != 'closest':

                        outname = 'Reconstruction/reconstruct_subtomograms_{:03d}_refmarker_{}.sh'.format(int(tomoindex),refid)
                        execfilename = os.path.join(self.subtomodir, outname)
                        qname, n_nodes, cores, time, modules = self.qparams['BatchSubtomoReconstruct'].values()

                        paramsCmd = [particleXML, folder_origin, bin_read, size, bin_subtomo, offx, offy, offz,
                                     self.subtomodir, weight, metafile,  str(cores*n_nodes), ppflag]

                        txt = extractParticles.format(d=paramsCmd)
                        jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder,
                                                               name='SubtomoRecon_{}'.format(nsj % num_nodes),
                                                               num_jobs_per_node=cores, time=time, partition=qname,
                                                               modules=modules, num_nodes=n_nodes, singleton=True) + txt

                    else:
                        outname = 'Reconstruction/reconstruct_subtomograms_{:03d}_refmarker_{}.sh'
                        outname = outname.format(int(tomoindex), refid)
                        execfilename = os.path.join(self.subtomodir, outname)

                        folder_origin = os.path.join(markerPath, folder_align_type, folder_origin)


                        tomodir = os.path.dirname(os.path.dirname(metafile))

                        reconAlg = None
                        for alg in ('WBP', 'INFR'):
                            if alg in particleXML:
                                reconAlg = alg

                        if not reconAlg:
                            print('FAIL: subtomogram reconstruction for {} failed. No INFR or WBP in xml path to particle.')
                            continue

                        # Find the most recent markerLocation file
                        end = 'reconstruction/{}/markerLocations*_irefmark_*.txt'
                        end = end.format(reconAlg, tomoindex, reconAlg)
                        logfilequery = os.path.join(tomodir, end)
                        logfiles = glob.glob(logfilequery)
                        logfiles.sort(key=os.path.getmtime)
                        logfile = logfiles[-1]

                        qname, n_nodes, cores, time, modules = self.qparams['BatchSubtomoReconstruct'].values()

                        paramsCmd = [particleXML, folder_origin, bin_read, size, bin_subtomo, offx, offy, offz,
                                     self.subtomodir, weight, metafile, logfile, str(cores*n_nodes), 'sorted_aligned']


                        txt = extractParticlesClosestMarker.format(d=paramsCmd)

                        jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder, singleton=True, num_nodes=n_nodes,
                                                               name='SubtomoRecon_{}'.format(nsj % num_nodes),
                                                               partition=qname,num_jobs_per_node=cores, time=time,
                                                               modules=modules) + txt


                    ID, num = self.submitBatchJob(execfilename, pid, jobtxt)
                    nsj += 1
                    qIDs.append(ID)

            if nsj:
                self.popup_messagebox('Info', 'Submission Status', f'Submitted {nsj} jobs to the queue.')
                self.addProgressBarToStatusBar(qIDs, key='QJobs', job_description='Subtom Recon Batch')

        except Exception as e:
            print('Submission failed due to following error: ')
            print(e)

    def massPolishParticles(self, pid, values):
        qIDs =[]
        num_nodes = int(self.num_nodes[pid].value())
        num_submitted_jobs = nsj = 0
        values = self.valuesBatchParticlePolishing
        for row in range(self.tables[pid].table.rowCount()):
            if self.tabPolish_widgets['widget_{}_1'.format(row)].isChecked():
                particleXML = values[row][0]  # [row][0]
                tomoindex = particleXML.split('_tomogram_')[-1][:3]
                origin = values[row][2][self.tabPolish_widgets['widget_{}_{}'.format(row, 2)].currentIndex()]
                folder_aligned = values[row][3][self.tabPolish_widgets['widget_{}_{}'.format(row, 3)].currentIndex()]
                template = values[row][4][self.tabPolish_widgets['widget_{}_{}'.format(row, 4)].currentIndex()]

                dimZ = self.tabPolish_widgets['widget_{}_{}'.format(row, 5)].text()
                bin_read = self.tabPolish_widgets['widget_{}_{}'.format(row, 6)].text()
                maxShift = self.tabPolish_widgets['widget_{}_{}'.format(row, 7)].text()
                offx = self.tabPolish_widgets['widget_{}_{}'.format(row, 8)].text()
                offy = self.tabPolish_widgets['widget_{}_{}'.format(row, 9)].text()
                offz = self.tabPolish_widgets['widget_{}_{}'.format(row, 10)].text()
                meta = values[row][11][self.tabPolish_widgets['widget_{}_{}'.format(row, 11)].currentIndex()]
                alignfile = values[row][12][self.tabPolish_widgets['widget_{}_{}'.format(row, 12)].currentIndex()]
                fscfile = values[row][13][self.tabPolish_widgets['widget_{}_{}'.format(row, 13)].currentIndex()]

                fscflag = f'--FSCPath {fscfile}' if fscfile else ''
                metaflag = f'--metaFile {meta}' if meta else ''
                alignflag = f'--alignmentResultsFile {alignfile}' if alignfile else ''

                refid = folder_aligned.split('marker_')[1].split('_')[0]

                # generate output directory for particle polishing results if output directory does not exist.
                destination = os.path.join(f'ParticlePolishing/tomogram_{tomoindex}')
                if not os.path.exists(os.path.join(self.subtomodir, destination)):
                    os.mkdir(os.path.join(self.subtomodir, destination))

                execfilename = os.path.join(self.subtomodir, destination, f'inputParticlePolishing.sh')
                qname, n_nodes, cores, time, modules = self.qparams['BatchParticlePolish'].values()

                paramsCmd = [self.subtomodir, str(cores * n_nodes), self.pytompath, particleXML, folder_aligned,
                             template, destination, bin_read, maxShift, offx, offy, offz, fscflag, metaflag, alignflag]

                txt = polishParticles.format(d=paramsCmd)
                jobtxt = guiFunctions.gen_queue_header(folder=self.logfolder,
                                                       name='PolishParticles_{}'.format(nsj % num_nodes),
                                                       num_jobs_per_node=cores, time=time, partition=qname,
                                                       modules=modules, num_nodes=n_nodes) + txt

                ID, num = self.submitBatchJob(execfilename, pid, jobtxt)
                nsj += 1
                qIDs.append(ID)


        if nsj:
            self.popup_messagebox('Info', 'Submission Status', f'Submitted {nsj} jobs to the queue.')
            self.addProgressBarToStatusBar(qIDs, key='QJobs', job_description='Particle Polish')

    def gen_average(self, params):
        key_particleList, key_filename_average, key_outputDir = params
        particleList = self.widgets[key_particleList].text()
        if not particleList:
            self.popup_messagebox('Warning', 'Averaging Failed', 'Averaging did not succeed. No particle list selected.')
            return
        folder = self.widgets[key_outputDir].text()
        if not os.path.exists(folder): os.mkdir(folder)
        output_name = os.path.join( folder, 'average.em')

        out = os.popen('cd {}; average.py -p {} -a {} -c 4'.format(self.subtomodir, particleList, output_name)).read()
        print(out)
        if not os.path.exists(particleList):
            self.popup_messagebox('Warning', 'Averaging Failed',
                                  'Averaging did not succeed, please try again.')
            return
        self.widgets[key_filename_average].setText(output_name)

    def gen_mask(self,params):
        maskfilename = CreateMaskFile(self, maskfname=params[-1])
        maskfilename.show()

    def selectTemplates(self, key):
        self.selectPLWindow.close()
        try: self.templateList.text()
        except: self.templateList = QLineEdit()

        self.selectTemplatesWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.em', '.mrc'],
                                                 outputline=self.templateList, id=key,
                                                 run_upon_complete=self.selectFSCFilterFiles,
                                                 title='Select templates')

    def selectFSCFilterFiles(self, key):
        self.selectTemplatesWindow.close()
        try: self.fscList.text()
        except: self.fscList = QLineEdit()

        self.selectFSCFilesWindow = SelectFiles(self, initdir=self.pickpartdir, search='file', filter=['.dat'],
                                                 outputline=self.fscList, id=key,
                                                 run_upon_complete=self.populatePartPolishingBatchTable,
                                                 title='Select FSC File for filtering cross-correlation (optional)')

    def updateFSCFlags(self, mode, id='all'):
        for n, (name, flag) in enumerate((('particleList', '--pl'), ('volume1','--v1'), ('volume2', '--v2'),
                                          ('mask', '--mask'))):
            if n != id and id != 'all': continue
            t = self.widgets[mode + name].text()

            if t: self.widgets[mode + 'flag' + name[:1].capitalize() + name[1:]].setText(f'{flag} {t}')
            else: self.widgets[mode + 'flag' + name[:1].capitalize() + name[1:]].setText('')


        di = {'-Even.':'-Odd.', '-Odd.':'-Even.'}
        if id == 1:
            t = self.widgets[mode + f'volume{id}'].text()
            for i in ('-Even.', '-Odd.'):
                if not i in t: continue
                otherHalfMap = t.replace(i,di[i])
                #print(otherHalfMap == self.widgets[mode + f'volume{3-id}'].text(), otherHalfMap, self.widgets[mode + f'volume{3-id}'].text())
                if os.path.exists(otherHalfMap) and self.widgets[mode + f'volume{3-id}'].text() != otherHalfMap:
                    self.widgets[mode + f'volume{3-id}'].setText(otherHalfMap)
                    return

    def updateFSCPlotFlag(self, mode):
        t = self.widgets[mode + 'plot'].isChecked()

        if t:
            self.widgets[mode + 'flagPlot'].setText('--plot ')
            self.widgets[mode + 'queue'].setChecked(False)
        else:
            self.widgets[mode + 'flagPlot'].setText('')

    def updateFSCCombResFlag(self, mode):
        t = self.widgets[mode + 'CombinedResolution'].isChecked()

        if t:
            self.widgets[mode + 'flagCombinedResolution'].setText('--combinedResolution')
        else:
            self.widgets[mode + 'flagCombinedResolution'].setText('')

    def gen_fsc_mask(self,params):
        maskfilename = CreateFSCMaskFile(self, params[-1])
        maskfilename.show()

    def updateLog(self, mode):
        if self.widgets[mode + 'queue'].isChecked():
            self.widgets[mode + 'plot'].setChecked(False)
