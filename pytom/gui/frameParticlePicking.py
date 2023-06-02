import os
import random
import glob
import numpy
from copy import deepcopy

# from ftplib import FTP_TLS, FTP
# from os.path import dirname, basename
# from multiprocessing import Manager, Event, Process

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# pytom gui imports
from pytom.gui.guiStructures import GuiTabWidget, SelectFiles, ParticlePicker, Worker
from pytom.gui.guiStructures import CreateMaskFile, ConvertEM2PDB, SelectTomogramDir
from pytom.gui.guiSupportCommands import templateTM, templateXML, templateExtractCandidates, createParticleList
from pytom.gui import guiFunctions
from pytom.basic.files import read
from pytom.convert.coords2PL import convertCoords2PL
from pytom.bin.updateParticleList import updatePL
from pytom.lib.pytom_numpy import vol2npy


class ParticlePick(GuiTabWidget):
    '''Collect Preprocess Widget'''

    # noinspection PyInterpreter
    def __init__(self, parent=None):
        super(ParticlePick, self).__init__(parent)
        self.stage='v03_'
        self.addGeneralVariables()

        headers = ["Manual Picking","Template Matching", "Create Particle List", "Alter Particle List"]
        subheaders  = [[],['Single', 'Batch Template Match', 'Batch Extract'], ['Single','Batch'], []]
        tabUIs = [[self.tab1UI],[self.tab21UI,self.tab22UI, self.tab23UI],[self.tab31UI,self.tab32UI],[self.tab4UI]]
        static_tabs = [[True],[True,False, False],[True, False],[True]]

        self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders, tabUIs=tabUIs,tabs=self.tabs_dict,
                     tab_actions=self.tab_actions, static_tabs=static_tabs)

    def tab1UI(self,  key=''):
        self.no_image = False
        parent = self.table_layouts[key]
        parent.setAlignment(self, Qt.AlignTop)
        mode = self.stage + 'manualpp_'

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        self.insert_label_line_push(parent, 'Tomogram', mode + 'tomogramFname', initdir=self.tomogramfolder,
                                    tooltip='Select the particle list.', mode='file', filetype=['mrc','em'],cstep=-1)
        self.insert_pushbutton(parent, 'Pick!', tooltip='Select the folder with the aligned tilt images.',
                               rstep=1, cstep=2, action=self.insert_image, params = [parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA,rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB,width=200)

    def tab21UI(self, key=''):

        mode, mode2, mode3 = self.stage + 'SingleTemplateMatch_', self.stage + 'SingleExtractCandidates_', self.stage + 'CreateMask_'

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        items += list(self.create_expandable_group(self.templateMatch, self.sizePolicyB, 'Template Matching',
                                                   mode=mode))
        items[-1].setVisible(False)

        items += list(self.create_expandable_group(self.extractCandidates, self.sizePolicyB, 'Extract Candidates',
                                                   mode=mode2))
        items[-1].setVisible(False)

        # items += list(self.create_expandable_group(self.createTemplateMask, self.sizePolicyB, 'Create Mask Template',
        #                                           mode=mode2))

        # items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab22UI(self, key=''):

        try:
            self.jobFiles.text()
        except:
            self.jobFiles = QLineEdit()

        try:
            self.templateFiles.text()
        except:
            self.templateFiles = QLineEdit()
        self.batchTM = SelectFiles(self, initdir=self.ccfolder, search='file', filter=['em', 'mrc'], id=key,
                                   outputline=self.templateFiles, run_upon_complete=self.getMaskFiles,
                                   title='Select Template')
        # self.batchTM = SelectFiles(self, initdir=self.tomogramfolder, search='file', filter=['em', 'mrc'],
        #                            outputline=self.jobFiles, run_upon_complete=self.getTemplateFiles, id=key,
        #                            title='Select Tomograms.')

        self.tomogramlist = []

        widget = SelectTomogramDir(self)
        widget.show()

    def tab23UI(self, key=''):

        try:
            self.jobFilesExtract.text()
        except:
            self.jobFilesExtract = QLineEdit()

        import glob
        query = os.path.join(self.ccfolder, '*/job*.xml')
        existing_files = self.jobFilesExtract.text().split('\n')
        files = [os.path.join(self.ccfolder, file) for file in sorted(glob.glob(query)) if file and not file in existing_files]
        files = existing_files + files
        if files:
            self.jobFilesExtract.setText('\n'.join(files))

        #else:
            #self.popup_messagebox('Info', 'No files found', 'No job files found from previous template matching jobs.')
        self.batchEC = SelectFiles(self, initdir=self.ccfolder, search='file', filter=['xml'],
                                   outputline=self.jobFilesExtract, run_upon_complete=self.populate_batch_extract_cand,
                                   id=key, title='Select job files.')

    def tab31UI(self, key=''):
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        items += list( self.create_expandable_group(self.createParticleList, self.sizePolicyB, 'Create Particle List',
                                                    mode=self.stage+'partlist_') )
        items[-1].setVisible(False)
        #items += list( self.create_expandable_group(self.createSubtomograms, self.sizePolicyB, 'Extract Subtomograms',
        #                                            mode=self.stage+'extract_') )
        #items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n+1, 0, Qt.AlignRight)

    def tab32UI(self, key=''):

        try:
            self.particleLists.text()
        except:
            self.particleLists = QLineEdit()

        self.b = SelectFiles(self, initdir=self.projectname, search='file', filter=['txt', 'xml'], id=key,
                             outputline=self.particleLists, run_upon_complete=self.populate_batch_create)

    def tab4UI(self,  key=''):

        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        self.modes = [self.stage + 'changeParams_', self.stage + 'extractParticlesFromXML_', self.stage + 'actions_']

        items = []

        items += list(self.create_expandable_group(self.changeXMLParameters, self.sizePolicyB, 'Change Parameters',
                                                   mode=self.modes[0], setVisible=False))
        items += list(
            self.create_expandable_group(self.extractParticles, self.sizePolicyB, 'Extract Particles from XML',
                                         mode=self.modes[1]))
        items += list(self.create_expandable_group(self.actions, self.sizePolicyB, 'Actions', mode=self.modes[2]))

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0)
        self.adjust_items = items
        for i in range(3):
            self.adjust_items[i * 2].stateChanged.connect(
                lambda d, m=self.modes[i], num=i: self.unsetCheckBoxes(m, num))

        pushbutton = QPushButton('Adjust!')
        pushbutton.setSizePolicy(self.sizePolicyB)
        pushbutton.setFixedWidth(100)
        pushbutton.clicked.connect(self.adjustParticleList)
        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(pushbutton, n + 1, 0)
        grid.addWidget(label, n + 2, 0)

    def templateMatch(self, mode, title=''):
        tooltip = 'Run pytom template matching routine.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 25, 20
        self.items = [['', ] * columns, ] * rows
        w=170

        self.ccfolder = os.path.join(self.templatematchfolder,'cross_correlation')

        self.insert_label_line_push(parent, 'Tomogram file', wname=mode+'tomoFname',width=w, initdir=self.tomogramfolder,
                                    tooltip='Select the tomogram file used for template matching.',
                                    filetype=['mrc', 'em'], mode='file')
        self.insert_label_line_push(parent, 'Template file', wname=mode + 'templateFname',width=w,initdir=self.ccfolder,
                                    tooltip='Select the tomogram file used for template matching.',
                                    filetype=['em', 'mrc'], mode='file', rstep=0,cstep=1)
        self.insert_pushbutton(parent, 'Create', action=self.pdb2em, params=[mode + 'templateFname', mode], cstep=-3, rstep=1)
        self.insert_label_line_push(parent, 'Mask file', wname=mode+'maskFname',width=w,initdir=self.ccfolder,
                                    tooltip='Select the tomogram file used for template matching.',
                                    filetype=['em', 'mrc'], mode='file',cstep=1,rstep=0)
        self.insert_pushbutton(parent,'Create',action=self.create_maskfile,params=[mode+'maskFname'],cstep=-3,rstep=1)
        self.insert_label_combobox(parent, 'Spherical mask', mode + 'sphericalMask', labels=['True', 'False'],
                                   tooltip='Set to true if the mask is spherical to save calculation time.')
        self.insert_label_spinbox(parent, mode + 'startZ', 'Start index of object in Z-dimension', width=w,
                                  value=0, stepsize=1, minimum=0, maximum=1000,
                                  tooltip='The tomogram might consist of empty space. Please determine the z-index where the object starts.')
        self.insert_label_spinbox(parent, mode + 'endZ', 'End index of object in Z-dimension', width=w,
                                  value=0, stepsize=1, minimum=0, maximum=1000,
                                  tooltip='The tomogram might consist of empty space. Please determine the z-index where the object ends.')

        self.insert_label_spinbox(parent, mode + 'Wedge1', 'Wedge Angle1 (degrees)',width=w,
                                  value=30, stepsize=1,minimum=0, maximum=90,
                                  tooltip='Angle between 90 and the highest tilt angle.')
        self.insert_label_spinbox(parent,  mode + 'Wedge2', 'Wedge Angle2 (degrees)',width=w, cstep=-1,
                                  value=30, stepsize=1, minimum=0, maximum=90,
                                  tooltip='Angle between -90 and the lowest tilt angle.')
        self.insert_label_combobox(parent,'Angular Sampling Filename',mode+'angleFname',
                                   labels=os.listdir(os.path.join( self.parent().pytompath, 'angles/angleLists') ),
                                   tooltip='Select the file that describes the angular sampling of the template model',
                                   width=w, cstep=-1)
        self.insert_label_spinbox(parent, mode + 'splitX', 'x split', width=w, value=1, stepsize=1,
                                  minimum=1, maximum=100,
                                  tooltip='Split the volume along the x dimensions this many times. '
                                          'Total subvolumes are determined by splitx * splity * splitz. '
                                          'Each subvolume will be assigned to an mpi process. So if you have 16 '
                                          'subvolumes, you will need 16 cores on your machine to run the job. '
                                          'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                                          'split the tomogram into 4 subvolumes. Then setup can have the following '
                                          'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3')
        self.insert_label_spinbox(parent, mode + 'splitY', 'y split', width=w, value=1, stepsize=1,
                                  minimum=1, maximum=100,
                                  tooltip='Split the volume along the y dimensions this many times. '
                                          'Total subvolumes are determined by splitx * splity * splitz. '
                                          'Each subvolume will be assigned to an mpi process. So if you have 16 '
                                          'subvolumes, you will need 16 cores on your machine to run the job. '
                                          'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                                          'split the tomogram into 4 subvolumes. Then setup can have the following '
                                          'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3')
        self.insert_label_spinbox(parent, mode + 'splitZ', 'z split', width=w, value=1, stepsize=1,
                                  minimum=1, maximum=100,
                                  tooltip='Split the volume along the z dimensions this many times. '
                                          'Total subvolumes are determined by splitx * splity * splitz. '
                                          'Each subvolume will be assigned to an mpi process. So if you have 16 '
                                          'subvolumes, you will need 16 cores on your machine to run the job. '
                                          'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                                          'split the tomogram into 4 subvolumes. Then setup can have the following '
                                          'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3')
        self.insert_label_line(parent, "GPU's", mode + 'gpuID', width=w, cstep=0,
                                  tooltip="Which GPU's do you want to reserve. If you want to use multiple GPUs separate them using a comma, e.g. 0,1,2 ")

        self.widgets[mode + 'widthZ'] = QLineEdit('0')
        self.widgets[mode + 'widthX'] = QLineEdit('0')
        self.widgets[mode + 'widthY'] = QLineEdit('0')
        self.widgets[mode + 'gpuString'] = QLineEdit('')
        self.widgets[mode + 'jobName'] = QLineEdit()
        self.widgets[mode + 'outfolderTM'] = QLineEdit(self.ccfolder)
        self.widgets[mode + 'numCores'] = QLineEdit(str(self.widgets[mode + 'splitX'].value() *
                                                    self.widgets[mode + 'splitY'].value() *
                                                    self.widgets[mode + 'splitZ'].value()))

        self.widgets[mode + 'tomoFname'].textChanged.connect(lambda d, m=mode: self.updateTM(m))
        self.widgets[mode + 'startZ'].valueChanged.connect(lambda d, m=mode: self.updateZWidth(m))
        self.widgets[mode + 'endZ'].valueChanged.connect(lambda d, m=mode: self.updateZWidth(m))
        self.widgets[mode + 'templateFname'].textChanged.connect(lambda d, m=mode: self.updateJobName(m))
        self.widgets[mode + 'gpuID'].textChanged.connect(lambda d, m=mode: self.updateGpuString(m))

        self.execfilenameTM = os.path.join(self.templatematchfolder, 'templateMatch.sh')
        self.xmlfilename  = os.path.join(self.templatematchfolder, 'job.xml')

        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'TemplateMatching'
        paramsSbatch['folder'] = self.logfolder
        paramsSbatch['id'] = 'SingleTemplateMatch'
        paramsSbatch['gpu'] = self.widgets[mode + 'gpuID']

        mandatory_fill = [mode + 'templateFname', mode + 'tomoFname', mode+'maskFname']

        self.updateTM(mode)

        self.insert_gen_text_exe(parent, mode, jobfield=True, exefilename=[mode+'outfolderTM','templateMatch.sh'], paramsSbatch=paramsSbatch,
                                 paramsXML=[mode + 'tomoFname', mode + 'templateFname', mode+'maskFname', mode + 'Wedge1',
                                            mode + 'Wedge2', mode+'angleFname', mode + 'outfolderTM', mode + 'startZ',
                                            mode + 'widthX', mode + 'widthY', mode + 'widthZ',
                                            mode + 'sphericalMask', templateXML],
                                 paramsCmd=[mode + 'outfolderTM', mode + 'numCores', self.pytompath, mode + 'jobName',
                                            mode + 'splitX', mode + 'splitY', mode + 'splitZ',
                                            mode + 'gpuString', templateTM],
                                 xmlfilename=[mode+'outfolderTM', mode + 'jobName'], mandatory_fill=mandatory_fill)

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA, rstep=1)

        self.updateTM(mode)

        setattr(self, mode + 'gb_TMatch', groupbox)
        return groupbox

    def extractCandidates(self, mode, title=''):
        tooltip = 'Run pytom extractCandidates.py'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)
        self.row, self.column = 0, 2
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w=170

        self.insert_label_line_push(parent, 'Job File', wname=mode + 'jobXML', width=w,
                                    tooltip='Localization job XML file.',
                                    filetype='xml', mode='file')
        self.insert_label_line_push(parent, 'Score File', wname=mode + 'scoreFile', width=w,
                                    tooltip='File with score coefficients (score.em).',
                                    filetype=['em','mrc'], mode='file')
        self.insert_label_line_push(parent, 'Angles File', wname=mode + 'anglesFile', width=w,
                                    tooltip='File with orientation indices (angles.em).',
                                    filetype=['em', 'mrc'], mode='file')
        self.insert_label_line_push(parent, 'Mask File', wname=mode + 'maskFile', width=w,
                                    tooltip='File with selected area used for finding hits.',
                                    filetype=['em', 'mrc'], mode='file')

        self.insert_label_line_push(parent, 'File Name Particle List', wname=mode + 'particleList', width=w,
                                    tooltip='Select the angles.em file produced by template matching.',
                                    filetype=['xml'], mode='file')
        self.insert_label_line(parent, 'Prefix', mode + 'particlePath', width=w,
                               tooltip='Path prepended to each particle.')

        self.insert_label_spinbox(parent, mode + 'Size', 'Radius particle (px)', width=w, cstep=-1,
                                  value=10, stepsize=1, minimum=0, maximum=90,
                                  tooltip='Radius around potential candidate that will be ignored during further processing.')
        self.insert_label_spinbox(parent, mode + 'NumberOfCandidates', 'Number of Candidates', width=w, cstep=-1,
                                  value=1000, stepsize=100, minimum=0, maximum=5000,
                                  tooltip='Number of candidates to extract.')
        self.insert_label_spinbox(parent, mode + 'MinimalScoreValue', 'Minimal Score Value', width=w, cstep=0,
                                  value=0., stepsize=0.1, minimum=0, maximum=1.,wtype=QDoubleSpinBox,
                                  tooltip='Minimal score value to which to extract.')

        self.widgets[mode + 'jobXML'].textChanged.connect(lambda d, m=mode: self.jobXMLChanged(m))
        self.widgets[mode + 'particleList'].textChanged.connect(lambda d, m=mode: self.particleListChanged(m))
        self.jobXMLChanged(mode)

        self.widgets[mode + 'maskGenerate'] = QLineEdit()
        self.widgets[mode + 'maskFile'].textChanged.connect(lambda d, m=mode: self.updateMaskGenerate(m))


        execfilename = os.path.join( self.templatematchfolder, 'extractCandidates.sh')
        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'ExtractCandidates'
        paramsSbatch[ 'folder' ] = self.logfolder #self.templatematchfolder
        #paramsSbatch['partition'] = 'fastq'
        paramsSbatch['time'] = 1
        paramsSbatch['num_jobs_per_node'] = 1

        paramsCmd=[self.templatematchfolder, self.pytompath, mode+'jobXML', mode+'scoreFile', mode+'anglesFile',
                   mode+'particleList', mode+'particlePath', mode+'Size', mode+'NumberOfCandidates',
                   mode+'MinimalScoreValue', mode+'maskGenerate', templateExtractCandidates]

        mandatory_fill = [mode + 'jobXML', mode + 'scoreFile', mode + 'anglesFile', mode + 'particleList',mode + 'particlePath' ]

        self.insert_gen_text_exe(parent, mode, jobfield=False, exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 paramsCmd=paramsCmd, action=self.activate_stage, paramsAction=3,
                                 mandatory_fill=mandatory_fill)

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA,rstep=1)

        setattr(self, mode + 'gb_extractCandidates', groupbox)
        return groupbox

    def populate_batch_templatematch(self, id='tab22'):
        print('multiple template matching job-submissions')
        self.batchTM.close()
        templateFiles = sorted(self.templateFiles.text().split('\n'))
        maskFiles     = sorted(self.maskFiles.text().split('\n'))

        if len(maskFiles) == 0 or len(templateFiles) == 0 or len(self.tomogramlist) == 0:
            print('\n\nPlease select at least one tomogram, template and mask file.\n\n')
            return

        headers = ['Filename Tomogram', 'Run', 'Mirrored', 'Optional Templates', 'Optional Masks', 'Spherical mask',
                   'Wedge Angle 1', 'Wedge Angle 2', 'Angle List', 'Start Z', 'End Z', 'x split', 'y split',
                   'z split', 'GPU ID', '']
        types = ['txt', 'checkbox', 'checkbox', 'combobox', 'combobox', 'combobox', 'lineedit', 'lineedit',
                 'combobox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'txt']
        sizes = [0, 0, 0, 80, 80, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        tooltip = ['Name of tomogram files.',
                   'Check this box if you want to do template matching using the optional settings.',
                   'Check this box if you want to do template matching using a mirrored version of the template. '
                   'If both run and mirrored are checked template matching will be deployed twice for that tomogram.',
                   'Optional templates.',
                   'Optional masks.',
                   'Set to true if the mask is spherical to save calculation time.',
                   'Angle between the lowest tilt angle and -90.',
                   'Angle between 90 and the highest tilt angle.',
                   'Optional angle lists.',
                   'At which Z-slice does your object start?', 'At which Z-slice does your object end?',
                   'Split the volume along the x dimensions this many times. '
                   'Total subvolumes are determined by splitx * splity * splitz. '
                   'Each subvolume will be assigned to an mpi process. So if you have 16 '
                   'subvolumes, you will need 16 cores on your machine to run the job. '
                   'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                   'split the tomogram into 4 subvolumes. Then setup can have the following '
                   'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3',
                   'Split the volume along the y dimensions this many times. '
                   'Total subvolumes are determined by splitx * splity * splitz. '
                   'Each subvolume will be assigned to an mpi process. So if you have 16 '
                   'subvolumes, you will need 16 cores on your machine to run the job. '
                   'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                   'split the tomogram into 4 subvolumes. Then setup can have the following '
                   'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3',
                   'Split the volume along the z dimensions this many times. '
                   'Total subvolumes are determined by splitx * splity * splitz. '
                   'Each subvolume will be assigned to an mpi process. So if you have 16 '
                   'subvolumes, you will need 16 cores on your machine to run the job. '
                   'If you have 4 gpus, each managed by 1 mpi proc, you can maximally '
                   'split the tomogram into 4 subvolumes. Then setup can have the following '
                   'parameters: splitx=2, splity=2, splitz=1, gpus=0,1,2,3',
                   'GPU ID']

        values = []

        angleLists = os.listdir(os.path.join(self.pytompath, 'angles/angleLists'))
        for n, tomogramFile in enumerate(self.tomogramlist):
            try:
                folder = os.path.dirname(os.popen(f'ls -alrt {tomogramFile}').read()[:-1].split(' ')[-1])
                widthfile = os.path.join(folder, 'z_limits.txt')
                if os.path.exists( widthfile):

                    start, end = map(int, list(numpy.loadtxt(widthfile)))
                else:
                    start, end = 0,0

                folderm = os.path.join(self.tomogram_folder,
                                       os.path.basename(tomogramFile)[:len('tomogram_000')], 'sorted')
                metafiles = [os.path.join(folderm, dir) for dir in os.listdir(folderm)
                             if not os.path.isdir(dir) and dir.endswith('.meta')] if os.path.exists(folderm) else []
                if metafiles:
                    from pytom.gui.guiFunctions import datatype, loadstar

                    metadata = loadstar(metafiles[0],dtype=datatype)
                    tiltAngles = metadata['TiltAngle']
                    w1 = int(numpy.round(tiltAngles.min()+90))
                    w2 = int(90-numpy.round(tiltAngles.max()))
                else:
                    w1,w2 = 30,30

            except Exception as e:
                print(e)
                start, end = 0, 0
                w1,w2 = 30,30
            values.append([tomogramFile, 1, 1, templateFiles, maskFiles, ['True', 'False'], w1, w2,
                           angleLists, start, end, 1, 1, 1, '', ''])

        try:
            self.num_nodes[id].setParent(None)
        except:
            pass

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, nn=True, wname=self.stage + 'BatchTemplateMatch_')

        self.tab22_widgets = self.tables[id].widgets

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.mass_submitTM(pid, v))

    def mass_submitTM(self, pid, values):
        num_nodes = int(self.num_nodes[pid].value())
        num_submitted_jobs = 0
        qIDs = []
        jobCode = {}
        execfilenames = {}

        for row in range(self.tables[pid].table.rowCount()):
            try:
                normal = self.tab22_widgets['widget_{}_{}'.format(row, 1)].isChecked()
                mirrored = self.tab22_widgets['widget_{}_{}'.format(row, 2)].isChecked()
                if not normal and not mirrored:
                    continue

                tomogramFile = values[row][0]
                templateFile = values[row][3][self.tab22_widgets['widget_{}_{}'.format(row, 3)].currentIndex()]
                maskFile = values[row][4][self.tab22_widgets['widget_{}_{}'.format(row, 4)].currentIndex()]
                mask_is_spherical = values[row][5][self.tab22_widgets['widget_{}_{}'.format(row, 5)].currentIndex()]
                w1 = float(self.tab22_widgets['widget_{}_{}'.format(row, 6)].text())
                w2 = float(self.tab22_widgets['widget_{}_{}'.format(row, 7)].text())
                angleList = self.tab22_widgets['widget_{}_{}'.format(row, 8)].currentText()
                start = self.tab22_widgets['widget_{}_{}'.format(row, 9)].text()
                end = self.tab22_widgets['widget_{}_{}'.format(row, 10)].text()
                splitx = int(self.tab22_widgets['widget_{}_{}'.format(row, 11)].text())
                splity = int(self.tab22_widgets['widget_{}_{}'.format(row, 12)].text())
                splitz = int(self.tab22_widgets['widget_{}_{}'.format(row, 13)].text())
                gpuID = self.tab22_widgets['widget_{}_{}'.format(row, 14)].text()

                if gpuID:
                    gpuIDFlag = f'--gpuID {gpuID}'
                else:
                    gpuIDFlag = ''

                widthX, widthY, widthZ = 0, 0, 0
                try:
                    start, end = int(start), int(end)
                    widthZ = end - start
                    if widthZ < 1:
                        widthZ = 0
                        start = 0
                    if widthZ:
                        from pytom.lib.pytom_volume import read
                        v = read(tomogramFile)
                        widthX = v.sizeX()
                        widthY = v.sizeY()
                        del v

                except:
                    continue

                tomofile, ext = os.path.splitext(tomogramFile)
                suffices = []
                if normal:
                    suffices.append('')
                if mirrored:
                    suffices.append('_Mirrored')
                    outDirectory = os.path.join(self.ccfolder, os.path.basename(tomofile))
                    if not os.path.exists(outDirectory): os.mkdir(outDirectory)

                    from pytom.lib.pytom_volume import read, vol, mirrorVolume
                    for n, input_filename in enumerate([templateFile, maskFile]):
                        v = read(input_filename)
                        res = vol(v)
                        mirrorVolume(v, res)

                        i = os.path.basename(input_filename)
                        base, ext = os.path.splitext(i)
                        i = base + '_Mirrored' + ext
                        output_filename = os.path.join(outDirectory, i)
                        if n == 0: templateFileMirrored = output_filename
                        if n == 1: maskFileMirrored = output_filename

                        res.write(output_filename)

                for n, suffix in enumerate(suffices):
                    outDirectory = os.path.join(self.ccfolder, os.path.basename(tomofile))
                    if not os.path.exists(outDirectory): os.mkdir(outDirectory)
                    if 'Mirror' in suffix:
                        templateFile = templateFileMirrored
                        maskFile = maskFileMirrored
                        mm = 'Mirrored'
                    else:
                        mm = ''

                    XMLParams = [tomogramFile, templateFile, maskFile, w1, w2, angleList, outDirectory,
                                 start, widthX, widthY, widthZ, mask_is_spherical]

                    jobxml = templateXML.format(d=XMLParams)
                    template = os.path.basename(templateFile).split('.')[0]
                    jobname = 'job_{}.xml'.format(template)
                    outjob = open(os.path.join(outDirectory, jobname), 'w')
                    outjob.write(jobxml)
                    outjob.close()

                    # check number of splits is not larger than number of gpus
                    if gpuIDFlag and len(list(map(int, gpuID.split(',')))) < splitx * splity * splitz:
                        self.popup_messagebox('Warning', 'Parameter issue',
                                              'Number of gpus is less than the number of subvolumes a '
                                              'tomogram is split into, this will crash.')
                        return

                    fname = 'TM_Batch_ID_{}'.format(num_submitted_jobs % num_nodes)
                    folder = outDirectory
                    cmd = templateTM.format(d=[outDirectory, splitx * splity * splitz, self.pytompath, jobname,
                                               splitx, splity, splitz, gpuIDFlag])
                    suffix = "_" + os.path.basename(outDirectory)
                    qname, n_nodes, cores, timed, modules, qcmd = self.qparams['BatchTemplateMatch'].values()

                    # TODO this should have some extra checks for combi of queue subm and gpus
                    if gpuID and not gpuID in jobCode.keys():
                        job = guiFunctions.gen_queue_header(folder=self.logfolder, name=fname, suffix=suffix,
                                                            time=timed, num_nodes=n_nodes, partition=qname,
                                                            modules=modules, cmd=qcmd,
                                                            num_jobs_per_node=cores, gpus=gpuID) * self.checkbox[
                                  pid].isChecked() + cmd
                        execfilenames[gpuID] = os.path.join(self.templatematchfolder,
                                                            f'templateMatchingBatch_{num_submitted_jobs}.sh')
                        jobCode[gpuID] = job

                        outjob = open(os.path.join(outDirectory, f'templateMatchingBatch{mm}.sh'), 'w')
                        outjob.write(job)
                        outjob.close()
                        num_submitted_jobs += 1

                    elif gpuID and gpuID in jobCode.keys():

                        jobCode[gpuID] += f'\nwait\n\n{cmd}\n'

                        job = guiFunctions.gen_queue_header(folder=self.logfolder, name=fname, suffix=suffix,
                                                            time=timed, num_nodes=n_nodes, partition=qname,
                                                            modules=modules, cmd=qcmd,
                                                            num_jobs_per_node=cores, gpus=gpuID) * self.checkbox[
                                  pid].isChecked() + cmd

                        outjob = open(os.path.join(outDirectory, f'templateMatchingBatch{mm}.sh'), 'w')
                        outjob.write(job)
                        outjob.close()
                        num_submitted_jobs += 1

                    elif not f'noGPU_{num_submitted_jobs % num_nodes}' in jobCode.keys():
                        job = guiFunctions.gen_queue_header(folder=self.logfolder, name=fname, suffix=suffix,
                                                            time=timed, num_nodes=n_nodes, partition=qname, cmd=qcmd,
                                                            modules=modules, num_jobs_per_node=cores) + cmd
                        execfilenames[f'noGPU_{num_submitted_jobs % num_nodes}'] = os.path.join(
                            self.templatematchfolder, f'templateMatchingBatch_{num_submitted_jobs}.sh')
                        jobCode[f'noGPU_{num_submitted_jobs % num_nodes}'] = job
                        outjob = open(os.path.join(outDirectory, f'templateMatchingBatch{mm}.sh'), 'w')
                        outjob.write(job)
                        outjob.close()
                        num_submitted_jobs += 1

                    else:
                        jobCode[f'noGPU_{num_submitted_jobs % num_nodes}'] += f'\nwait\n\n{cmd}\n'
                        job = guiFunctions.gen_queue_header(folder=self.logfolder, name=fname, suffix=suffix,
                                                            time=timed, num_nodes=n_nodes, partition=qname,
                                                            modules=modules, cmd=qcmd,
                                                            num_jobs_per_node=cores, gpus=gpuID) * self.checkbox[
                                  pid].isChecked() + cmd

                        outjob = open(os.path.join(outDirectory, f'templateMatchingBatch{mm}.sh'), 'w')
                        outjob.write(job)
                        outjob.close()
                        num_submitted_jobs += 1

            except Exception as e:
                print(e)

        qIDs, num_submitted_jobs = [], 0

        for key, values in jobCode.items():
            ID, num = self.submitBatchJob(execfilenames[key], pid, values)
            qIDs.append(ID)
            num_submitted_jobs += 1

        if num_submitted_jobs > 0:
            self.popup_messagebox('Info', 'Submission Status', f'Submitted {num_submitted_jobs} jobs to the queue.')
            self.addProgressBarToStatusBar(qIDs, key='QJobs', job_description='Templ. Match',
                                           num_submitted_jobs=num_submitted_jobs)
        else:
            self.popup_messagebox('Info', 'Submission Status', 'Failed at starting template matching jobs. Check if '
                                                               'tomograms exist.')

    def populate_batch_extract_cand(self,id='tab23'):
        self.batchEC.close()
        jobfiles = sorted(self.jobFilesExtract.text().split('\n'))

        jobfiles += glob.glob(os.path.join(self.ccfolder, '*/job*flipped.xml'))
        jobfiles = numpy.unique(jobfiles)

        if len(jobfiles) == 0:
            print('\n\nPlease select at least one job file template and mask file.\n\n')
            return

        headers = ["Job name", "Extract", "File Name Particle List", "Ouput Dir Subtomograms", "Radius (px)",
                   "# Candidates", 'Min. Score', '']
        types = ['txt', 'checkbox', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'txt']
        sizes = [0, 0, 80, 150, 0, 0, 0, 0]

        tooltip = ['Name of job files.',
                   'Select if you want to extract particles using this job file.',
                   'File name of particle list in which a number of candidates are written.',
                   'Prefix indicating in which folder particles will be saved when subtomograms are extracted',
                   'Radius of Particle of interest in pixels',
                   'Number of candidates extracted from tomogram.',
                   'Minimum cross-correlation coefficient for a particle to be selected.']

        values = []

        for n, jobFile in enumerate(jobfiles):
            print(jobFile)
            if not jobFile: continue
            folder = os.path.basename(os.path.dirname(jobFile))
            suffix = '_' + os.path.basename(jobFile)[4:-4]

            particleList = os.path.join(self.pickpartfolder, f'particleList_TM_{folder}{suffix}.xml')
            if not os.path.exists(os.path.dirname(particleList)): os.mkdir(os.path.dirname(particleList))
            p = os.path.basename(particleList)[:-4]
            prefix = f'Subtomograms/particleList_TM_{folder}'
            values.append([jobFile, True, particleList, prefix, 16, 1000, 0.001, ''])

        try:
            self.num_nodes[id].setParent(None)
        except:
            pass

        self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip, nn=True, wname=self.stage + 'BatchExtractCandidates_')

        self.tab23_widgets = self.tables[id].widgets

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.mass_submitEC(pid, v))

    def mass_submitEC(self, pid, values):
        num_nodes = int(self.num_nodes[pid].value())
        num_submitted_jobs = 0
        qIDs = []

        for row in range(self.tables[pid].table.rowCount()):
            try:
                if self.tab23_widgets['widget_{}_{}'.format(row, 1)].isChecked():
                    jobFile       = values[row][0]
                    suffix        = os.path.basename(jobFile)[3:-4]
                    scoresFile    = os.path.join(os.path.dirname(jobFile), f'scores{suffix}.em')
                    anglesFile    = os.path.join(os.path.dirname(jobFile), f'angles{suffix}.em')
                    particleList  = os.path.join(self.pickpartfolder, self.tab23_widgets['widget_{}_{}'.format(row, 2)].text())
                    particlePath  = self.tab23_widgets['widget_{}_{}'.format(row, 3)].text()
                    particleSize  = self.tab23_widgets['widget_{}_{}'.format(row, 4)].text()
                    numCandidates = self.tab23_widgets['widget_{}_{}'.format(row, 5)].text()
                    minCCScore    = self.tab23_widgets['widget_{}_{}'.format(row, 6)].text()

                    paramsCmd     = [self.templatematchfolder, self.pytompath, jobFile, scoresFile, anglesFile,
                                    particleList, particlePath, particleSize, numCandidates, minCCScore, '']

                    fname         = 'EC_Batch_ID_{}'.format(num_submitted_jobs % num_nodes)
                    cmd           = templateExtractCandidates.format(d=paramsCmd)
                    qname, n_nodes, cores, time, modules, qcmd = self.qparams['BatchExtractCandidates'].values()
                    job           = guiFunctions.gen_queue_header(folder=self.logfolder, name=fname, singleton=True,
                                                                  time=time, num_nodes=n_nodes, partition=qname,
                                                                  modules=modules, num_jobs_per_node=cores,
                                                                  cmd=qcmd) * self.checkbox[pid].isChecked() + cmd

                    execfilename = os.path.join(os.path.dirname(jobFile), f'extractCandidatesBatch_{self.ECCounter}.sh')
                    ID, num = self.submitBatchJob(execfilename, pid, job)
                    num_submitted_jobs += 1
                    self.ECCounter += 1
                    qIDs.append(ID)

            except Exception as e:
                print('Error: ', e)
                print(f'Failed to setup job for {os.path.basename(values[row][0])}')

        if num_submitted_jobs:
            self.activate_stage(3)

            self.popup_messagebox('Info', 'Submission Status', f'Submitted {num_submitted_jobs} jobs to the queue.')
            self.addProgressBarToStatusBar(qIDs, key='QJobs', job_description='Extract Cand. Batch')


    def createParticleList(self, mode='', title=''):
        tooltip = 'Tick this box to create a particle list from the selected coordinate file.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 2
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_label(parent, cstep=1, sizepolicy=self.sizePolicyB, width=200)

        self.insert_label_line_push(parent, 'Coordinate File', mode + 'CoordinateFile',mode='file',filetype='txt',
                                    initdir=self.pickpartfolder,
                                    tooltip='Select the coordinate file which contains the x,y,z coordinates of the particles.')

        self.insert_label_line(parent, 'Prefix Subtomograms', mode + 'PrefixSubtomo', value='Subtomograms/',
                               tooltip='Define a prefix that is prepended to the particle name. A unique identifier makes it'+
                               'easier to track the respective particles.')

        self.insert_label_spinbox(parent, mode + 'Wedge1', 'Wedge degrees(positive)',
                               'Angle between 90 and the highest tilt angle.', stepsize=1, value=30, minimum=0, maximum=90)

        self.insert_label_spinbox(parent, mode + 'Wedge2','Wedge degrees(negative)',
                                  'Angle between -90 and the lowest tilt angle.',
                                  stepsize=1,minimum=0,maximum=90,value=30, cstep=-1,rstep=1)

        self.insert_label_line_push(parent, 'Filename particleList', mode+'FnameParticleList', cstep=-2, rstep=1,
                                    pushtext='Browse', width=150, action=self.setFnameParticleList)

        self.insert_label_checkbox(parent, mode+'Randomize', 'Randomize Particle Orientation', cstep=0, rstep=1,
                                   tooltip='Randomize the orientation angles of the particles in the particle list.')

        self.widgets[mode+'flagRandomize'] = QLineEdit()
        self.widgets[mode+'Randomize'].stateChanged.connect(lambda dummy, m=mode: self.update_flag(m) )
        self.widgets[ mode + 'CoordinateFile'].textChanged.connect(lambda d, mm=mode: self.update_prefix_pL(mm))


        execfilename = os.path.join( self.pickpartfolder, 'createParticleList.sh')
        paramsSbatch = guiFunctions.createGenericDict()
        paramsSbatch['fname'] = 'createPrtclLst'
        paramsSbatch['folder'] = self.logfolder #self.pickpartfolder

        paramsCmd=[mode+'CoordinateFile', mode+'PrefixSubtomo', mode+'Wedge1', mode+'Wedge2', mode+'FnameParticleList',
                   mode+'flagRandomize', createParticleList]

        self.insert_gen_text_exe(parent,mode,paramsCmd=paramsCmd, exefilename=execfilename,paramsSbatch=paramsSbatch,
                                 action=self.createFolder,paramsAction=mode)

        setattr(self, mode + 'gb_create_particle_list', groupbox)
        return groupbox

    def populate_batch_create(self, id='tab32'):
        self.b.close()
        coordinateFiles = sorted( self.particleLists.text().split('\n') )

        headers = ["Filename Coordinate List", "Prefix Subtomograms", 'Wedge Angle 1', 'Wedge Angle 2', "Filename Particle List", 'Randomize Angles', '']
        types = ['txt', 'lineedit', 'lineedit', 'lineedit', 'lineedit', 'checkbox', 'txt']
        sizes = [0, 420, 80, 0, 420, 0]


        tooltip = ['Name of coordinate files',
                   'Prefix used for subtomograms',
                   'Angle between 90 and the highest tilt angle.',
                   'Angle between -90 and the lowest tilt angle.',
                   'Filename of generate particle list file (xml)',
                   'Randomize angles in particleList']

        values = []

        for n, coordinateFile in enumerate( coordinateFiles ):
            if not coordinateFile: continue
            outfolder = os.path.join( self.projectname, '04_Particle_Picking/Picked_Particles' )

            try:
                tomogramNUM = int(coordinateFile.split('tomogram_')[-1][:3])
            except:
                tomogramNUM = n


            prefix = 'Subtomograms/{}/particle_'.format(os.path.basename(coordinateFile[:-4]))
            ff = os.path.join(self.subtomofolder,os.path.dirname(prefix))
            if not os.path.exists(ff): os.mkdir(ff)
            fname_plist = 'particleList_{}.xml'.format(os.path.basename(coordinateFile[:-4]))
            active = True
            #if coordinateFile.endswith('.xml'):
                #prefix = 'UNCHANGED'
                #fname_plist = 'UNCHANGED'
                #active=False
            values.append( [coordinateFile, prefix, 30, 30, fname_plist, active, ''] )

        if values: self.fill_tab(id, headers, types, values, sizes, tooltip=tooltip)
        else: return
        self.tab32_widgets = tw = self.tables[id].widgets

        for row in range(self.tables[id].table.rowCount()):
            val = values[row][1]
            widget1 = self.tab32_widgets['widget_{}_{}'.format(row,1)]
            widget2 = self.tab32_widgets['widget_{}_{}'.format(row,4)]
            widget1.textChanged.connect(lambda d, w1=widget1, w2=widget2: self.update_change_tab32_table(w1, w2))

        self.pbs[id].clicked.connect(lambda dummy, pid=id, v=values: self.mass_convert_txt2xml(pid, v))

    def mass_convert_txt2xml(self,pid,values):
        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.pickpartfolder, filter='*.xml')[0])
        if not fname: return
        if not fname.endswith('.xml'):fname+='.xml'
        randomize = False
        AL = False
        conf = [[],[],[],[],[]]

        fnamesPL = []
        wedges = ''
        angleListDefault = os.path.join(self.pytompath, 'angles/angleLists/angles_18_3040.em')
        vol = read(angleListDefault)
        angleListDefaultData = deepcopy(vol2npy(vol)).T.astype(numpy.float32)
        for row in range(self.tables[pid].table.rowCount()):
            try:
                c = values[row][0]
                if c.endswith('.xml'):
                    prefix  = self.tab32_widgets['widget_{}_{}'.format(row, 1)].text()
                    w1      = self.tab32_widgets['widget_{}_{}'.format(row, 2)].text()
                    w2      = self.tab32_widgets['widget_{}_{}'.format(row, 3)].text()
                    outname = self.tab32_widgets['widget_{}_{}'.format(row, 4)].text()
                    try:
                        r   = self.tab32_widgets['widget_{}_{}'.format(row, 5)].isChecked()
                    except:
                        r   = None

                    outname = os.path.join(os.path.dirname(c), outname)
                    wedges = '{},{},'.format(w1,w2)
                    angleListData = angleListDefaultData if r else None
                    updatePL([c], [outname], wedgeangles=[w1,w2], anglelist=angleListData)
                    fnamesPL.append(outname)
                    continue

                p = self.tab32_widgets['widget_{}_{}'.format(row, 1)].text()
                w1 = float(self.tab32_widgets['widget_{}_{}'.format(row,2)].text() )
                w2 = float(self.tab32_widgets['widget_{}_{}'.format(row,3)].text() )
                pl = self.tab32_widgets['widget_{}_{}'.format(row,4)].text()
                pl = os.path.join(self.pickpartfolder, pl)
                r = self.tab32_widgets['widget_{}_{}'.format(row,5)].isChecked()
                #print(createParticleList.format(d=[c, p, w1, w2, pl]))
                wedge = [w1,w2]
                for n, inp in enumerate((c, pl, p, w1, w2)):
                    if n==4: n=3
                    conf[n].append(inp)
                pFolder = os.path.join(self.subtomofolder, os.path.basename(p))
                if not os.path.exists(pFolder): os.mkdir(pFolder)

                angleList = angleListDefault*r


                if angleList:
                    if not os.path.exists(angleList):
                        raise Exception('Angle List is not existing.')
                        angleList = None
                    elif not angleList.split('.')[-1] in ('em', 'mrc'):
                        raise Exception('File format should be mrc or em.')
                        angleList = None
                    else:
                        vol = read(angleList)
                        angleList = deepcopy(vol2npy(vol)).T.astype(numpy.float32)
                        rotation = random.choice(angleList)
                        if len(rotation) != 3 or type(rotation[0]) != numpy.float32:
                            raise Exception('AngleList should contain three floats per row.')
                            angleList = None

                if r and not randomize:
                    randomize = True
                    AL = deepcopy(angleList)

                convertCoords2PL([c], pl, subtomoPrefix=[p], wedge_angles=wedge, angleList=angleList, projDir=self.projectname)
                #os.system(createParticleList.format(d=[c, p, wedge, pl]))

            except Exception as e:
                print(e)
                print('Writing {} failed.'.format(os.path.basename(fname)))
                return


        convertCoords2PL(conf[0], fname, subtomoPrefix=conf[2], wedge_angles=conf[3], angleList=AL, projDir=self.projectname)

        fnamesPL2 = fnamesPL + [fname]

        if wedges: wedges = wedges[:-1]

        if len(fnamesPL2) > 1:
            os.system('combineParticleLists.py -f {} -o {} '.format(",".join(fnamesPL2), fname))
            for remove_fname in fnamesPL:
                os.remove(remove_fname)


    def changeXMLParameters(self, mode='', title=''):
        tooltip = 'Tick this box to update specific fields in an particle list (XML format).'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.insert_label_line_push(parent, 'Particle List', 'particleList0', 'Select your particleList (XML) file.',
                                    initdir=self.projectname, mode='file', cstep=-4, filetype=['xml'],
                                    action2=self.countNumberOfParticles)
        self.insert_checkbox_label_line(parent, mode+'adjustSuffix', 'Suffix', mode+'suffix', width=0, enabled=True)
        self.insert_checkbox_label_line(parent, mode + 'adjustDir', 'Change Directory', mode + 'directory',enabled=True)
        self.insert_checkbox_label_spinbox(parent, mode + 'adjustWedgeAngles', 'Adjust Wedge Angle 1', mode + 'wedge_angle1',
                                           cstep=-1, rstep=1, value=30, stepsize=1, wtype=QSpinBox)
        self.insert_label_spinbox(parent, mode + 'wedge_angle2', 'Wedge Angle 2', cstep=-2, value=30, wtype=QSpinBox)
        self.insert_checkbox_label_spinbox(parent, mode + 'adjustBinning', 'Multiply Pick Positions', mode + 'binning',
                                           value=1, stepsize=1, minimum=0, maximum=32, wtype=QDoubleSpinBox)
        self.insert_checkbox_label_spinbox(parent, mode + 'multiplyShifts', 'Multiply Shifts',
                                           mode + 'factorMultiplyShifts', value=1, stepsize=1, wtype=QDoubleSpinBox,
                                           minimum=0, rstep=0, cstep=4)
        self.insert_label(parent,'', sizepolicy=self.sizePolicyA, cstep=-6,rstep=1)




        self.widgets['particleList0'].textChanged.connect(lambda d, pl='particleList0': self.updatePLs(pl))
        setattr(self, mode + 'gb_changeXMLparameters', groupbox)
        return groupbox

    def extractParticles(self, mode='', title=''):
        tooltip = 'Tick this box to create a particle list including only specific particles.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.insert_label_line_push(parent, 'Particle List', 'particleList1', 'Select your particleList (XML) file.',
                                    initdir=self.projectname, mode='file', cstep=-4, filetype=['xml'],
                                    action2=self.countNumberOfParticles)
        self.insert_checkbox_label_line(parent, mode + 'extractByTomoName', 'By Tomogram Name', mode + 'tomoname',
                                        value='all', width=0, enabled=True)
        self.insert_checkbox_label_line(parent, mode + 'extractByClass', 'By Class', mode + 'classes', cstep=4, rstep=0,
                                        enabled=True, value='')
        self.insert_label(parent, '', sizepolicy=self.sizePolicyA)

        self.widgets['particleList1'].textChanged.connect(lambda d, pl='particleList1': self.updatePLs(pl))
        setattr(self, mode + 'gb_extractParticles', groupbox)
        return groupbox

    def actions(self, mode='', title=''):
        tooltip = 'Tick this box to update specific fields in an particle list using actions such as rotation.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        self.insert_label_line_push(parent, 'Particle List', 'particleList2', 'Select your particleList (XML) file.',
                                    initdir=self.projectname, mode='file', cstep=-4, filetype=['xml'],
                                    action2=self.countNumberOfParticles)
        self.insert_label(parent,'', rstep=1)
        self.insert_checkbox_label(parent, mode + 'moveShiftToPickPos', 'Shift To Pick Position', width=200,
                                   rstep=0, cstep=5)
        self.insert_label(parent, '', sizepolicy=self.sizePolicyA, cstep=-5, rstep=1)

        self.insert_label_spinbox(parent, mode + 'binFactorRecon', 'Binning Factor Reconstruction', value=8, minimum=1)
        self.insert_label_spinbox(parent, mode + 'binFactorSubtomo', 'Binning Factor Subtomogram', value=3, minimum=1, cstep=-2)
        self.insert_label(parent, '', rstep=1)

        self.insert_checkbox_label_spinbox(parent, mode + 'recenterParticles', 'Recenter Subtomogram -- X (px)',
                                           mode + 'cX', value=0, wtype=QSpinBox, rstep=1, cstep=-1, minimum=-100)
        self.insert_label_spinbox(parent, mode + 'cY', 'Y (px)', value=0, wtype=QSpinBox, rstep=1, minimum=-100)
        self.insert_label_spinbox(parent, mode + 'cZ', 'Z (px)', value=0, wtype=QSpinBox, rstep=1, minimum=-100)
        self.insert_label_spinbox(parent, mode + 'sizeSubtomo', 'Size subtomograms (px)', value=128, wtype=QSpinBox,
                                  minimum=1, maximum=1000, stepsize=10, cstep=-3, rstep=1)

        self.insert_label(parent,'', rstep=1, cstep=1)
        self.insert_checkbox_label_spinbox(parent, mode + 'reorientParticles', 'Reorient Subtomogram -- Z (deg)',
                                           mode + 'Z1', value=0, wtype=QDoubleSpinBox, decimals=2, rstep=1,cstep=-1)
        self.insert_label_spinbox(parent, mode + 'X',  'X (deg)', value=0, wtype=QDoubleSpinBox, decimals=2,rstep=1)
        self.insert_label_spinbox(parent, mode + 'Z2', 'Z (deg)', value=0, wtype=QDoubleSpinBox, decimals=2,rstep=1,
                                  cstep=-2)
        self.insert_label(parent,'',rstep=1)
        self.insert_checkbox_label(parent, mode + 'randomizeAngles', 'Randomize Angles', width=200, rstep=1)
        self.insert_checkbox_label(parent, mode + 'mirrorCoordinates', 'Mirror Pick Positions', width=200, rstep=1)

        self.widgets['particleList2'].textChanged.connect(lambda d, pl='particleList2': self.updatePLs(pl))


        setattr(self, mode + 'gb_actions', groupbox)
        return groupbox

    def adjustParticleList(self):
        from pytom.bin.extractClassesFromParticleList import extractClassesFromPL
        from pytom.bin.extractTomoNameFromXML import extractParticleListsByTomoNameFromXML

        mode0, mode1, mode2 = self.modes
        particleList = self.widgets['particleList0'].text()
        if not particleList:
            self.popup_messagebox('Warning', 'No particle list selectected.', 'Please select a particle list.')
            return

        if self.adjust_items[0].isChecked():
            mode = mode0
            suffix = self.widgets[mode + 'suffix'].text()
            dir = self.widgets[mode + 'directory'].text()

            w = [self.widgets[mode+'wedge_angle1'].text(), self.widgets[mode+'wedge_angle2'].text()]
            bin = self.widgets[mode + 'binning'].text()
            fmShifts = self.widgets[mode + 'factorMultiplyShifts'].text()
            outputName = self.requestOutputName(folder=self.pickpartfolder)

            values = [suffix, dir, w, bin, fmShifts]
            print(values)
            for n, adj in enumerate( ('adjustSuffix','adjustDir','adjustWedgeAngles','adjustBinning','multiplyShifts')):
                if not self.widgets[mode + adj].isChecked():
                    values[n] *= 0

            print(values)
            suffix, dir, w, bin, fm = values

            if dir:
                folders = dir.split('/')
                tempfolder = '' if dir.startswith('/') else self.subtomofolder
                for folder in folders:
                    tempfolder = os.path.join(tempfolder, folder)
                    if tempfolder == '':
                        tempfolder = '/'
                        continue
                    if not os.path.exists(tempfolder): os.mkdir(tempfolder)

            if bin == '': bin = 1
            if fm == '' or self.widgets[mode + 'multiplyShifts'].isChecked() == False:
                fm = None
            else:
                fm = float(fm)

            if outputName:
                print('Update {}. Output saved as {}.'.format(os.path.basename(particleList), outputName) )
                print(fm)
                updatePL(particleList, outputName, directory=dir, wedgeangles=w, suffix=suffix, multiplypickpos=float(bin), multiplyshift=fm)
            else:
                self.popup_messagebox('Warning', '', 'No valid output filename provided. No file saved.')

        if self.adjust_items[2].isChecked():
            mode=self.modes[1]
            tomoname = self.widgets[mode + 'tomoname'].text()
            classes  = self.widgets[mode + 'classes'].text()
            if self.widgets[mode + 'extractByTomoName'].isChecked():
                directory = self.requestOutputDirectory(self.pickpartfolder)
                if directory: extractParticleListsByTomoNameFromXML(particleList, directory=directory, query=tomoname)
                else: self.popup_messagebox('Warning', '', 'No valid output filename provided. No file saved.')

            if self.widgets[mode + 'extractByClass'].isChecked() and classes:
                outputName = self.requestOutputName(folder=self.pickpartfolder)
                if outputName: extractClassesFromPL(particleList, classes, outputName)
                else: self.popup_messagebox('Warning', '', 'No valid output filename provided. No file saved.')

        if self.adjust_items[4].isChecked():
            mode = self.modes[2]
            values =[]
            for name in ('mirrorCoordinates', 'randomizeAngles', 'moveShiftToPickPos', 'recenterParticles',
                         'reorientParticles'):
                values.append(self.widgets[mode + name].isChecked())
            mirror, randomize, moveShift, recenter, reorient = values

            values, new_center, rotation, angleListDefaultData = [], [], [], False
            for name in ('cX', 'cY', 'cZ', 'Z1', 'X', 'Z2', 'sizeSubtomo', 'binFactorRecon', 'binFactorSubtomo'):
                values.append(self.widgets[mode + name].value())
            cx, cy, cz, z1, x, z2, sizeSubtomo, binRecon, binSubtomo = map(float, values)

            if recenter:
                new_center = [cx, cy, cz]

            if reorient:
                rotation = [z1, x, z2]

            if randomize:
                angleListDefault = os.path.join(self.pytompath, 'angles/angleLists/angles_18_3040.em')
                vol = read(angleListDefault)
                angleListDefaultData = deepcopy(vol2npy(vol)).T.astype(numpy.float32)

            if recenter or reorient or mirror or randomize or moveShift:
                outputName = self.requestOutputName(folder=self.pickpartfolder)
                if outputName:
                    updatePL(particleList, outputName, new_center=new_center, rotation=rotation, mirror=mirror,
                             anglelist=angleListDefaultData, move_shift=moveShift, tomogram_dir=self.tomogramfolder,
                             binSubtomo=binSubtomo, binRecon=binRecon, sizeSubtomo=sizeSubtomo)
                else:
                    self.popup_messagebox('Warning', 'No valid output name provided', 'No valid output filename provided. No file saved.')


    # Update Methods

    def updatePLs(self, name):
        for key in ('particleList0', 'particleList1', 'particleList2'):
            if not key == name: self.widgets[key].setText(self.widgets[name].text())

    def updateGpuString(self, mode):
        gpu_ids = self.widgets[mode + 'gpuID'].text()
        try:
            a = list(map(int, [el for el in gpu_ids.split(',') if el != '']))
        except:
            self.widgets[mode + 'gpuID'].setText('')
            self.popup_messagebox('Warning', 'Invalid value in field', 'Impossible to parse gpu IDs, field has been cleared.')
            return

        if len(a) > 0:
            self.widgets[mode + 'gpuString'].setText(f'--gpuID {gpu_ids}')
        else:
            self.widgets[mode + 'gpuString'].setText('')

    def updateJobName(self, mode):
        template = os.path.basename(self.widgets[mode+'templateFname'].text()).split('.')[0]
        self.widgets[mode + 'jobName'].setText('job_{}.xml'.format(template))

    def updateTM(self,mode):
        tomoname = os.path.basename(self.widgets[mode + 'tomoFname'].text())
        if not tomoname: return
        filename, file_extension = os.path.splitext(tomoname)
        if not os.path.exists(os.path.join(self.templatematchfolder, 'cross_correlation', filename)):
            os.mkdir(os.path.join(self.templatematchfolder, 'cross_correlation', filename))
        self.execfilenameTM = os.path.join( self.templatematchfolder, 'cross_correlation', filename, 'templateMatch.sh')
        self.xmlfilename = os.path.join(self.templatematchfolder, 'cross_correlation', filename, 'job.xml')
        self.widgets[mode + 'outfolderTM'].setText(os.path.dirname(self.xmlfilename))
        self.widgets[mode + 'jobName'].setText(self.xmlfilename)

        try:
            folder = os.path.dirname(os.popen(f'ls -alrt {self.widgets[mode+"tomoFname"].text()}').read()[:-1].split(' ')[-1])
            widthfile = os.path.join(folder, 'z_limits.txt')
            if os.path.exists(widthfile):
                start, end = map(int, list(numpy.loadtxt(widthfile)))
            else:
                start, end = 0, 0
        except Exception as e:
            print(e)
            start, end = 0, 0

        self.widgets[mode + 'endZ'].setValue(end)
        self.widgets[mode + 'startZ'].setValue(start)
        self.updateZWidth(mode)

        try:
            folder = os.path.join(self.tomogram_folder, tomoname[:len('tomogram_000')], 'sorted')
            print(folder)
            metafiles = [ os.path.join(folder, dir) for dir in os.listdir(folder) if not os.path.isdir(dir) and dir.endswith('.meta')]
            if metafiles:
                from pytom.gui.guiFunctions import datatype, loadstar

                metadata = loadstar(metafiles[0],dtype=datatype)
                tiltAngles = metadata['TiltAngle']
                self.widgets[mode + 'Wedge1'].setValue(int(numpy.round(tiltAngles.min()+90)))
                self.widgets[mode + 'Wedge2'].setValue(int(90-numpy.round(tiltAngles.max())))
        except Exception as e:
            print(e)

    def updateZWidth(self, mode):
        start, end = self.widgets[mode + 'startZ'].value(), self.widgets[mode + 'endZ'].value()
        width = end-start
        if end < start:
            width=0
            self.widgets[mode + 'endZ'].setValue(start)
        self.widgets[mode + 'widthZ'].setText(str(int(round(width))))

        from pytom.basic.files import read
        tomogramFile = self.widgets[mode + 'tomoFname'].text()
        if os.path.exists(tomogramFile):
            v = read(tomogramFile)
            widthX = v.sizeX() if width else 0
            widthY = v.sizeY() if width else 0
        else:
            widthX = widthY = width
        self.widgets[mode + 'widthX'].setText(str(int(round(widthX))))
        self.widgets[mode + 'widthY'].setText(str(int(round(widthY))))

    def createTemplateMask(self, mode):
        title = "Create Mask for Template Matching"
        tooltip = 'Run pytom template matching routine.'
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.ccfolder = os.path.join(self.templatematchfolder, 'cross_correlation')

        self.insert_label_line_push(parent, 'Tomogram', mode + 'tomogramFname', initdir=self.ccfolder,
                                    tooltip='Select the particle list.', mode='file', filetype=['em','mrc'],cstep=-1)
        self.insert_pushbutton(parent, 'Pick!', tooltip='Select the folder with the aligned tilt images.',
                               rstep=1, cstep=2, action=self.cm, params = [parent])

        self.insert_label(parent, cstep=-self.column, sizepolicy=self.sizePolicyA,rstep=1)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyB,width=200)

    def updateMaskGenerate(self,mode):
        if self.widgets[mode + 'maskFile'].text():
            self.widgets[mode+'maskGenerate'].setText(' \\\n    --mask {}'.format(self.widgets[mode + 'maskFile'].text()))

    def update_prefix_pL(self, mode):
        coords, prefixS = mode + 'CoordinateFile', mode + 'PrefixSubtomo'
        prefix = os.path.basename(self.widgets[coords].text())[:-4].replace('coords_', '')
        self.widgets[prefixS].setText(os.path.join('Subtomograms', prefix, 'particle_'))
        self.setFnameParticleList(['file',self.widgets[mode + 'FnameParticleList']], ask=False)

    def update_flag(self,mode):
        w1 = self.widgets[mode+'flagRandomize']
        w2 = self.widgets[mode+'Randomize']
        w1.setText(w2.isChecked()*'-r')

    def update_prefix(self,params):
        prefix = os.path.basename(self.widgets[params[0]].text())[:-4].replace('coords_','')
        self.widgets[params[1]].setText( os.path.join('Subtomograms', prefix) )

    def update_change_tab32_table(self,widget1, widget2):
        widget2.setText('particleList_{}.xml'.format(widget1.text().replace('/particle_','')))


    # Other methods

    def cm(self):
        if self.no_image == False:
            self.partpick = ParticlePicker(self)
            self.partpick.show()

    def jobXMLChanged(self, mode):
        jobXML = self.widgets[mode+'jobXML'].text()
        if not jobXML: return
        folder = os.path.basename(os.path.dirname(jobXML))

        suffix = os.path.basename(jobXML)[3:-4]

        scores = os.path.join( os.path.dirname(jobXML), f'scores{suffix}.em')
        angles = os.path.join( os.path.dirname(jobXML), f'angles{suffix}.em')

        if '_Mirrored' in os.path.basename(jobXML):
            s = '_Mirrored'
        else:
            s = ''
        particleList = os.path.join(self.pickpartfolder, f'particleList_TM_{folder}{s}.xml')
        if not os.path.exists(os.path.dirname(particleList)): os.mkdir(os.path.dirname(particleList))
        self.widgets[mode + 'particleList'].setText(particleList)
        if os.path.exists(scores): self.widgets[mode + 'scoreFile'].setText(scores)
        if os.path.exists(angles): self.widgets[mode + 'anglesFile'].setText(angles)

    def particleListChanged(self, mode):
        particleList = self.widgets[mode + 'particleList'].text()
        if not particleList: return
        folder = os.path.basename(particleList)[:-4]
        self.widgets[mode + 'particlePath'].setText('Subtomograms/{}'.format(folder))

    def create_maskfile(self,params):
        maskfilename = CreateMaskFile(self,maskfname=params[-1])

    def pdb2em(self, params):
        print(params)
        ConvertEM2PDB(self, emfname=params[0],folder=self.widgets[params[1]+'outfolderTM'].text())

    def getMaskFiles(self, key=''):
        try: self.maskFiles.text()
        except: self.maskFiles = QLineEdit()
        self.batchTM.close()
        self.batchTM = SelectFiles(self, initdir=self.ccfolder, search='file', filter=['em','mrc'], id=key,
                                   outputline=self.maskFiles, run_upon_complete=self.populate_batch_templatematch,
                                   title='Select Masks.')

    def createFolder(self,mode):
        prefix=self.widgets[mode+'PrefixSubtomo'].text()
        a = os.path.join(self.subtomofolder, os.path.dirname(prefix))
        if not os.path.exists(a): os.mkdir(a)

    def unsetCheckBoxes(self, mode, num):
        if self.adjust_items[num*2].isChecked():
            for n, m in enumerate(self.modes):
                if m != mode and self.adjust_items[n*2].isChecked():
                    self.toggle_groupbox_visible(self.adjust_items[n*2+1], self.adjust_items[n*2])
                    self.adjust_items[n * 2].setChecked(False)

    def setFnameParticleList(self, params, ask=True):
        t = os.path.basename(self.widgets[self.stage + 'partlist_' + 'CoordinateFile'].text())
        t = t.replace('coords_', '').replace('.txt', '')
        initdir = os.path.join(self.projectname, '04_Particle_Picking/Picked_Particles',
                               'particleList_{}.xml'.format(t))

        if ask:
            particleListFname = str(
                QFileDialog.getSaveFileName(self, 'Save particle list.', initdir, filter='*.xml')[0])
        else:
            particleListFname = initdir

        if particleListFname:
            if particleListFname.endswith('.xml'):
                params[1].setText(particleListFname)
            else:
                params[1].setText('')

    def requestOutputName(self, ext='xml', folder='', title='Save particle list.'):
        fname = str(QFileDialog.getSaveFileName(self, title, folder, filter='*.'+ext)[0])
        if not fname: return
        ext = '.'+ext
        if not fname.endswith(ext): fname += ext
        return fname

    def requestOutputDirectory(self, folder='', title='Output Folder'):
        return QFileDialog.getExistingDirectory(self, title, folder)

    def insert_image(self,params):
        if self.no_image == False:
            #self.partpick = CreateMaskTM(self)
            self.partpick = ParticlePicker(self)
            self.partpick.show()
            #params[0].addWidget(self.partpick,self.row+1,0,4,4)
