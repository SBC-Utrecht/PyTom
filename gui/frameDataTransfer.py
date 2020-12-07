import sys
import os
import random
import glob
import numpy
import time
import atexit

from os.path import dirname, basename
from multiprocessing import Manager, Event, Process
from ftplib import FTP_TLS, FTP

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets

from pytom.gui.fiducialAssignment import FiducialAssignment
from pytom.gui.guiStyleSheets import *
from pytom.gui.guiStructures import *
from pytom.gui.guiFunctions import avail_gpu
from pytom.gui.guiSupportCommands import *
import pytom.gui.guiFunctions as guiFunctions

class BrowseWindowRemote(QMainWindow):
    '''This class creates a new windows for browsing'''
    def __init__(self, parent=None, initdir='/',filter=[''],search='file',credentials=['','',''],outputline=''):
        super(BrowseWindowRemote, self).__init__(parent)
        self.setGeometry(50, 50, 400, 300)
        self.pathdisplay = QLineEdit()
        self.pathdisplay.setEnabled(False)
        splitter0 = QSplitter(Qt.Vertical)
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand
        self.topleft = QListWidget()
        self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)
        self.topright = QListWidget()

        if search == 'file': self.topright.itemDoubleClicked.connect(self.select_file)
        splitter1 = QSplitter(Qt.Horizontal)
        splitter1.addWidget(self.topleft)
        splitter1.addWidget(self.topright)
        splitter1.setSizes([200,200])

        bottom_layout = QHBoxLayout()
        ok_button = QPushButton(text='OK')
        ok_button.clicked.connect(self.select_item)
        cancel_button = QPushButton(text='Cancel')
        cancel_button.clicked.connect(self.cancel)

        bottom = QWidget()
        bottom_layout.addWidget(ok_button)
        bottom_layout.addWidget(cancel_button)
        bottom_layout.addStretch(1)
        bottom.setLayout(bottom_layout)

        splitter0.addWidget(self.pathdisplay)
        splitter0.addWidget(splitter1)
        splitter0.addWidget(bottom)

        self.success = False
        self.ftps=None
        self.initdir = initdir
        self.folderpath = initdir
        self.filters = filter
        self.search = search
        self.outputline = outputline
        self.activeProcesses = {}

        try:
            self.servername = str(credentials[0])
            self.username   = str(credentials[1])
            self.password   = str(credentials[2])
        except:
            pass
            #self.servername,self.username,self.password = 'emsquare1.science.uu.nl','emuser','#99@3584cg'

        try:
            self.connect_ftp_server(self.servername, self.username, self.password)
            self.setCentralWidget(splitter0)
            self.success = True
            self.add_folders()
            self.show()
        except:
            QMessageBox().critical(self, "Credentials not valid.",
                                   "Please provide a valid combination of servername, username, and password",
                                   QMessageBox.Ok)
            #self.destroy(True,True)
            self.close()





    def select_item(self):
        if self.outputline: self.outputline.setText(self.folderpath)
        self.close()

    def select_file(self):
        self.folderpath = os.path.join(self.folderpath, self.topright.currentItem().text())
        self.outputline.setText(self.folderpath)
        self.close()

    def cancel(self):
        #self.pathselect = self.initdir
        self.close()

    def connect_ftp_server(self,servername, username, password):
        self.ftps = FTP_TLS(servername,username,password)
        self.ftps.prot_p()

    def append(self,line):
        words = line.split(None, 3)
        datestamp = words[0]
        timestamp = words[1]
        folder = words[2]
        filename = words[3]
        size = None
        if folder == '<DIR>':
            folder = True
        else:
            size = int(folder)
            folder = False
        if folder:
            self.subdirs.append(filename)
        else:
            self.matchingfiles.append(filename)

    def add_folders(self, folders=[]):
        self.data()

        self.topleft.clear()
        self.topright.clear()

        self.topleft.insertItems(0,self.subdirs)
        self.topright.insertItems(0,self.matchingfiles)

        #if self.search == 'file': self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)

    def data(self):
        self.subdirs = [os.pardir]
        self.matchingfiles = []


        self.ftps.cwd(self.folderpath)
        self.ftps.retrlines('LIST', self.append)

        matchingfiles = []
        for pat in self.filters:
            matchingfiles += [line for line in self.matchingfiles if line.endswith(pat.split('.')[-1])]
        self.matchingfiles = matchingfiles

    def repopulate_folder_list(self):

        extra = self.topleft.currentItem().text()
        if extra != '..':
            self.folderpath = os.path.join(self.folderpath, extra)

        elif len(self.folderpath) >1:
            self.folderpath = os.path.dirname(self.folderpath)
        self.pathdisplay.setText(self.folderpath)
        self.add_folders()

class CollectPreprocess(GuiTabWidget):
    '''Collect Preprocess Widget'''
    def __init__(self, parent=None):
        super(CollectPreprocess, self).__init__(parent)
        self.stage = 'v01_'
        self.projectname = self.parent().projectname
        self.logfolder = self.parent().logfolder
        self.tomogram_folder = self.parent().tomogram_folder
        self.rawnanographs_folder = self.parent().rawnanographs_folder
        self.motioncor_folder = self.parent().motioncor_folder
        self.qtype = self.parent().qtype
        self.qcommand = self.parent().qcommand
        self.widgets = {}
        self.qparams = self.parent().qparams
        self.localqID = {}
        self.activeProcesses = {}
        self.workerID = 0
        self.threadPool = self.parent().threadPool

        self.widgets['pytomPath'] = QLineEdit()
        self.widgets['pytomPath'].setText(self.parent().pytompath)

        self.pytompath = self.parent().pytompath
        self.tabs_dict, self.tab_actions = {}, {}

        headers = ["Data Collection", "Motion Correction"]
        subheaders = [[], []]
        tabUIs = [self.tab1UI,
                  self.tab2UI]
        static_tabs = [[True],[True]]
        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, tabUIs=tabUIs, tabs=self.tabs_dict, tab_actions=self.tab_actions)

        self.widgets = {}
        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.num_nodes = {}
        self.modes = {}
        self.checkbox= {}

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if static_tabs[i][j]:  # tt in ('tab11', 'tab31', 'tab32', 'tab41', 'tab42'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()
                self.tables[tt] = QWidget()
                self.pbs[tt] = QWidget()
                self.ends[tt] = QWidget()
                self.ends[tt].setSizePolicy(self.sizePolicyA)

                if not static_tabs[i][j]:  # tt in ('tab12'):
                    button = QPushButton('Refresh Tab')
                    button.setSizePolicy(self.sizePolicyC)
                    button.clicked.connect(lambda d, k=tt, a=self.tab_actions[tt]: a(k))

                    self.table_layouts[tt].addWidget(button)
                    self.table_layouts[tt].addWidget(self.ends[tt])
                else:  # if not tt in ('tab12'):
                    self.tab_actions[tt](tt)

                tab = self.tabs_dict[tt]
                tab.setLayout(self.table_layouts[tt])


    def tab0UI(self):
        '''Fill out the second tab of the collect & preprocess widget.'''
        gridLayout = QtWidgets.QGridLayout()
        gridLayout.setContentsMargins(10, 10, 10, 10)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = gridLayout
        self.row, self.column = 0, 0

        self.tab2.setLayout(gridLayout)

    def tab1UI(self, key):
        '''Fill out the second tab of the collect & preprocess widget.'''
        print('tabb1UI')
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'DataCollection_', self.stage + 'InsertCredentialsCD_'
        self.modes[t1] = 0
        self.modes[1] = t1

        items += list(self.create_expandable_group(self.createCollectGroup, self.sizePolicyB,
                                                   'Data Collection', mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

        return

    def tab2UI(self, key=''):
        '''Fill out the second tab of the collect & preprocess widget.'''
        print('tabb2UI')
        grid = self.table_layouts[key]
        grid.setAlignment(self, Qt.AlignTop)

        items = []

        t0, t1 = self.stage + 'MotionCorrection', self.stage + 'InsertCredentialsMC_'
        self.modes[t1] = 1
        self.modes[0] = t1

        items += list(self.create_expandable_group(self.createMotioncorGroup, self.sizePolicyB,
                                                   'Motion Correction', mode=t0))
        items[-1].setVisible(False)

        for n, item in enumerate(items):
            grid.addWidget(item, n, 0, 1, 3)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        grid.addWidget(label, n + 1, 0, Qt.AlignRight)

    def tab3UI(self):
        '''Fill out the second tab of the collect & preprocess widget.'''
        gridLayout = QtWidgets.QGridLayout()
        gridLayout.setContentsMargins(10, 10, 10, 10)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = gridLayout
        self.row, self.column = 0, 0

        self.insert_checkbox(parent, 'v01_batch_collect_data', cstep=1, logvar=True)
        self.insert_label(parent, text='Collect Data', columnspan=2, cstep=2)
        self.insert_checkbox(parent, 'v01_batch_onthefly', text='On the fly', cstep=0, rstep=1, logvar=True)
        self.insert_label(parent, rstep=1)
        self.insert_label(parent, text='Search Remotely', rstep=1, cstep=self.column * -1 + 2,
                          tooltip='By checking "Search remotely" one can search on a remote server.\n' +
                                  'For this option, please fill out the server credentials below.')

        # NANOGRAPH ROW
        self.insert_label(parent, text=' - Nanograph Folder', cstep=1,
                          tooltip='Select the folder and the data format in which the nanographs are stored.')
        self.insert_checkbox(parent, 'v01_batch_remote_nanograph', cstep=1, alignment=QtCore.Qt.AlignHCenter,
                             logvar=True)
        self.insert_combobox(parent, 'v01_batch_filetype_nanographs', ['tif', 'mrc', 'em'], cstep=1, logvar=True)
        self.insert_lineedit(parent, 'v01_batch_folder_nanographs', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1 + 1, rstep=1, text='Browse',
                               action=self.browse,
                               params=['folder', self.widgets['v01_batch_folder_nanographs'],
                                       self.widgets['v01_batch_filetype_nanographs'],
                                       self.widgets['v01_batch_remote_nanograph']])

        # MDOC LINE
        self.insert_checkbox(parent, 'v01_batch_collect_mdoc', cstep=1, alignment=QtCore.Qt.AlignHCenter, logvar=True)
        self.insert_label(parent, text=' - Mdoc Folder', cstep=1,
                          tooltip='If you have accompanying mdoc files, select the folder where these files are stored')
        self.insert_checkbox(parent, 'v01_batch_remote_mdoc', cstep=1, logvar=True)
        self.insert_combobox(parent, 'batch_filtype_mdoc', ['mdoc'], cstep=1, logvar=True)
        self.insert_lineedit(parent, 'v01_batch_folder_mdoc', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1, rstep=1, text='Browse',
                               action=self.browse,
                               params=['folder', self.items[self.row][self.column - 1],
                                       self.items[self.row][self.column - 2], self.items[self.row][self.column - 3]])

        # EMPTY LINE
        self.insert_label(parent, rstep=1)

        # Motion Correction
        self.insert_checkbox(parent, 'v01_batch_motioncorrection', cstep=1, logvar=True)
        self.insert_label(parent, text='Do Motion Correction', columnspan=2, rstep=1,
                          tooltip='Execute Motion correction xsusing motioncor2')

        # Gain file
        self.gain_correct = self.insert_checkbox(parent, '01_batch_gaincorrection', cstep=1, logvar=True)
        self.insert_label(parent, text=' - Gain File', cstep=1,
                          tooltip='Apply gain correction using the file specified in entry line.')
        self.insert_checkbox(parent, 'v01_batch_remote_gainfile', cstep=1, logvar=True)
        self.insert_combobox(parent, 'v01_batch_filetype_gaincorrection', ['dm4', 'mrc'], cstep=1, logvar=True)
        self.insert_lineedit(parent, 'v01_batch_gainfile', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1 + 1, rstep=1, text='Browse', action=self.browse,
                               params=['file', self.v01_batch_gainfile,
                                       self.v01_batch_remote_gainfile, self.v01_batch_remote_gainfile])

        # Patch Size
        self.insert_checkbox(parent, 'v01_batch_use_patch', cstep=1,logvar=True)
        self.insert_label(parent, text=' - Patch Size', cstep=3,
                          tooltip='Sets the size of the patches used for motioncorrection (NxN).')
        self.insert_lineedit(parent, 'v01_batch_patchsize', cstep=self.column * -1 + 1, rstep=1,
                             value='5', validator=QIntValidator(), logvar=True)

        # FtBin Size
        self.insert_checkbox(parent, 'v01_batch_use_ftbin', cstep=1,logvar=True)
        self.insert_label(parent, text=' - FtBin Factor', cstep=3,
                          tooltip='Sets the binning factor for Fourier binning of your tilt images.')
        self.insert_lineedit(parent, 'v01_batch_ftbin_factor', cstep=self.column * -1 + 1, rstep=1,
                             value='2', validator=QIntValidator(), logvar=True)

        # EMPTY LINE
        self.insert_label(parent, rstep=1)
        self.insert_label(parent, text='Credentials Remote Login', rstep=1, cstep=1, columnspan=2)

        # SERVERNAME, USERNAME, PASSWORD
        self.insert_label(parent, text=' - Server Name', cstep=3, tooltip='Sets the name of the remote server.')
        self.insert_lineedit(parent, 'servername1', cstep=self.column * -1 + 2, rstep=1, att=True)
        self.insert_label(parent, text=' - User Name', cstep=3,
                          tooltip='Sets the username used for login into the remote server.')
        self.insert_lineedit(parent, 'username1', cstep=self.column * -1 + 2, rstep=1, att=True)
        self.insert_label(parent, text=' - Password', cstep=3,
                          tooltip='Sets the password used for login into the remote server.')
        self.insert_lineedit(parent, 'password1', rstep=2, password=True, att=True)
        self.insert_pushbutton(parent, text='Run', action=self.download, params=[5],wname='CandP')


        self.widgets['username1'].setText('emuser')
        self.widgets['servername1'].setText('emsquare1.science.uu.nl')
        self.statusBar = self.parent().statusBar

        self.tab2.setLayout(gridLayout)

    def createCollectGroup(self,mode='v01_DataCollection_'):
        print('InsertDataCollectionGB', mode)
        title = "Data Collection"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, rstep=1, cstep=-1)#, text='Search Remotely',rstep=1,cstep=-1,
                          #tooltip='By checking "Search remotely" one can search on a remote server.\n' +
                          #        'For this option, please fill out the server credentials below.')

        # Nanpgraph row
        self.insert_label(parent, text=' - Nanograph Folder', cstep=1, width=160,
                          tooltip='Select the folder and the data format in which the nanographs are stored.')
        self.widgets[mode+'remote_nanograph'] = QCheckBox()
        self.insert_combobox(parent, mode+'filetype_nanographs', ['tif', 'mrc', 'st'], cstep=1, logvar=True)
        self.insert_lineedit(parent, mode+'folder_nanographs', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1 , rstep=1, text='Browse',
                               action=self.browse, width=100,
                               params=['folder', self.items[self.row][self.column - 1],
                                                 self.items[self.row][self.column - 2],
                                                 self.widgets[mode + 'remote_nanograph']])


        # MDOC LINE
        #self.insert_checkbox(parent, mode +'collect_mdoc', cstep=1, alignment=QtCore.Qt.AlignHCenter,
        #                     logvar=True, width=20)
        self.insert_label(parent, text=' - Mdoc Folder', cstep=1, width=120,
                      tooltip='If you hae accompanying mdoc files, select the folder where these files are stored')
        self.widgets[mode + 'remote_mdoc'] = QCheckBox()
        self.insert_combobox(parent, mode + 'filtype_mdoc', ['mdoc'], cstep=1)
        self.insert_lineedit(parent, mode + 'folder_mdoc', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=0, rstep=1, text='Browse', width=100,
                               action=self.browse,
                               params=['folder', self.items[self.row][self.column - 1],
                                                 self.items[self.row][self.column - 2],
                                                 self.widgets[mode + 'remote_mdoc']])
        self.insert_label(parent,rstep=1)
        self.insert_pushbutton(parent, text='Run', rstep=1, action=self.download, params=mode, width=100)

        setattr(self,mode+'gb_collect', groupbox)
        return groupbox

    def createMotioncorGroup(self,mode='v01_single_'):

        title = "Motion Correction"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 170

        self.insert_label(parent, rstep=1, cstep=-1)

        # Nanograph row
        self.insert_label(parent, text=' - Nanograph Folder', cstep=1, width=160,
                          tooltip='Select the folder and the data format in which the nanographs are stored.')
        self.widgets[ mode + 'remote_nanograph'] = QCheckBox()
        self.insert_combobox(parent, mode + 'filetype_nanographs', ['tif', 'mrc'], cstep=1, logvar=True)
        self.insert_lineedit(parent, mode + 'folder_nanographs', cstep=3, logvar=True, columnspan=3,
                             value=self.rawnanographs_folder)
        self.insert_pushbutton(parent, cstep=self.column * -1, rstep=1, text='Browse',
                               action=self.browse, width=100,
                               params=['folder', self.items[self.row][self.column - 3],
                                       self.items[self.row][self.column - 4],
                                       self.widgets[mode + 'remote_nanograph']])

        #self.gain_correct = self.insert_checkbox(parent, mode+'gaincorrection', cstep=1, logvar=True, width=20)
        self.insert_label(parent, text=' - Gain File', cstep=1, alignment=QtCore.Qt.AlignLeft,width=120,
                          tooltip='Apply gain correction using the file specified in entry line.')
        self.widgets[mode+'remote_gainfile'] = QCheckBox()
        self.insert_combobox(parent, mode+'filetype_gaincorrection', ['dm4', 'mrc'], cstep=1)
        self.insert_lineedit(parent, mode+'gainFile', cstep=3, logvar=True, columnspan=3)
        self.insert_pushbutton(parent, cstep=self.column * -1 , rstep=1, text='Browse', action=self.browse,
                               params=['file',  self.items[self.row][self.column - 3],
                                                self.items[self.row][self.column - 4],
                                                self.widgets[mode + 'remote_gainfile']])

        # Gpu ID's
        self.insert_label(parent, text=" - Gpu IDs", cstep=2, alignment=QtCore.Qt.AlignLeft, width=120,
                          tooltip='Gpu ID"s of the gpus you want to use. Multiple gpu IDs need to be separated by a space.')
        self.insert_lineedit(parent, mode + 'gpuId', cstep=self.column * -1, rstep=1, columnspan=3, value='0',
                             logvar=True)

        # Patch Size
        self.insert_label(parent, text=' - Patch Size', cstep=2, alignment=QtCore.Qt.AlignLeft, width=120,
                          tooltip='Sets the size of the patches used for motioncorrection (NxN).')
        self.insert_lineedit(parent, mode + 'patchSize', cstep=self.column * -1 , rstep=1, columnspan=3,
                             validator=QIntValidator(), logvar=True)

        # FtBin Size
        self.insert_label(parent, text=' - FtBin Factor', cstep=2, alignment=QtCore.Qt.AlignLeft, width=120,
                          tooltip='Sets the binning factor for Fourier binning of your tilt images.')
        self.insert_lineedit(parent, mode+'ftBin', cstep=2 , rstep=1, columnspan=3,
                             validator=QIntValidator(), logvar=True)

        self.widgets[mode + 'fileTypeCapitalized'] = QLineEdit('Tiff')
        self.widgets[mode + 'gainFileFlag'] = QLineEdit('')
        self.widgets[mode + 'patchSizeFlag'] = QLineEdit('')
        self.widgets[mode + 'ftBinFlag'] = QLineEdit('')
        self.widgets[mode + 'gpuIdFlag'] = QLineEdit('')

        for flag, name in (('-Gain ', 'gainFileFlag'), ('-Patch ', 'patchSizeFlag'), ('-FtBin ', 'ftBinFlag'), ('-Gpu ', 'gpuIdFlag')):
            self.widgets[mode + name[:-4]].textChanged.connect(lambda d, n=mode+name, f=flag: self.updateFlag(n, f))
            self.updateFlag(mode+name, flag)

        execfilename = os.path.join(self.motioncor_folder, 'MotionCorrection.sh')
        paramsSbatch = guiFunctions.createGenericDict(fname='motionCorrection', folder=self.logfolder,
                                                      id='MotionCorrection')
        paramsCmd = [self.motioncor_folder, mode + 'fileTypeCapitalized', mode + 'folder_nanographs',
                     self.motioncor_folder, mode + 'gainFileFlag', mode + 'patchSizeFlag', mode + 'ftBinFlag',
                     mode + 'gpuIdFlag', templateMotionCorrection]

        self.insert_gen_text_exe(parent, mode, paramsCmd=paramsCmd, exefilename=execfilename, paramsSbatch=paramsSbatch,
                                 cs=3, queue=True)

        setattr(self, mode + 'gb_mcor', groupbox)
        return groupbox

    def createCredentialsGroup(self,mode=''):
        title = "Credentials Remote Login"
        tooltip = ''
        sizepol = self.sizePolicyB
        groupbox, parent = self.create_groupbox(title, tooltip, sizepol)

        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        w = 250

        # SERVERNAME, USERNAME, PASSWORD
        self.insert_label(parent, text=' - Server Name', cstep=1, tooltip='Sets the name of the remote server.', width=160)
        self.insert_lineedit(parent, mode + 'servername', cstep=1 , rstep=1,width=w)
        self.insert_label(parent, cstep=-2,columnspan=2)

        self.insert_label(parent, text=' - User Name', cstep=1,
                          tooltip='Sets the username used for login into the remote server.')
        self.insert_lineedit(parent, mode + 'username', cstep=1, width=w)
        self.insert_label(parent, cstep=-2, rstep=1, width=20)

        self.insert_label(parent, text=' - Password', cstep=1,
                          tooltip='Sets the password used for login into the remote server.')
        self.insert_lineedit(parent, mode + 'password', rstep=2, password=True, width=w)

        for name in ('servername', 'username', 'password'):
            self.widgets[mode + name].textChanged.connect(lambda d, m=mode: self.updateCredentials(m))

        setattr(self, mode + 'gb_cred', groupbox)
        return groupbox

    def download(self, mode):
        if not self.widgets[mode + 'folder_nanographs'].text():
            self.popup_messagebox('Info', 'No Folder Selected', 'Please select a nanograph folder')
            return

        step = 1
        widget = QWidget(self)
        layout = QHBoxLayout()
        widget.setLayout(layout)

        self.progressBar_Collect = QProgressBar()
        self.progressBar_Collect.setSizePolicy(self.sizePolicyA)
        self.progressBar_Collect.setStyleSheet(DEFAULT_STYLE_PROGRESSBAR)
        self.progressBar_Collect.setFormat('Collect Data: %v/%m')
        self.progressBar_Collect.setSizePolicy(self.sizePolicyB)
        layout.addWidget(self.progressBar_Collect)
        self.statusBar.addPermanentWidget(self.progressBar_Collect,stretch=1)

        self.em_fscollect(mode)

    def retrieve_files(self, ftps, flist, folder_local, value=0, only_missing=True):
        existing_files = glob.glob('%s/*' % folder_local)
        total = 0
        num_files_downloaded = 0
        for fname in flist:
            if not os.path.exists(fname):
                outname = "%s/%s" % (folder_local, fname.split('/')[-1])
                if not outname in existing_files:
                    gFile = open("%s.temp" % outname, "wb")
                    ftps.retrbinary('RETR %s' % fname, gFile.write)
                    gFile.close()
                    if '.tif' == outname[-4:]:
                        os.system('mv %s.temp %s' % (outname, outname))
                    else:
                        os.system('mv {}.temp {}'.format(outname, outname))
                    num_files_downloaded += 1
            else:
                outname = "{}/{}".format(folder_local, os.path.basename(fname) )
                if os.path.exists(fname) and not os.path.exists(outname):

                    if '.tif' == outname[-4:]:
                        try:
                            os.link(fname, outname)
                        except:
                            os.system('cp -f {} {}'.format(fname, outname))

                    elif outname[-3:] == '.st':


                        for end in ('.tlt', '.prexg', '.fid'):
                            outname_temp = outname.replace('.st', end)
                            fname_temp = fname.replace('.st', end)

                            if os.path.exists(fname_temp):
                                try:
                                    os.link(fname_temp, outname_temp)
                                except:
                                    os.system('cp -f {} {}'.format(fname_temp, outname_temp))

                        try:
                            os.link(fname, outname)
                        except:
                            os.system('cp -f {} {}'.format(fname, outname))


                    else:
                        try:
                            os.link(fname, outname)
                        except:
                            os.system('cp -f {} {}'.format(fname, outname))
                    num_files_downloaded += 1

        return num_files_downloaded

    def meta_retrieve(self, ftps, flist, mdoc_list, folder_remote, folder_local, e, value, only_missing=True):
        import time

        num_flist = self.retrieve_files(ftps, flist, folder_local, value, only_missing=True)
        num_mdoc = self.retrieve_files(ftps, mdoc_list, folder_local, value, only_missing=True)

        e.set()
        if ftps: ftps.close()

    def connect_ftp_server(self, username, servername, password):
        ftps = FTP_TLS(servername, username, password)
        ftps.prot_p()
        return ftps

    def update_progress_collect(self,total):
        self.progressBar_Collect.setValue(total)

    def update_progress_motioncor(self,total):
        self.progressBar_Motioncor.setValue(total)

    def delete_progressbar_motioncor(self):
        self.statusBar.removeWidget(self.progressBar_Motioncor)
        self.popup_messagebox("Info", "Completion", 'Successfully finished motion correction')
        self.widgets['CandP'].setEnabled(True)

    def delete_progressbar_collect(self):
        self.statusBar.removeWidget(self.progressBar_Collect)

        mdocs = glob.glob("{}/*.mdoc.temp".format(self.local_data_folder))
        for mdoc in mdocs:
            os.system("mv {} {}".format(mdoc, mdoc[:-5]))

        guiFunctions.createMetaDataFiles(self.running_folder, '')
        self.popup_messagebox("Info", "Completion", 'Successfully finished data transfer')
        #self.widgets['CandP'].setEnabled(True)

    def em_fscollect(self,mode='v01_batch_'):
        self.number_collection_subprocesses = 1

        if not self.widgets[mode + 'folder_nanographs'].text():
            self.popup_messagebox('Info', 'No Folder Selected', 'Please select a nanograph folder')
            return

        self.remote_data_folder = self.widgets[mode + 'folder_nanographs'].text()
        self.search_remote_nanograph = self.widgets[mode+'remote_nanograph'].isChecked()
        self.data_file_type = self.widgets[mode+'filetype_nanographs'].currentText()
        self.local_data_folder = self.rawnanographs_folder
        self.remote_mdoc_folder = self.widgets[mode+'folder_mdoc'].text()
        collect_mdoc = True if self.remote_mdoc_folder else False

        mdoc_list = []
        flist = []

        procs = []
        flist = []

        e = Event()
        finish_collect = Event()

        try:
            ftps = self.connect_ftp_server(self.servername, self.username, self.password)
        except:
            ftps = 0


        # Create  User Input.

        # Check list of nanograph files (flist), either from remote or local files.
        if self.search_remote_nanograph:
            try:
                if not os.path.splitext(self.remote_data_folder)[-1]:
                    ftps = self.connect_ftp_server(self.servername, self.username, self.password)
                    flist_all = ftps.nlst(self.remote_data_folder)
                else:
                    flist_all = [self.remote_data_folder]
                flist = [tiffile for tiffile in flist_all if tiffile.endswith(self.data_file_type)]
            except:
                self.popup_messagebox("Error", "Error", "Remote data folder does not exist.")
                return

        else:
            try:
                if not os.path.isfile(self.remote_data_folder):
                    flist_all = os.listdir(self.remote_data_folder)
                else:
                    flist_all = [self.remote_data_folder]
                flist = [os.path.join(self.remote_data_folder, tiffile) for tiffile in flist_all if
                         tiffile.endswith(self.data_file_type)]

            except:
                self.popup_messagebox("Error", "Error", "Data folder does not exist.")
                return

        if len(flist) == 0:
            self.popup_messagebox("Error", "Error", "Data folder is empty.")
            if ftps: ftps.close()
            return

        # Check mdoc input.
        if collect_mdoc:
            try:
                if not os.path.splitext(self.remote_mdoc_folder)[-1]:
                    all_files = ftps.nlst(self.remote_mdoc_folder)
                else:
                    all_files = [self.remote_mdoc_folder]
                mdoc_list = [mdocfile for mdocfile in all_files if mdocfile[-4:] == 'mdoc']
            except:
                try:
                    if not os.path.isfile(self.remote_mdoc_folder):
                        all_files = os.listdir(self.remote_mdoc_folder)
                    else:
                        all_files = [self.remote_mdoc_folder]

                    mdoc_list = [os.path.join(self.remote_mdoc_folder, mdocfile) for mdocfile in all_files if
                                 mdocfile[-4:] == 'mdoc']
                except:
                    self.popup_messagebox("Error", "Error", "*.mdoc folder does not exist.")
                    if ftps: ftps.close()
                    return

            if len(mdoc_list) == 0:
                self.popup_messagebox("Error", "Error", "No *.mdoc files in mdoc folder.")
                if ftps: ftps.close()
                return
        self.progressBar_Collect.setMaximum(len(mdoc_list+flist))


        for ni in range(10000):
            ff = f'{self.rawnanographs_folder}/import_{ni:05d}'
            if not os.path.exists(ff):
                os.mkdir(ff)
                break
        self.running_folder = ff
        # Download or copy selected files.
        for i in range(self.number_collection_subprocesses):
            proc = Process(target=self.meta_retrieve,
                           args=(ftps, flist, mdoc_list[i::self.number_collection_subprocesses],
                                 self.remote_data_folder[i::self.number_collection_subprocesses],
                                 self.running_folder, e, 1))
            procs.append(proc)
            proc.start()
            atexit.register(guiFunctions.kill_proc, proc)
        self.flist = flist
        proc = Worker( fn=self.check_run, args=([e], self.flist, self.running_folder) )

        proc.signals.result1.connect(self.update_progress_collect)
        proc.signals.finished_collect.connect(self.delete_progressbar_collect)

        proc.start()

    def em_fscollect2(self, mode='v01_batch_'):

        self.c = self.collect
        self.m = self.motioncor
        self.collection_finished = not self.collect
        self.number_collection_subprocesses = 1
        self.number_motioncor_subprocesses = 4
        self.continuous = False  # self.on_the_fly_data_retrieval.get()

        self.search_remote_nanograph = self.widgets[mode + 'remote_nanograph'].isChecked()
        self.data_file_type = self.widgets[mode + 'filetype_nanographs'].currentText()
        self.remote_data_folder = self.widgets[mode + 'folder_nanographs'].text()
        self.local_data_folder = self.rawnanographs_folder
        self.remote_mdoc_folder = self.widgets[mode + 'folder_mdoc'].text()

        self.search_remote_dm4 = self.widgets[mode + 'remote_gainfile'].isChecked()
        self.dm4_file = self.widgets[mode + 'gainfile'].text()
        self.apply_gaincorrection = self.widgets[mode + 'gaincorrection'].isChecked()
        self.ftbinfactor = int(self.widgets[mode + 'ftbin_factor'].text())
        self.apply_ftbinning = self.widgets[mode + 'use_ftbin'].isChecked()
        self.use_patches = self.widgets[mode + 'use_patch'].isChecked()
        self.patchsize = int(self.widgets[mode + 'patchsize'].text())

        mdoc_list = []
        flist = []

        procs = []
        flist = []
        motioncor_flist = []
        e = Event()
        finish_collect = Event()
        finish_motioncor = Event()
        try:
            ftps = self.connect_ftp_server(self.servername, self.username, self.password)
        except:
            ftps = 0

        # Create  User Input.
        if self.c:
            # Check list of nanograph files (flist), either from remote or local files.
            if self.search_remote_nanograph:
                try:
                    if not os.path.splitext(self.remote_data_folder)[-1]:
                        ftps = self.connect_ftp_server(self.servername, self.username, self.password)
                        flist_all = ftps.nlst(self.remote_data_folder)
                    else:
                        flist_all = [self.remote_data_folder]
                    flist = [tiffile for tiffile in flist_all if tiffile[-3:] == self.data_file_type[-3:]]
                except:
                    self.popup_messagebox("Error", "Error", "Remote data folder does not exist.")
                    return

            else:
                try:
                    if not os.path.isfile(self.remote_data_folder):
                        flist_all = os.listdir(self.remote_data_folder)
                    else:
                        flist_all = [self.remote_data_folder]
                    flist = [os.path.join(self.remote_data_folder, tiffile) for tiffile in flist_all if
                             tiffile[-3:] == self.data_file_type[-3:]]
                except:
                    self.popup_messagebox("Error", "Error", "Data folder does not exist.")
                    return

            if len(flist) == 0:
                self.popup_messagebox("Error", "Error", "Data folder is empty.")
                if ftps: ftps.close()
                return

            # Check mdoc input.
            if self.widgets[mode + 'collect_mdoc'].isChecked():
                try:
                    if not os.path.splitext(self.remote_mdoc_folder)[-1]:
                        all_files = ftps.nlst(self.remote_mdoc_folder)
                    else:
                        all_files = [self.remote_mdoc_folder]
                    mdoc_list = [mdocfile for mdocfile in all_files if mdocfile[-4:] == 'mdoc']
                except:
                    try:
                        if not os.path.isfile(self.remote_mdoc_folder):
                            all_files = os.listdir(self.remote_mdoc_folder)
                        else:
                            all_files = [self.remote_mdoc_folder]

                        mdoc_list = [os.path.join(self.remote_mdoc_folder, mdocfile) for mdocfile in all_files if
                                     mdocfile[-4:] == 'mdoc']
                    except:
                        self.popup_messagebox("Error", "Error", "*.mdoc folder does not exist.")
                        if ftps: ftps.close()
                        return

                if len(mdoc_list) == 0:
                    self.popup_messagebox("Error", "Error", "No *.mdoc files in mdoc folder.")
                    if ftps: ftps.close()
                    return

        # Check user input for gain correction
        if self.m and self.apply_gaincorrection:
            if len(self.dm4_file) < 5:
                self.popup_messagebox("Error", "Error", "Invalid *.dm4 file.")
                return

            if self.search_remote_dm4:
                try:
                    self.retrieve_files(ftps, [self.dm4_file], self.motioncor_folder)
                    self.dm4_file = os.path.join(self.motioncor_folder, os.path.basename(self.dm4_file))
                    if not os.path.exists(self.dm4_file):
                        self.popup_messagebox("Error", "Error", "{} does not exist.".format(self.dm4_file))
                        ftps.close()
                        return
                except:
                    self.popup_messagebox("Error", "Error", "*.dm4 file does not exist.")
                    ftps.close()
                    return
            else:
                if not os.path.exists(self.dm4_file):
                    self.popup_messagebox("Error", "Error", "*.dm4 file does not exist.")
                    return

        # Download or copy selected files.
        if self.c or self.m:

            for i in range(self.number_collection_subprocesses):
                proc = Process(target=self.meta_retrieve,
                               args=(ftps, flist, mdoc_list[i::self.number_collection_subprocesses],
                                     self.remote_data_folder[i::self.number_collection_subprocesses],
                                     self.local_data_folder, e, 1))
                procs.append(proc)
                proc.start()
                atexit.register(guiFunctions.kill_proc, proc)
            self.flist = flist + mdoc_list
            # self.progressBar_Collect.setRange(0, len(self.flist))

        self.motioncor_flist = []

        self.mcor_procs = []

        if self.m:
            # Sleep 100 ms in order to ensure processes are started.
            time.sleep(0.1)
            if not self.c:
                e.set()
                motioncor_flist = [line.split('/')[-1] for line in
                                   glob.glob("{}/*.{}".format(self.local_data_folder, self.data_file_type[-3:]))]

            else:
                motioncor_flist = [line.split('/')[-1] for line in flist]

            if self.apply_gaincorrection and (len(self.dm4_file) < 5 or not os.path.exists(self.dm4_file)):
                print('Error', 'Error', 'no or invalid *.dm4 file submitted')
                return

            if len(motioncor_flist) == 0:
                print('Error', 'Error', 'no *.tif files exist in %s' % self.local_data_folder)
                return

            self.progressBar_Motioncor.setRange(0, len(motioncor_flist))
            self.motioncor_flist = motioncor_flist

            availableGPU = avail_gpu()
            if len(availableGPU) == 0:
                self.popup_messagebox('Warning', 'No Motioncor', 'No available GPUs, motioncor is not running.')
            else:
                num_jobs = min(len(availableGPU), self.number_motioncor_subprocesses)
                for i in range(num_jobs):
                    finish_motioncor = Event()
                    self.mcor_procs.append(finish_motioncor)
                    proc = Process(target=self.apply_motioncor, args=(motioncor_flist[i::num_jobs],
                                                                      i + 1, e, availableGPU[i], self.completedMC,
                                                                      finish_motioncor))
                    procs.append(proc)
                    proc.start()
                    atexit.register(guiFunctions.kill_proc, proc)

        if self.c or self.m:
            self.widgets['CandP'].setEnabled(False)
            proc = Worker(fn=self.check_run, args=(procs, [e] + self.mcor_procs, self.c, self.m, self.flist,
                                                   self.motioncor_flist, self.continuous, self.rawnanographs_folder,
                                                   self.motioncor_folder))

            proc.signals.result1.connect(self.update_progress_collect)
            proc.signals.result2.connect(self.update_progress_motioncor)
            proc.signals.finished_collect.connect(self.delete_progressbar_collect)
            proc.signals.finished_mcor.connect(self.delete_progressbar_motioncor)

            proc.start()

    def check_run(self, events, flist, folder1, signals):


        while len(events):
            events = [proc for proc in events if not proc.is_set()]

            a = [os.path.basename(line) for line in glob.glob( os.path.join(folder1, '*')) ]

            total = 0

            # Check if filename f is in the downloaded folder
            for f in flist:
                if os.path.basename(f) in a: total += 1
            signals.result1.emit(total)

            if total >= len(flist):
                signals.finished_collect.emit()
                break
            time.sleep(0.1)

    def apply_motioncor(self, flist, process_id, e, gpu, value,finish_mcor):
        mc_dir = '{}/02_Preprocessed_Nanographs/Motion_corrected'.format(os.path.dirname(self.local_data_folder))

        if self.apply_gaincorrection:
            gainfile = os.path.join(self.local_data_folder, self.dm4_file) + '.mrc'

        quit = False
        while True:

            for fname in flist:

                logfile = os.path.join(mc_dir, '{}.log'.format(fname[:-4]))
                name, ext = os.path.splitext(os.path.basename(fname))
                if not os.path.exists(os.path.join(self.local_data_folder, os.path.basename(fname))):
                    continue
                if not os.path.exists(os.path.join(mc_dir, name + '.mrc')):
                    intiff = os.path.join(self.local_data_folder, os.path.basename(fname))
                    outmrc = os.path.join(mc_dir, name + '.mrc')
                    gain, patch, ftbin = '', '', ''

                    if self.apply_gaincorrection:
                        gain = '-Gain {}'.format(gainfile)
                    if self.apply_ftbinning:
                        ftbin = '-FtBin {}'.format(self.ftbinfactor, self.ftbinfactor)
                    if self.use_patches:
                        patch = '-Patch {} {}'.format(self.patchsize, self.patchsize)

                    logfile2 = "{}/logfile.motioncor.process{}.txt".format(mc_dir, process_id)

                    mcor2cmd = "-InTiff {} -OutMrc {} {} {} {} -Gpu {} -LogFile {} >> {}".format( intiff, outmrc, gain,
                                                                                                  patch, ftbin, gpu,
                                                                                                  logfile, logfile2)
                    if fname[-3:] == 'mrc':
                        mcor2cmd = mcor2cmd.replace('-InTiff','-InMrc')

                    os.system('motioncor2 ' + mcor2cmd)

                continue

            if e.is_set() and not quit:
                quit = True
                continue
            elif e.is_set() and quit:
                break
            else:
                time.sleep(1)

        finish_mcor.set()

    def updateFlag(self, name, flag):

        text = self.widgets[name[:-4]].text()

        if text:
            if flag == '-Patch ': text = f'{text} {text}'
            if flag == '-Gpu ': text = text.replace(',', ' ')
            self.widgets[name].setText(flag + text)
        else:
            self.widgets[name].setText('')

    def updateFileType(self, mode):
        text = self.widgets[mode + 'filetype_nanograph'].currentText()
        if text == 'tif': text = 'tiff'
        self.widget[mode + 'fileTypeCapitalized'].setText(text.capitalize())

    def updateCredentials(self, mode):
        print(self.modes.keys())
        for name in ('servername', 'username', 'password'):
            setattr(self, name, self.widgets[mode + name].text())
            self.widgets[self.modes[self.modes[mode]] + name].setText(self.widgets[mode + name].text())

def main():
    app = QApplication(sys.argv)
    ex = CollectPreprocess()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
