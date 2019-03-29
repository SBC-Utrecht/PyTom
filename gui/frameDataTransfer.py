import sys
import os
import random
import glob
import numpy
import time

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
        headers = ["Individual Assessment", "Batch"]
        offx, offy, dimx, dimy=0, 0, 900, 721
        self.tomogram_folder = self.parent().tomogram_folder
        self.rawnanographs_folder = self.parent().rawnanographs_folder
        self.motioncor_folder = self.parent().motioncor_folder
        self.addTabs(headers)
        #self.tab1.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        #self.tab2.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        self.tab1UI()
        self.tab2UI()


    def tab0UI(self):
        '''Fill out the second tab of the collect & preprocess widget.'''
        gridLayout = QtWidgets.QGridLayout()
        gridLayout.setContentsMargins(10, 10, 10, 10)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = gridLayout
        self.row, self.column = 0, 0

        self.tab2.setLayout(gridLayout)

    def tab1UI(self):
        '''Fill out the second tab of the collect & preprocess widget.'''

        self.grid_single = QtWidgets.QGridLayout()
        self.grid_single.setAlignment(self, Qt.AlignTop)

        self.csa, self.gsa = self.create_expandable_group(self.createCollectGroup, self.sizePolicyB,
                                                          'Collect Data', mode='v01_single_')
        self.csc, self.gsc = self.create_expandable_group(self.createCredentialsGroup,
                                                          self.sizePolicyB,'Credentials Remote Login')
        self.csm, self.gsm = self.create_expandable_group(self.createMotioncorGroup, self.sizePolicyB,
                                                          'Do Motion Correction', mode='v01_single_')

        self.gsa.setChecked(True)
        self.grid_single.addWidget(self.gsa, 0, 0)
        self.grid_single.addWidget(self.csm, 2, 0)
        self.grid_single.addWidget(self.gsm, 3, 0)
        self.grid_single.addWidget(self.csc, 4, 0)
        self.grid_single.addWidget(self.gsc, 5, 0)
        self.tab1.setLayout(self.grid_single)
        self.tab1.setSizePolicy(self.sizePolicyB)

        p = QPushButton(text='Run')
        p.clicked.connect(lambda ignore, p=[self.gsa,self.csm]: self.download(p))
        p.setFixedWidth(100)
        self.widgets['CandP'] = p
        self.grid_single.addWidget(p,6,0)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.grid_single.addWidget(label, 7, 0, Qt.AlignRight)

        #self.grid_single.addWidget(PlotCanvas(self, width=15, height=4), 8, 0)

    def tab2UI(self):
        '''Fill out the second tab of the collect & preprocess widget.'''
        self.grid_batch = QtWidgets.QGridLayout()
        self.tab2.setLayout(self.grid_batch)

        self.cba, self.gba = self.create_expandable_group(self.createCollectGroup, self.sizePolicyB,
                                                          'Collect Data', mode='v01_batch_')
        self.cbc, self.gbc = self.create_expandable_group(self.createCredentialsGroup, self.sizePolicyB,
                                                          'Credentials Remote Login',mode='1' )
        self.cbm, self.gbm = self.create_expandable_group(self.createMotioncorGroup, self.sizePolicyB,
                                                          'Do Motion Correction', mode='v01_batch_')

        self.grid_batch.addWidget(self.cba, 0, 0)
        self.grid_batch.addWidget(self.gba, 1, 0)
        self.grid_batch.addWidget(self.cbm, 2, 0)
        self.grid_batch.addWidget(self.gbm, 3, 0)
        self.grid_batch.addWidget(self.cbc, 4, 0)
        self.grid_batch.addWidget(self.gbc, 5, 0)

        p = QPushButton(text='Run')
        p.clicked.connect(lambda ignore, p=[self.cba,self.cbm]: self.download(p))
        p.setFixedWidth(100)
        self.grid_batch.addWidget(p, 6, 0, Qt.AlignRight)

        label = QLabel()
        label.setSizePolicy(self.sizePolicyA)
        self.grid_batch.addWidget(label, 7, 0)

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

    def createCollectGroup(self,mode='v01_single_'):
        groupbox = QGroupBox("Collect Data")
        groupbox.setToolTip('Collect file from local or remote file system')
        groupbox.setCheckable(False)
        groupbox.setChecked(True)
        groupbox.setSizePolicy(self.sizePolicyB)
        if mode !='v01_single_':
            groupbox.setEnabled(True)
            groupbox.setVisible(False)
            groupbox.setCheckable(True)
            groupbox.setChecked(False)

        gridLayout = QtWidgets.QGridLayout()
        #gridLayout.setContentsMargins(10, 10, 10, 10)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        parent = gridLayout

        self.row, self.column = 0, 3
        #self.insert_checkbox(parent, 'v01_single_collect', cstep=1, alignment=QtCore.Qt.AlignHCenter)
        #self.v01_single_collect.setChecked(True)
        #self.v01_single_collect.setEnabled(False)
        #self.insert_label(parent, text='Collect Data', columnspan=2, cstep=2,rstep=1)

        self.insert_label(parent, text='Search Remotely', rstep=1, cstep=self.column * -1 + 2,
                          tooltip='By checking "Search remotely" one can search on a remote server.\n' +
                                  'For this option, please fill out the server credentials below.')

        # Nanpgraph row
        self.insert_label(parent, text=' - Nanograph Folder', cstep=1,
                          tooltip='Select the folder and the data format in which the nanographs are stored.')
        self.insert_checkbox(parent, mode+'remote_nanograph', cstep=1, alignment=QtCore.Qt.AlignHCenter,
                             logvar=True)
        self.insert_combobox(parent, mode+'filetype_nanographs', ['tif', 'mrc'], cstep=1, logvar=True)
        self.insert_lineedit(parent, mode+'folder_nanographs', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1 + 1, rstep=1, text='Browse',
                               action=self.browse,
                               params=['folder', self.items[self.row][self.column - 1],
                                       self.v01_single_filetype_nanographs, self.items[self.row][self.column - 3]])


        # MDOC LINE
        self.insert_checkbox(parent, 'v01_batch_collect_mdoc', cstep=1, alignment=QtCore.Qt.AlignHCenter,
                             logvar=True)
        self.insert_label(parent, text=' - Mdoc Folder', cstep=1,
                      tooltip='If you hae accompanying mdoc files, select the folder where these files are stored')
        self.insert_checkbox(parent, 'v01_batch_remote_mdoc', cstep=2)
        # self.insert_combobox(parent, 'batch_filtype_mdoc', ['dm4','mrc'], cstep=1)
        self.insert_lineedit(parent, 'v01_batch_folder_mdoc', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1, rstep=1, text='Browse',
                               action=self.browse,
                               params=['folder', self.items[self.row][self.column - 1],
                                       self.items[self.row][self.column - 2],
                                       self.items[self.row][self.column - 3]])
        groupbox.setLayout(gridLayout)
        setattr(self,mode+'gb_collect', groupbox)
        return groupbox

    def createMotioncorGroup(self,mode='v01_single_'):
        groupbox = QGroupBox("Do Motion Correction")
        groupbox.setToolTip('Execute Motion correction using motioncor2')
        groupbox.setVisible(False)
        groupbox.setCheckable(True)
        groupbox.setChecked(False)

        groupbox.setSizePolicy(self.sizePolicyB)

        gridLayout = QtWidgets.QGridLayout()
        gridLayout.setContentsMargins(10, 10, 10, 10)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = gridLayout
        self.row, self.column = 0, 2
        self.insert_label(parent, text='Search Remotely', rstep=1, cstep=self.column * -1,
                          tooltip='By checking "Search remotely" one can search on a remote server.\n' +
                                  'For this option, please fill out the server credentials below.')
        # Motion Correction
    #    self.insert_checkbox(parent, '01_single_motioncorrection', cstep=1, logvar=True)
     #   self.insert_label(parent, text='Do Motion Correction', columnspan=2, rstep=1, cstep = -1*self.column+1,
      #                    tooltip='Execute Motion correction xsusing motioncor2')

        # Gain file
        self.gain_correct = self.insert_checkbox(parent, mode+'gaincorrection', cstep=1, logvar=True, width=20)
        self.insert_label(parent, text=' - Gain File', cstep=1, alignment=QtCore.Qt.AlignLeft,width=120,
                          tooltip='Apply gain correction using the file specified in entry line.')
        self.insert_checkbox(parent, mode+'remote_gainfile', cstep=1, logvar=True, alignment=QtCore.Qt.AlignHCenter)
        self.insert_combobox(parent, mode+'filetype_gaincorrection', ['dm4', 'mrc'], cstep=1)
        self.insert_lineedit(parent, mode+'gainfile', cstep=1, logvar=True)
        self.insert_pushbutton(parent, cstep=self.column * -1 , rstep=1, text='Browse', action=self.browse,
                               params=['file', self.items[self.row][self.column - 1],
                                       self.items[self.row][self.column - 2], self.items[self.row][self.column - 3]])

        # Patch Size
        self.insert_checkbox(parent, mode+'use_patch', cstep=1,width=20,alignment=Qt.AlignHCenter)
        self.insert_label(parent, text=' - Patch Size', cstep=3, alignment=QtCore.Qt.AlignLeft,
                          tooltip='Sets the size of the patches used for motioncorrection (NxN).')
        self.insert_lineedit(parent, mode+'patchsize', cstep=self.column * -1 , rstep=1,
                             value='5', validator=QIntValidator(), logvar=True)

        # FtBin Size
        self.insert_checkbox(parent, mode+'use_ftbin', cstep=1,width=20)
        self.insert_label(parent, text=' - FtBin Factor', cstep=3, alignment=QtCore.Qt.AlignLeft,
                          tooltip='Sets the binning factor for Fourier binning of your tilt images.')
        self.insert_lineedit(parent, mode+'ftbin_factor', cstep=self.column * -1 , rstep=1,
                             value='2', validator=QIntValidator(), logvar=True)
        groupbox.setLayout(gridLayout)
        return groupbox

    def createCredentialsGroup(self,mode=''):
        groupbox = QGroupBox("Credentials Remote Login")
        groupbox.setVisible(False)
        groupbox.setCheckable(True)
        groupbox.setChecked(False)
        gridLayout = QtWidgets.QGridLayout()
        gridLayout.setContentsMargins(10, 10, 10, 10)
        groupbox.setSizePolicy(self.sizePolicyB)

        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = gridLayout
        self.row, self.column = 0, 0

        # SERVERNAME, USERNAME, PASSWORD
        self.insert_label(parent, text=' - Server Name', cstep=3, tooltip='Sets the name of the remote server.')
        self.insert_label(parent, cstep=1)
        self.insert_lineedit(parent, 'servername'+mode, cstep=self.column * -1 , rstep=1, att=True,width=150)

        self.insert_label(parent, text=' - User Name', cstep=3,
                          tooltip='Sets the username used for login into the remote server.')
        self.insert_label(parent, cstep=1)
        self.insert_lineedit(parent, 'username'+mode, cstep=1, att=True,width=150)
        self.insert_label(parent, cstep=self.column * -1, rstep=1, width=20)

        self.insert_label(parent, text=' - Password', cstep=3,
                          tooltip='Sets the password used for login into the remote server.')
        self.insert_label(parent, cstep=1)
        self.insert_lineedit(parent, 'password'+mode, rstep=2, password=True, att=True,width=150)

        groupbox.setLayout(gridLayout)
        return groupbox

    def download(self, params=None):
        

        step = 1
        widget = QWidget(self)
        layout = QHBoxLayout()
        widget.setLayout(layout)

        self.collect, self.motioncor = params[0].isChecked(), params[1].isChecked()

        if self.collect:

            self.progressBar_Collect = QProgressBar()
            self.progressBar_Collect.setSizePolicy(self.sizePolicyA)
            self.progressBar_Collect.setStyleSheet(DEFAULT_STYLE_PROGRESSBAR)
            self.progressBar_Collect.setFormat('Collect Data: %v/%m')
            self.progressBar_Collect.setSizePolicy(self.sizePolicyB)
            layout.addWidget(self.progressBar_Collect)
            self.statusBar.addPermanentWidget(self.progressBar_Collect,stretch=1)

        if self.motioncor:

            self.progressBar_Motioncor = QProgressBar()
            self.progressBar_Motioncor.setSizePolicy(self.sizePolicyA)
            self.progressBar_Motioncor.setStyleSheet(DEFAULT_STYLE_PROGRESSBAR)
            self.progressBar_Motioncor.setFormat('Motion Correct: %v/%m')
            self.progressBar_Motioncor.setSizePolicy(self.sizePolicyB)
            layout.addWidget(self.progressBar_Motioncor)
            self.statusBar.addPermanentWidget(self.progressBar_Motioncor,stretch=1)
            self.completedMC = QLineEdit()
            self.completedMC.setText('0')
            self.completedMC.textChanged.connect(lambda ignore, value=self.completedMC, pb=self.progressBar_Motioncor:
                                           self.update_pb(pb, value))


        self.em_fscollect()

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
                            os.system('cp {} {}'.format(fname, outname))

                    else:
                        try:
                            os.link(fname, outname)
                        except:
                            os.system('cp {} {}'.format(fname, outname))
                    num_files_downloaded += 1

        return num_files_downloaded

    def meta_retrieve(self, ftps, flist, mdoc_list, folder_remote, folder_local, e, value, only_missing=True):
        sleeptime = 0

        while 1:
            num_flist = self.retrieve_files(ftps, flist, folder_local, value, only_missing=True)
            num_mdoc = self.retrieve_files(ftps, mdoc_list, folder_local, value, only_missing=True)

            if self.continuous and num_mdoc == 0 and num_flist == 0:
                time.sleep(30)
                sleeptime += 30

                if (self.search_remote_nanograph) == 1:
                    flist_all = ftps.nlst(self.remote_data_folder)
                else:
                    flist_all = os.listdir(self.remote_data_folder)

                flist = [tiffile for tiffile in flist_all if tiffile[-4:] == self.data_file_type[-4:]]

            if sleeptime > 3600 * 5 or not self.continuous:
                break
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

        guiFunctions.createMetaDataFiles(self.rawnanographs_folder, '')
        self.popup_messagebox("Info", "Completion", 'Successfully finished data transfer')
        self.widgets['CandP'].setEnabled(True)

    def em_fscollect(self,mode='v01_batch_'):

        self.widgets['CandP'].setEnabled(False)
        os.system( 'module load motioncor2/1.2.1' )

        self.c = self.collect
        self.m = self.motioncor
        self.collection_finished = not self.collect
        self.number_collection_subprocesses = 1
        self.number_motioncor_subprocesses = 4
        self.continuous = False #self.on_the_fly_data_retrieval.get()
        self.search_remote_nanograph = self.widgets[mode+'remote_nanograph'].isChecked()
        self.search_remote_dm4 = self.widgets[mode + 'remote_gainfile'].isChecked()
        self.data_file_type = self.widgets[mode+'filetype_nanographs'].currentText()
        self.remote_data_folder = self.widgets[mode+'folder_nanographs'].text()
        self.local_data_folder = self.rawnanographs_folder
        self.dm4_file = self.widgets[mode+'gainfile'].text()
        self.apply_gaincorrection = self.widgets[mode+'gaincorrection'].isChecked()
        self.remote_mdoc_folder = self.widgets[mode+'folder_mdoc'].text()
        self.ftbinfactor = int(self.widgets[mode+'ftbin_factor'].text())
        self.apply_ftbinning = self.widgets[mode+'use_ftbin'].isChecked()
        self.use_patches = self.widgets[mode+'use_patch'].isChecked()
        self.patchsize = int(self.widgets[mode+'patchsize'].text())
        mdoc_list = []
        flist = []
        #self.username, self.servername, self.password = 'emsquare1.science.uu.nl','emuser','#99@3584cg'

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
                    ftps = self.connect_ftp_server(self.servername, self.username, self.password)
                    flist_all = ftps.nlst(self.remote_data_folder)
                    flist = [tiffile for tiffile in flist_all if tiffile[-3:] == self.data_file_type[-3:]]
                except:
                    self.popup_messagebox("Error", "Error", "Remote data folder does not exist.")
                    return

            else:
                try:
                    flist_all = os.listdir(self.remote_data_folder)
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
            if self.widgets[mode+'collect_mdoc'].isChecked():
                try:
                    all_files = ftps.nlst(self.remote_mdoc_folder)
                    mdoc_list = [mdocfile for mdocfile in all_files if mdocfile[-4:] == 'mdoc']
                except:
                    try:
                        all_files = os.listdir(self.remote_mdoc_folder)
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
            self.flist = flist + mdoc_list
            #self.progressBar_Collect.setRange(0, len(self.flist))

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


            self.progressBar_Motioncor.setRange(0,len(motioncor_flist))
            self.motioncor_flist = motioncor_flist

            availableGPU = avail_gpu()

            for i in range(min(len(availableGPU), self.number_motioncor_subprocesses)):
                finish_motioncor = Event()
                self.mcor_procs.append(finish_motioncor)
                proc = Process(target=self.apply_motioncor,args=(motioncor_flist[i::self.number_motioncor_subprocesses],
                                                                 i + 1, e, availableGPU[i],self.completedMC,
                                                                 finish_motioncor))
                procs.append(proc)
                proc.start()

        if self.c or self.m:
            self.widgets['CandP'].setEnabled(False)
            proc = Worker( fn=self.check_run, args=(procs,[e]+self.mcor_procs, self.c, self.m, self.flist,
                                                    self.motioncor_flist, self.continuous, self.rawnanographs_folder,
                                                    self.motioncor_folder) )

            proc.signals.result1.connect(self.update_progress_collect)
            proc.signals.result2.connect(self.update_progress_motioncor)
            proc.signals.finished_collect.connect(self.delete_progressbar_collect)
            proc.signals.finished_mcor.connect(self.delete_progressbar_motioncor)

            proc.start()



    def check_run(self, procs, events, collect, mcor, flist, motioncor_flist, continuous, folder1, folder2, signals):

        if not continuous and (mcor or collect):
            while len(events):

                events = [proc for proc in events if not proc.is_set()]

                a = [os.path.basename(line) for line in glob.glob( os.path.join(folder1, '*')) ]
                b = [os.path.basename(line) for line in glob.glob( os.path.join(folder2, '*')) ]

                if collect:
                    total = 0
                    for f in flist:
                        if os.path.basename(f) in a: total += 1
                    signals.result1.emit(total)

                    if total >= len(flist):
                        signals.finished_collect.emit()

                if mcor:
                    total = 0
                    for f in motioncor_flist:
                        print (f)
                        if os.path.basename(f).replace('.tif','.mrc') in b: total += 1
                    signals.result2.emit(total)
                    if total >= len(motioncor_flist):
                        signals.finished_mcor.emit()
                time.sleep(.1)


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

def main():
    app = QApplication(sys.argv)
    ex = CollectPreprocess()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()
