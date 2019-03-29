import os
import copy

import pyqtgraph as pg

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtWidgets, QtCore, QtGui

from pytom.basic.functions import initSphere
from pytom.basic.files import pdb2em
from pytom.gui.guiStyleSheets import *
from pytom.gui.mrcOperations import *
import pytom.gui.guiFunctions as guiFunctions
import traceback
from numpy import zeros, meshgrid, arange, sqrt
import numpy as np

from scipy.ndimage.filters import gaussian_filter
from ftplib import FTP_TLS, FTP


class BrowseWindowRemote(QMainWindow):
    '''This class creates a new windows for browsing'''
    def __init__(self, parent=None, initdir='/',filter=[''],search='file',credentials=['','',''],outputline='',
                 validate=True):
        super(BrowseWindowRemote, self).__init__(parent)
        self.setGeometry(50, 50, 400, 300)
        self.pathdisplay = QLineEdit()
        self.pathdisplay.setStyleSheet(WHITE)
        self.pathdisplay.setEnabled(False)
        self.splitter0 = splitter0 = QSplitter(Qt.Vertical)
        self.topleft = QListWidget()
        self.topleft.setStyleSheet(WHITE)
        self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)
        self.topright = QListWidget()
        self.topright.setStyleSheet(WHITE)
        if search == 'file': self.topright.itemDoubleClicked.connect(self.select_file)
        self.splitter1 = splitter1 = QSplitter(Qt.Horizontal)
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
            return


        #self.servername,self.username,self.password = 'emsquare1.science.uu.nl','emuser','#99@3584cg'

        print(self.servername,self.username,self.password)

        try:
            self.connect_ftp_server(self.servername, self.username, self.password)
            self.setCentralWidget(splitter0)
            self.success = True
            self.add_folders()

        except:
            if validate:
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


class SelectFiles(BrowseWindowRemote):
    def __init__(self,parent=None, initdir='/',filter=[''],search='file',credentials=['','',''],outputline='',
                 validate=False):
        self.finished = False
        super(SelectFiles,self).__init__(parent,initdir=initdir,filter=[''],search='file',credentials=['','',''],outputline=outputline,
                                     validate=validate)


        self.filters = filter
        self.matchingfiles = []
        self.selectedfiles = []
        if self.outputline.text(): self.selectedfiles += [ll for ll in self.outputline.text().split('\n')]
        self.setCentralWidget(self.splitter0)
        self.add_folders()
        self.show()

    def select_item(self):
        self.outputline.setText('\n'.join(self.selectedfiles))
        self.finished = True
        self.parent().populate_batch_create()

    def select_file(self):
        item = self.topright.currentItem().text()

        for n, file in enumerate(self.selectedfiles):
            if file == item:
                self.selectedfiles = self.selectedfiles[:n]+self.selectedfiles[n+1:]
                break
        self.add_folders()

    def data(self):
        self.subdirs = [os.pardir]
        self.matchingfiles = []
        self.append(os.listdir(self.folderpath))

        matchingfiles = []
        for pat in self.filters:
            matchingfiles += [line for line in self.matchingfiles if line.endswith(pat.split('.')[-1])]
        self.matchingfiles = matchingfiles

    def append(self,line):
        for filename in line:
            folder = os.path.isdir(os.path.join(self.folderpath,filename))

            if folder:
                self.subdirs.append(filename)
            else:
                selected = False
                for file in self.selectedfiles:
                    if filename == os.path.basename(file):
                        selected = True
                        break
                if not selected: self.matchingfiles.append(filename)

        self.matchingfiles = sorted(self.matchingfiles)
        self.subdirs = sorted(self.subdirs)

    def add_folders(self, folders=[]):
        self.data()

        self.topleft.clear()
        self.topright.clear()

        self.topleft.insertItems(0,self.subdirs+self.matchingfiles)
        self.topright.insertItems(0,self.selectedfiles)

    def repopulate_folder_list(self):

        extra = self.topleft.currentItem().text()
        potential = os.path.join(self.folderpath, extra)
        if extra != '..' and os.path.isdir(potential):
            self.folderpath = potential

        elif len(self.folderpath) >1 and extra == '..':
            self.folderpath = os.path.dirname(self.folderpath)
        else:
            self.selectedfiles.append(potential)

        self.pathdisplay.setText(self.folderpath)
        self.add_folders()


class WorkerSignals(QObject):
    '''
    Defines the signals available from a running worker thread.
    Supported signals are:

    finished
        No data

    error
        `tuple` (exctype, value, traceback.format_exc() )

    result
        `object` data returned from processing, anything

    '''

    finished_mcor = pyqtSignal()
    finished_collect = pyqtSignal()
    error = pyqtSignal(tuple)
    result1 = pyqtSignal(object)
    result2 = pyqtSignal(object)


class Worker(QRunnable):
    def __init__(self, fn=None, args=[]):
        #super(ProcessRunnable,self).__init__(parent)
        super(Worker, self).__init__()
        self.signals = WorkerSignals()
        self.fn = fn
        self.args = tuple(list(args)+[self.signals])

    def run(self):
        self.fn(*self.args)

    def start(self):
        QThreadPool.globalInstance().start(self)


class CommonFunctions():
    def __init__(self,parent=None):
        pass

    def size_policies(self):

        self.sizePolicyA = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.sizePolicyA.setHorizontalStretch(1)
        self.sizePolicyA.setVerticalStretch(1)

        self.sizePolicyB = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.sizePolicyB.setHorizontalStretch(1)
        self.sizePolicyB.setVerticalStretch(0)

        self.sizePolicyC = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sizePolicyC.setHorizontalStretch(0)
        self.sizePolicyC.setVerticalStretch(0)

    def add_toolbar(self, decider, new=False, open=True, save=True):
        self.setStyleSheet(MAINC)
        bar = self.menuBar()
        tb = QToolBar()
        tb.setStyleSheet(BARS)
        self.addToolBar(tb)

        if new:
            new = QAction(QIcon("{}/gui/Icons/new_project4.png".format(self.pytompath)), "Create New Project", self)
            tb.addAction(new)

        if open:
            load = QAction(QIcon("{}/gui/Icons/open_project4.png".format(self.pytompath)), "Load Project", self)
            tb.addAction(load)

        if save:
            s = QAction(QIcon("{}/gui/Icons/save_project4.png".format(self.pytompath)), "Save Project", self)
            tb.addAction(s)
            #tb.actionTriggered[QAction].connect(self.save_particleList)

        tb.actionTriggered[QAction].connect(decider)

    def insert_checkbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,
                        alignment=Qt.AlignCenter, tooltip='', logvar=False,width=0):

        widget = QtWidgets.QCheckBox(self)
        widget.setStyleSheet("QCheckBox{background:none;}")
        if logvar:
            if wname in self.logbook.keys():
                widget.setChecked(self.logbook[wname])
            else:
                self.logbook[wname] = False
            widget.setObjectName(wname)
            widget.stateChanged.connect(lambda dummy, name=wname, text='checkbox':
                                        self.update_logbook(name=name, text=text))

        if text:
            widget.setText(text)
        else:
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
        if tooltip: widget.setToolTip(tooltip)
        if width: widget.setFixedWidth(width)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan, alignment)
        setattr(self, wname, widget)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label(self, parent, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,transparent=True,
                     tooltip='', alignment=Qt.AlignLeft, value=None,width=0,sizepolicy=''):
        widget = QtWidgets.QLabel(self)

        if transparent:
            widget.setStyleSheet("QLabel{background:transparent;} \nQToolTip { background-color: white; }")

        if text:
            widget.setText(text)
        elif width == 0:
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)
        if tooltip: widget.setToolTip(tooltip)
        if sizepolicy:widget.setSizePolicy(sizepolicy)
        if width: widget.setFixedWidth(width)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan, alignment)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_combobox(self, parent, wname, labels, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,
                        width=80, tooltip='', logvar=False):
        widget = QtWidgets.QComboBox(self)
        widget.addItems(labels)

        if logvar:
            if wname in self.logbook.keys():
                for n, label in enumerate(labels):
                    if label == self.logbook[wname]:
                        widget.setCurrentIndex(n)
            else:
                self.logbook[wname] = labels[0]
            widget.setObjectName(wname)
            widget.currentIndexChanged.connect(lambda dummy, name=wname, text='combobox':
                                               self.update_logbook(name=name, text=text))

        if tooltip: widget.setToolTip(tooltip)
        widget.setStyleSheet("QLineEdit{background:white; selection-background-color: #1989ac;}")
        widget.setFixedWidth(width)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Maximum, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
        widget.setSizePolicy(sizePolicy)

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget

        self.row += rstep
        self.column += cstep


    def insert_pgSpinbox(self, parent, wname, text='',rowspan=1, columnspan=1, rstep=0, cstep=0, decimals=0,width=0,
                         height=0, tooltip='', logvar=False, value=None, maximum=None, minimum=0, stepsize=1,suffix=''):

        widget = pg.SpinBox(parent=parent)
        self.widgets[wname] = widget

        if maximum: widget.setMaximum(maximum)
        if minimum: widget.setMinimum(minimum)
        if suffix: widget.setSuffix(suffix)
        if stepsize: widget.setSingleStep(stepsize)
        if value: widget.setValue(value)
        if decimals>=0: widget.setDecimals(decimals)
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)

        if logvar:
            try:
                widget.setText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.valueChanged.connect(lambda ignore, name=wname, text='lineedit':
                                       self.update_logbook(name=name, text=text))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_spinbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,width=0, height=0,
                       validator=None, password=False, tooltip='', logvar=False, value=None, att=False, enabled=True,
                       maximum=0, minimum=0, stepsize=0, widgetType=QSpinBox):

        widget = widgetType(self)
        widget.setLocale(QLocale('.'))
        self.widgets[wname] = widget
        widget.setEnabled(enabled)
        widget.setAutoFillBackground(True)
        widget.setStyleSheet("background-color: rgb(255,255,255);")

        widget.setPalette(QPalette(Qt.white ) )
        widget.setKeyboardTracking(False)

        if logvar:
            try:
                widget.setValue(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.valueChanged.connect(lambda ignore, name=wname, text='spinbox':
                                       self.update_logbook(name=name, text=text))

        if password: widget.setEchoMode(QLineEdit.Password)
        if validator: widget.setValidator(validator)
        if tooltip: widget.setToolTip(tooltip)
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)
        if minimum: widget.setMinimum(minimum)
        if maximum: widget.setMaximum(maximum)
        if stepsize: widget.setSingleStep(stepsize)
        if value: widget.setValue(value)
        if att and not logvar: widget.textChanged.connect(lambda ignore, name=wname: self.set_attribute(wname))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)

        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_lineedit(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=0, cstep=0,width=0, height=0,
                        validator=None, password=False, tooltip='', logvar=False, value=None, att=False, enabled=True):
        widget = QtWidgets.QLineEdit(self)

        if logvar:
            try:
                widget.setText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            widget.textChanged.connect(lambda ignore, name=wname, text='lineedit':
                                       self.update_logbook(name=name, text=text))
        if password: widget.setEchoMode(QLineEdit.Password)
        if validator: widget.setValidator(validator)
        if tooltip: widget.setToolTip(tooltip)
        if att and not logvar:
            widget.textChanged.connect(lambda ignore, name=wname:
                                       self.set_attribute(wname))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        setattr(self, wname, widget)


        widget.setStyleSheet("QLineEdit{background:white;}; QLineEdit:disabled{background:white; color:black;} ")
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        if width: widget.setFixedWidth(width)
        if height: widget.setFixedHeight(height)
        if value: widget.setText(str(value))
        widget.setEnabled(enabled)
        self.row += rstep
        self.column += cstep

    def set_attribute(self, wname):
        txt = self.widgets[wname].text()

        if wname[-1] == '1':
            wname = wname[:-1]
            setattr(self, wname, txt)
            self.widgets[wname].setText(txt)
        else:
            setattr(self, wname, txt)
            wname = wname + '1'
            self.widgets[wname].setText(txt)

    def insert_pushbutton(self, parent, text='', rowspan=1, columnspan=1, rstep=0, cstep=0, action=None,
                          tooltip='', width=0, params=['file', '', '', False], iconpath='', wname='', state=True):
        widget = QtWidgets.QPushButton(self)
        if wname:
            print (wname)
            self.widgets[wname] = widget
        widget.setEnabled(state)
        if tooltip: widget.setToolTip(tooltip)
        if text:
            widget.setText(text)
        elif width==0:
            sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
            sizePolicy.setHorizontalStretch(0)
            sizePolicy.setVerticalStretch(0)
            sizePolicy.setHeightForWidth(widget.sizePolicy().hasHeightForWidth())
            widget.setSizePolicy(sizePolicy)

        if action:
            if params:
                widget.clicked.connect(lambda ignore, p=params: action(p))
            else:
                widget.clicked.connect(lambda ignore: action())
        if iconpath:
            icon = QtGui.QPixmap(iconpath)
            widget.setIcon(icon)
        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        if width: widget.setFixedWidth(width)
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_textfield(self, parent, wname, width=500, height=400, logvar=True,
                         rstep=0, cstep=0, rowspan=1, columnspan=1):

        widget = QPlainTextEdit(self)
        widget.setSizePolicy(self.sizePolicyB)
        widget.setStyleSheet("QPlainTextEdit{background:white;}")
        widget.setObjectName(wname)



        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label_spinbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=1, cstep=-1, width=0, height=0,
                             validator=None, tooltip='', value=None, att=False, enabled=True, maximum=0, minimum=0,
                             wtype=QSpinBox, stepsize=0,logvar=True):

        self.insert_label(parent, text=text, cstep=1, alignment=QtCore.Qt.AlignRight, tooltip=tooltip)
        self.insert_spinbox(parent, wname, validator=validator, width=width, enabled=enabled,
                            maximum=maximum, minimum=minimum, cstep=cstep, rstep=rstep, value=value, widgetType=wtype,
                            stepsize=stepsize, logvar=logvar)

    def insert_label_line_push(self, parent, textlabel, wname, tooltip='', text='', cstep=-2, rstep=1, validator=None,
                               mode='folder', remote=False, pushtext='Browse', width=150, filetype='',action='',
                               initdir=''):

        if action == '':
            action = self.browse


        self.insert_label(parent, text=textlabel, cstep=1, alignment=QtCore.Qt.AlignRight, tooltip=tooltip)
        self.insert_lineedit(parent, wname, cstep=1, logvar=True, text='',validator=validator,width=width,enabled=False)

        if initdir:
            params = [mode, self.widgets[wname], filetype, remote, initdir]
        else:
            params = [mode, self.widgets[wname], filetype, remote]

        self.insert_pushbutton(parent, cstep=cstep, rstep=rstep, text=pushtext,width=100,
                               action=action, params=params)

    def insert_label_line(self, parent, textlabel, wname, tooltip='', value='', cstep=-1, rstep=1, validator=None,
                          width=150,logvar=True):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_lineedit(parent, wname, cstep=cstep, rstep=rstep, value=value, logvar=logvar, validator=validator,
                             width=width)

    def insert_label_checkbox(self, parent, textlabel, wname, tooltip='', cstep=-1, rstep=1, logvar=True):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_checkbox(parent, wname, logvar=logvar, cstep=cstep, rstep=rstep,alignment=Qt.AlignLeft)

    def insert_label_combobox(self, parent, textlabel, wname, labels, rowspan=1, columnspan=1, rstep=1, cstep=-1,
                        width=150, tooltip='', logvar=False ):
        self.insert_label(parent, text=textlabel, cstep=1, alignment=Qt.AlignRight, tooltip=tooltip)
        self.insert_combobox(parent, wname, labels, rowspan=rowspan, columnspan=columnspan, rstep=rstep, cstep=cstep,
                             width=width, logvar=logvar)

    def insert_label_action_label(self, parent, actiontext, action='', cstep=0, rstep=1, sizepolicy='', params=''):
        self.insert_label(parent,alignment=Qt.AlignRight,sizepolicy=sizepolicy,rstep=1)
        self.insert_pushbutton(parent, text=actiontext, rstep=rstep, cstep=cstep, action=action, params=params)
        self.insert_label(parent,sizepolicy=sizepolicy,rstep=1)

    def insert_gen_text_exe(self, parent, mode, gen_action='', action='', paramsAction=[], paramsXML=[], paramsCmd=[],
                            paramsSbatch={}, xmlfilename='', exefilename='exe.sh', jobfield=False):

        self.insert_label_action_label(parent, 'Generate command', cstep=1, rstep=-1, sizepolicy=self.sizePolicyB,
                                       action=self.gen_action,
                                       params=[[mode + 'XMLText']+paramsXML,
                                               [mode+'CommandText'] + paramsCmd,
                                               paramsSbatch
                                               ])


        self.insert_checkbox(parent,mode + 'queue',text='sbatch',cstep=-3,rstep=1,logvar=True,alignment=Qt.AlignLeft)

        if jobfield:
            self.insert_textfield(parent, mode + 'XMLText', columnspan=3, rstep=1, cstep=2, width=600,logvar=True)
            self.insert_label(parent, alignment=Qt.AlignRight, rstep=1, cstep=-2, sizepolicy=self.sizePolicyB)
        self.insert_textfield(parent, mode + 'CommandText', columnspan=3, rstep=1, cstep=2, width=600, logvar=True)
        self.insert_label_action_label(parent, 'Execute command', rstep=1, action=self.exe_action,
                                       params=[exefilename, mode+'CommandText', xmlfilename, mode+'XMLText', action,
                                               paramsAction])

    def exe_action(self, params):

        if params[4]: params[4](params[5])

        try:
            # Check if one needs to write pytom related XML file.
            if params[2]:
                jobfile = open(params[2],'w')
                jobfile.write(self.widgets[params[3]].toPlainText())
                jobfile.close()

            # Write executable file
            if len(params[0][0]) > 1:
                params[0] = os.path.join(self.widgets[params[0][0]].text(), params[0][1])
            exefile = open(params[0], 'w')
            exefile.write(self.widgets[params[1]].toPlainText())
            exefile.close()

            if len(self.widgets[params[1]].toPlainText().split('SBATCH') ) > 2:
                os.system('sbatch {}'.format(params[0]))
            else:
                os.system('sh {}'.format(params[0]))
        except:
            print ('Please check your input parameters. They might be incomplete.')

    def gen_action(self, params):
        mode = params[1][0][:-len('CommandText')]

        for i in range(2):
            if params[i][0] in self.widgets.keys():
                text = params[i][-1]
                if params[i][1:-1]:
                    d = []
                    for a in params[i][1:-1]:

                        if a in self.widgets.keys():
                            datatype = type(self.widgets[a])
                            if datatype in (QSpinBox, QDoubleSpinBox, QLineEdit): d.append(self.widgets[a].text())
                            elif datatype == QComboBox: d.append(self.widgets[a].currentText())
                        elif type(str(a)) == type(''):
                            d.append(a)

                    text = text.format( d=d )
                if i==0: self.widgets[params[i][0]].setPlainText(text)
        # Check if user wants to submit to queue. If so, add queue header.
        if self.widgets[mode+'queue'].isChecked():
            d = params[2]
            folder = d['folder']

            if type(folder) == type([]):
                folder = os.path.dirname( os.path.join(self.widgets[mode + 'tomofolder'].text(),
                                                       self.widgets[folder[0]].text(), folder[1]) )


            text = guiFunctions.gen_queue_header(d['fname'], folder, d['cmd'], d['modules']) + text

        self.widgets[params[i][0]].setPlainText(text)

    def browse(self, params):
        mode = params[0]
        line = params[1]
        filetype = params[2]

        try:
            remote = params[3].isChecked()
        except:
            remote = False


        if line.text():
            if not remote:
                initdir = os.path.join('/', line.text())
            else:
                initdir = line.text()
        else:
            try:
                initdir = params[4]
            except:
                try:
                    initdir = self.projectname
                except:
                    initdir='./'


        try: initdir = initdir.text()
        except: pass


        if filetype:
            try:
                filters = "Image files (*.{})".format( filetype.currentText() )
            except:
                if type(filetype) == str:
                    filetype=[filetype]
                    filters = "Image files (*.{})".format( filetype )
                if type(filetype) == list:
                    filters = ''
                    for f in filetype:
                        filters += "{} files (*.{});; ".format(f, f)
                    filters = filters[:-3]
        else:
            filters = "Image files (*.*)"

        if not remote:
            if mode == 'folder':
                line.setText(QFileDialog.getExistingDirectory(self, 'Open file', initdir))
            if mode == 'file':

                filename = QFileDialog.getOpenFileName(self, 'Open file', initdir, filters)
                if str(filename[0]): line.setText(str(filename[0]))
        if remote:
            credentials = [self.servername, self.username, self.password]
            if not initdir: initdir = '/'
            filters = ['*.{}'.format(filetype.currentText())]

            b = BrowseWindowRemote(self, initdir=initdir, search=mode, filter=filters, outputline=line,
                                   credentials=credentials)

    def update_logbook(self,name='',text=''):

        if name:
            if text == 'checkbox':
                self.logbook[name] = self.widgets[name].isChecked()
            if text == 'lineedit':
                self.logbook[name] = self.widgets[name].text()
            if text == 'combobox':
                self.logbook[name] = self.widgets[name].currentText()
            elif text == 'spinbox':
                self.logbook[name] = self.widgets[name].value()

    def create_expandable_group(self, action, sizeP, text, mode=''):
        a = QCheckBox(text=text)
        a.setSizePolicy(sizeP)
        b = action(mode=mode)
        a.clicked.connect(lambda ignore,  w=a, w2=b: self.toggle_groupbox_visible(w, w2))
        b.clicked.connect(lambda ignore, w2=a,  w=b: self.toggle_groupbox_visible(w, w2))
        return a, b

    def toggle_groupbox_visible(self,widget,widget2):
        widget.setVisible(False)
        widget2.setChecked(True==widget.isChecked())
        widget2.setVisible(True)

    def popup_messagebox(self, messagetype, title, message):
        if messagetype == 'Error':
            QMessageBox().critical(self, title, message, QMessageBox.Ok)

        elif messagetype == 'Warning':
            QMessageBox().warning(self, title, message, QMessageBox.Ok)

    def update_pb(self, pb, value):
        print ( 'update: ', value.text() )
        pb.setValue( int(value.text()) )

    def fill_tab(self, id, headers, types, values, sizes, tooltip=[],wname='v02_batch_aligntable_', connect=0):
        try:
            self.tables[id].setParent(None)
            self.pbs[id].setParent(None)
            self.ends[id].setParent(None)
        except:

            pass

        self.tables[id] = SimpleTable(headers, types, values, sizes, tooltip=tooltip,connect=connect)
        self.widgets['{}{}'.format(wname,id)] = self.tables[id]
        self.pbs[id] = QPushButton('Run')
        self.pbs[id].setSizePolicy(self.sizePolicyC)
        self.ends[id] = QWidget()
        self.ends[id].setSizePolicy(self.sizePolicyA)

        for a in (self.tables[id], self.pbs[id], self.ends[id]):
            self.table_layouts[id].addWidget(a)

    def create_groupbox(self,title,tooltip,sizepol):
        groupbox = QGroupBox(title)
        groupbox.setToolTip(tooltip)
        groupbox.setCheckable(True)
        groupbox.setChecked(False)
        groupbox.setSizePolicy(sizepol)
        parent = QtWidgets.QGridLayout()
        groupbox.setLayout(parent)
        return groupbox, parent


class CreateMaskFile(QMainWindow, CommonFunctions):
    def __init__(self,parent=None,maskfname=''):
        super(CreateMaskFile,self).__init__(parent)
        w = QWidget(self)
        l = QGridLayout()

        self.logbook = self.parent().logbook
        w.setLayout(l)
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = l
        self.insert_label_line(parent,'Size of Template',wname='size_template',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent, 'Radius Spherical Mask', wname='radius',rstep=1, cstep=-1,
                               validator=QIntValidator())
        self.insert_label_spinbox(parent, 'smooth_factor', text='Smoothing Factor', rstep=1, cstep=0, value=0,
                                  wtype=QDoubleSpinBox,stepsize=0.1)
        self.insert_pushbutton(parent, 'Create', action=self.generate_mask,
                               params=['size_template', 'radius','smooth_factor',maskfname])
        self.setCentralWidget(w)
        self.show()

    def generate_mask(self,params):
        radius = int(self.widgets['radius'].text())
        smooth = float(self.widgets['smooth_factor'].value())
        size = int(self.widgets['size_template'].text())
        fname = os.path.join(self.parent().frmdir, 'FRM_mask.em')
        maskfilename = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.em')[0])
        if maskfilename and not maskfilename.endswith('.em'): maskfilename += '.em'
        try:
            initSphere(size, size, size, radius=radius, smooth=smooth, filename=maskfilename)
            self.parent().widgets[params[-1]].setText(maskfilename)
        except:
            self.popup_messagebox('Error','Mask Generation Failed', 'Generation of the mask failed. Please select an existing mask, or generate a mask yourself.')

        self.close()


class ConvertEM2PDB(QMainWindow, CommonFunctions):
    def __init__(self,parent,emfname='',folder='./'):
        super(ConvertEM2PDB, self).__init__(parent)
        self.folder = folder
        w = QWidget(self)
        l = QGridLayout()
        self.logbook = self.parent().logbook
        w.setLayout(l)
        self.widgets = {}
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = l

        self.insert_label_line(parent, 'PDB ID', wname='pdb_id', cstep=-1, rstep=1,width=150)

        self.insert_label_line(parent, 'Chain', wname='chain', rstep=1, cstep=-1, value='A')

        self.insert_label_spinbox(parent, 'volume_size', text='Volume size', rstep=1, cstep=-1,
                                  stepsize=4,minimum=8,maximum=4096,value=128,
                                  tooltip='Volume length(size) in all dimensions')

        self.insert_label_spinbox(parent, 'pixel_size', text='Pixel Size', rstep=1, cstep=0, value=3.5,
                                  wtype=QDoubleSpinBox, stepsize=0.1, minimum=1.0, maximum=100,
                                  tooltip='Pixel size of output volume (in Angstrom)')

        self.insert_pushbutton(parent, 'Create', action=self.pdb2em,
                               params=['size_template', 'radius', 'smooth_factor', emfname])

        self.setCentralWidget(w)
        self.show()

    def pdb2em(self,params):
        pdbid = self.widgets['pdb_id'].text()
        if not pdbid: return
        chain = self.widgets['chain'].text()
        cube_size = int(self.widgets['volume_size'].text())
        pixel_size = float(self.widgets['pixel_size'].value())
        os.system('getpdb {} {}'.format(pdbid, self.folder))
        out_fname = str(QFileDialog.getSaveFileName(self, 'Save EM model.', self.folder, filter='*.em')[0])
        if not out_fname: return
        if not out_fname.endswith('.em'): out_fname += '.em'
        pdb2em('{}/{}.pdb'.format(self.folder,pdbid), pixel_size, cube_size, chain=chain, fname=out_fname)
        self.parent().widgets[params[-1]].setText(out_fname)
        self.close()


class MyCircleOverlay(pg.EllipseROI):
    def __init__(self, pos, size, label='', **args):
        pg.ROI.__init__(self, pos, size, **args)
        self.aspectLocked = True
        colorsHEX = ['00ffff', 'f5f5dc', '0000ff', 'a52a2a', '7fff00', 'd2691e', 'daa520', 'ff7f50', '00ffff',
                     'dc143c', '00008b', '006400', '7fffd4', 'ff00ff', 'ffd700', '008000', '4b0082', 'f0e68c', 'add8e6',
                     '90ee90', 'ff00ff', '800000', '808000', 'ffa500', 'ff4500', 'da70d6', 'ffc0cb', '800080', 'dda0dd',
                     'fa8072', 'fa8072', '008080', 'ffff00', '40e0d0', '9acd32',]*15

        html = '<div style="text-align: center"><span style="color: #{}; font-size: 20pt;">{}</span></div>'

        if label:
            self.label = QtGui.QGraphicsTextItem(label, self)
            #for d in dir(self.label): print (d)
            self.label.setHtml(html.format(colorsHEX[int(label)], int(label)))
            self.label.setPos(QtCore.QPoint( self.boundingRect().center().x() - (self.label.boundingRect().width()/2)*0,
                                             self.state['size'][1]/1.5 ))
            #self.label = pg.TextItem(text=label,anchor=(),angle=180)


def circle(pos, size = 2, label='',color=Qt.blue):
    pensize = 0.02*40./size
    return MyCircleOverlay(pos=pos, pen=QtGui.QPen(color, pensize), size=size, movable=False, label=label)


class SimpleTable(QMainWindow):

    def __init__(self, headers, types, values, sizes=[], tooltip=[] , connect=0):


        super(SimpleTable, self).__init__()
        central_widget = QWidget(self)
        #self.setGeometry(10, 10, 700, 300)  # Create a central widget
        grid_layout = QGridLayout(self)
        grid_layout.setContentsMargins(0,0,0,0) # Create QGridLayout
        central_widget.setLayout(grid_layout)  # Set this layout in central widget
        table = QTableWidget(self)
        table.setColumnCount(len(headers))  # Set three columns
        table.setRowCount(len(values))  # and one row
        header = table.horizontalHeader()
        self.widgets = {}
        # Set the table headers
        table.setHorizontalHeaderLabels(headers)
        table.setVerticalHeaderLabels(['', ] * len(values))
        table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        #table.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # Set the tooltips to headings
        self.headers = headers
        hh = table.horizontalHeader()
        table.verticalHeader().hide()
        # Set the alignment to the headers
        for i in range(len(headers)):

            if i+1 == len(headers):
                hh.setResizeMode(i, QtGui.QHeaderView.Stretch)
            elif sizes[i] == 0:
                hh.setSectionResizeMode(i,QHeaderView.ResizeToContents)
            else:
                hh.setSectionResizeMode(i,sizes[i])

            if i+1 != len(headers):
                table.horizontalHeaderItem(i).setTextAlignment(Qt.AlignHCenter)
            else:
                table.horizontalHeaderItem(i).setTextAlignment(Qt.AlignLeft)

            if len(tooltip) > i: table.horizontalHeaderItem(i).setToolTip(tooltip[i])

            for v in range(len(values)):
                # Fill the first line
                if types[i] == 'txt':
                    widget = QWidget()
                    data= "{}".format(values[v][i].split('/')[-1])
                    cb = QLabel( data )
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(10, 0, 10, 0)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v, i)] = cb
                    widget.setStyleSheet('background: white;')

                elif types[i] == 'checkbox' and values[v][i]:
                    widget = QWidget()
                    cb = QCheckBox()
                    if connect: cb.stateChanged.connect(lambda d, ID=id: connect(ID))
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    if i+1==len(headers): layoutCheckBox.setAlignment(Qt.AlignLeft)
                    layoutCheckBox.setContentsMargins(10, 0, 10, 0)
                    table.setCellWidget(v,i,widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = cb
                    widget.setStyleSheet('background: white;')
                elif types[i] == 'lineedit':
                    widget = QWidget()
                    layoutCheckBox = QHBoxLayout(widget)
                    le = QLineEdit()
                    le.setFixedWidth(table.columnWidth(i))
                    le.setText(str(values[v][i]))
                    layoutCheckBox.addWidget(le)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = le
                    widget.setStyleSheet('background: white;')

                elif types[i] == 'spinbox':
                    widget = QSpinBox()
                    widget.setValue(values[v][i])

                elif types[i] == 'combobox':
                    widget = QWidget()
                    cb = QComboBox()
                    l = QVBoxLayout(widget)
                    l.addWidget(cb)
                    cb.setContentsMargins(0, 0, 0, 0)
                    l.setContentsMargins(0, 0, 0, 0)
                    for t in values[v][i]:
                        cb.addItem(t)
                    table.setCellWidget(v, i, widget)
                    self.widgets['widget_{}_{}'.format(v,i)] = cb
                    widget.setStyleSheet('background: white; selection-background-color: #1989ac;')
                else:
                    widget = QWidget()
                    cb = QLabel('')
                    layoutCheckBox = QHBoxLayout(widget)
                    layoutCheckBox.addWidget(cb)
                    layoutCheckBox.setAlignment(Qt.AlignCenter)
                    layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                    table.setCellWidget(v, i, widget)
                    widget.setStyleSheet('background: white;')

        #if not sizes: table.resizeColumnsToContents()
        self.sizePolicyC = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sizePolicyC.setHorizontalStretch(0)
        self.sizePolicyC.setVerticalStretch(0)
        self.sizePolicyB = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.sizePolicyB.setVerticalStretch(0)



        self.types = types
        self.general_widgets = []
        self.table2 = QTableWidget()
        self.table2.setColumnCount(len(headers))  # Set three columns
        self.table2.setRowCount(1)
        self.table2.verticalHeader().hide()
        self.table2.horizontalHeader().hide()
        self.table2.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table2.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        applyOptionSum = 0

        for n, t in enumerate(types):
            if t == 'checkbox':
                widget = QWidget()
                cb = QCheckBox()
                l = QVBoxLayout(widget)
                l.addWidget(cb)
                #cb.setStyleSheet('QCheckBox::indicator{ width: 10px; height: 10px; }')
                cb.setFixedWidth(table.columnWidth(n))
                widget.setSizePolicy(self.sizePolicyC)
                cb.setToolTip('{} All'.format( headers[n].capitalize() ))
                #l.setAlignment(Qt.AlignCenter)
                cb.setContentsMargins(10, 0, 10, 0)
                l.setContentsMargins(10, 0, 10, 0)
                #widget.setStyleSheet('background: white;')
                #glayout.addWidget(cb)
                self.general_widgets.append(cb)
                self.table2.setCellWidget(0, n, widget)
                cb.stateChanged.connect(lambda dummy, rowIndex=n, c=t: self.on_changeItem(rowIndex,c))
                applyOptionSum +=1
            elif t == 'combobox':
                widget = QWidget()
                cb = QComboBox()
                l = QVBoxLayout(widget)
                l.addWidget(cb)

                cb.setContentsMargins(0, 0, 0, 0)
                l.setContentsMargins(0, 0, 0, 0)
                cb.setFixedWidth(table.columnWidth(n))
                for value in values[0][n]:
                    cb.addItem(value)
                #glayout.addWidget(cb)
                widget.setStyleSheet('selection-background-color: #1989ac;')
                self.general_widgets.append(cb)
                self.table2.setCellWidget(0, n, widget)
                cb.currentTextChanged.connect(lambda dummy, rowIndex=n, c=t: self.on_changeItem(rowIndex, c))
                applyOptionSum +=1

            elif t == 'lineedit':
                widget = QWidget()
                layoutCheckBox = QHBoxLayout(widget)
                le = QLineEdit()
                le.setText(str(values[0][n]))
                layoutCheckBox.addWidget(le)
                layoutCheckBox.setAlignment(Qt.AlignCenter)
                layoutCheckBox.setContentsMargins(0, 0, 0, 0)
                self.table2.setCellWidget(0, n, widget)
                self.general_widgets.append(le)
                le.textChanged.connect(lambda dummy, rowIndex=n, c=t: self.on_changeItem(rowIndex, c))
            else:
                widget = QWidget()
                a = QLabel('')
                l = QVBoxLayout(widget)
                l.setAlignment(Qt.AlignCenter)
                if n == 0:
                    a.setText('Apply to all')
                    l.setAlignment(Qt.AlignCenter)
                l.addWidget(a)
                a.setFixedWidth(table.columnWidth(n))
                widget.setSizePolicy(self.sizePolicyC)
                l.setContentsMargins(0, 0, 0, 0)
                self.table2.setCellWidget(0, n, widget)
                self.general_widgets.append(a)

        #glayout.addWidget(self.table2)
        self.table = table
        self.table2.setMaximumHeight(self.table2.rowHeight(0))
        self.table.setMinimumHeight(min(400, self.table2.rowHeight(0)*(len(values)+.8)))

        if applyOptionSum: grid_layout.addWidget(self.table2, 0, 0)
        grid_layout.addWidget(table, 1, 0)
        #self.table2.horizontalHeaderItem(len(headers) - 1).setTextAlignment(Qt.AlignLeft)

        for i in range(self.table.columnCount()):
            self.table2.setColumnWidth(i, self.table.columnWidth(i))
        self.setCentralWidget(central_widget)

    def on_changeItem(self, rowIndex,widgetType):
        for i in range(self.table.rowCount()):
            if 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='combobox':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setCurrentText(self.general_widgets[rowIndex].currentText())
            elif 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='checkbox':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setChecked(self.general_widgets[rowIndex].isChecked())
            elif 'widget_{}_{}'.format(i,rowIndex) in self.widgets.keys() and widgetType=='lineedit':
                self.widgets['widget_{}_{}'.format(i,rowIndex)].setText(self.general_widgets[rowIndex].text())

        self.table.horizontalHeader().setResizeMode(len(self.headers)-1, QtGui.QHeaderView.Stretch)

        for i in range(self.table.columnCount()):
            self.table2.setColumnWidth(i, self.table.columnWidth(i))


class GuiTabWidget(QWidget, CommonFunctions):
    def __init__(self, parent=None, headers=[],offx=0,offy=0,dimx=900,dimy=721,logbook=[]):

        super(GuiTabWidget, self).__init__(parent)

        self.addTabs(headers=headers, offx=offx, offy=offy, dimx=dimx, dimy=dimy,soff=50)

    def addTabs(self, headers, widget=QWidget, subheaders=[], offx=0,offy=0,dimx=900,dimy=721,soff=0):
        self.size_policies()
        self.scrollarea = QScrollArea(self)
        if soff: self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        #self.scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scrollarea.setSizePolicy(self.sizePolicyA)
        self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setContentsMargins(0,0,0,0)
        #self.scrollarea.setGeometry(offx, offy, dimx, dimy)

        self.scrollarea.setGeometry(0,0,900,700)

        self.tabWidget = QtWidgets.QTabWidget(self.scrollarea)
        self.tabWidget.setContentsMargins(0,0,0,0)
        self.tabWidget.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)

        self.scrollarea.setWidget(self.tabWidget)
        self.tabs = []

        try:
            self.scrollareas.append(self.scrollarea)
            self.scrolloffset.append(soff)
        except:
            self.scrollareas=[self.scrollarea]
            self.scrolloffset = [soff]


        try:
            self.logbook = self.parent().logbook
            self.widgets = {}
            self.statusBar = self.parent().sbar
            self.scrollarea.setStyleSheet('background: #{};'.format(self.parent().middlec))
        except:
            pass


        for n, header in enumerate(headers):
            if len(subheaders) >= len(headers) and len(subheaders[n]) >1:
                subtab = widget(headers=subheaders[n],dimx=dimx-50)
                self.scrollareas += subtab.scrollareas
                self.scrolloffset += subtab.scrolloffset
                tab = QWidget(self.tabWidget)
                tab.setSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
                layout = QVBoxLayout()
                layout.addWidget(subtab)
                tab.setLayout(layout)
                for m, subheader in enumerate(subheaders[n]):
                    setattr(self, 'tab{}{}'.format(n + 1, m + 1), subtab.tabs[m])

            else:
                tab = QWidget()
            self.tabs.append(tab)
            tab.setObjectName(header)

            setattr(self,'tab{}'.format(n+1),tab)
            self.tabWidget.addTab(tab, header)


class KeyPressGraphicsWindow(pg.GraphicsWindow):
    sigKeyPress = QtCore.pyqtSignal(object)

    def __init__(self, *args, **kwargs):
        super(KeyPressGraphicsWindow,self).__init__(*args, **kwargs)

    def keyPressEvent(self, ev):
        self.scene().keyPressEvent(ev)
        self.sigKeyPress.emit(ev)


class ParticlePicker(QMainWindow, CommonFunctions):
    def __init__(self,parent=None, fname=''):
        super(ParticlePicker,self).__init__(parent)
        self.size_policies()
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.layout = QGridLayout(self)
        self.cw = QWidget(self)
        self.cw.setSizePolicy(self.parent().sizePolicyB)
        self.cw.setLayout(self.layout)
        self.setCentralWidget(self.cw)
        self.setGeometry(0, 0, 1000, 1000)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.add_toolbar(self.open_load)

        self.radius = 20
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0,0)

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(250, 750), border=True)
        self.leftimage  = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(750, 750), border=True)
        self.centimage  = w.addViewBox(row=0, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(750, 250), border=True)
        self.bottomimage  = w2.addViewBox(row=0, col=0 )
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]



        self.layout.addWidget(w1,0,0)
        self.layout.addWidget(w,0,1)
        self.layout.addWidget(w2,1,1)
        self.layout.addWidget(self.operationbox,1,0)
        self.title = parent.widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle( "Manual Particle Selection From: {}".format(self.title ) )
        self.centcanvas.wheelEvent = self.wheelEvent
        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        #self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centcanvas.sigMouseReleased.connect(self.empty)
        self.load_image()
        self.leftimage.setXRange(0, self.vol.shape[0])

        self.add_controls(self.layout_operationbox)

        self.subtomo_plots = PlotterSubPlots(self)
        self.subtomo_plots.show()
        pg.QtGui.QApplication.processEvents()

    def wheelEvent(self, event):

        step = event.angleDelta().y()/120
        increment = int(self.widgets['step_size'].text())*step
        if self.slice+increment < self.vol.shape[0] and self.slice+increment > -1:
            self.update_circles(increment)
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked()==False)

        if Qt.Key_D == evt.key():
            w = self.widgets['delete_item']
            w.setChecked(w.isChecked()==False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = -1*int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.subtomo_plots.close()
            self.close()

    def add_controls(self,prnt):
        vDouble = QtGui.QDoubleValidator()
        vInt = QtGui.QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows


        self.insert_label(prnt, text='Size Selection: ',cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ',cstep=-2,rstep=1, value=self.radius*2, minimum=1, maximum=100)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_checkbox(prnt, 'delete_item', cstep=1)
        self.insert_label(prnt, text='Delete Particle', alignment=Qt.AlignLeft,rstep=1)

        self.insert_label(prnt, sizepolicy=self.sizePolicyB, rstep=1, cstep=0)

        self.insert_label(prnt,text='Step Size',cstep=1,alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt,wname='step_size',cstep=-2, value=10, rstep=1,
                            minimum=1, maximum=int(self.vol.shape[0]/4))

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Apply Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt,'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.')
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def open_load(self, q):
        if q.text() == 'Save Project': self.save_particleList()
        elif q.text() == 'Load Project': self.load_particleList()

    def save_particleList(self):

        ext = '.'+os.path.basename(self.title).split('.')[-1]
        fname = os.path.join(self.parent().pickpartfolder,'coords_' + os.path.basename(self.title.replace(ext,'')))

        fname = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.txt')[0])

        if fname:
            out = open(fname, 'w')
            out.write('#FNAME \t\t{}\n'.format(os.path.basename(self.title)))

            folder = os.path.dirname(self.title)
            if 'INFR' in os.path.basename(self.title):
                key = 'INFR'
            else:
                key='WBP'

            inputJobName = "cat {}/{}_reconstruction.sh | grep 'referenceMarkerIndex' | awk '{print $2}'"

            try:
                refMarkIndex = os.popen(inputJobName.format(folder,key)).read()[:-1]
            except:
                refMarkIndex = 1

            out.write('#MARKERINDEX \t{}\n'.format(refMarkIndex))

            for x, y, z in self.particleList:
                outtext =  '{:8.0f} {:8.0f} {:8.0f}\n'.format(x,y,z)
                out.write(outtext)

    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "Coordinate files (*.{})".format(filetype))[0])
        if not filename: return

        particlePos = [map(float, line.split()) for line in open(filename).readlines() if line != '' and not '#' in line]

        self.remove_all()
        self.slice = int(self.vol.shape[0]/2)

        for x,y,z in particlePos:

            self.add_points(self.pos, int(round(x)), int(round(y)), int(round(z)), self.radius, add=True)
        self.replot()
        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

        pass

    def empty(self):
        pass

    def stateGaussianChanged(self):
        w = self.widgets['apply_gaussian_filter']
        width = self.widgets['width_gaussian_filter'].text()

        if w.isChecked():
            if len(width) > 0 and abs(self.current_width - float(width)) > 0.01:
                self.vol = self.volg = gaussian_filter(self.backup, float(width))
                self.current_width = float(width)
            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def sizeChanged(self):
        a = self.widgets['size_selection: '].text()
        if not a: return


        plist = copy.deepcopy(self.particleList)
        self.remove_all()

        self.radius = int(self.widgets['size_selection: '].text()) / 2


        for nn, (x,y,z) in enumerate(plist):
            self.add_points(QPoint(0,0), x, y, z, self.radius, add=True)

        self.subtomo_plots.size_subtomo = self.radius*2

        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=1))
        self.img1b.setImage(image=self.vol.sum(axis=2).T)

    def replot(self):
        crop = self.vol[int(self.slice), :, :]
        self.img1m.setImage(image=crop.T)

        #self.centcanvas.removeItem(self.hist)
        #self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(numpy.median(crop)-crop.std()*3,numpy.median(crop)+crop.std()*3)
        #self.centcanvas.addItem(self.hist)

    def mouseHasMoved(self, evt):

        try:
            pos = self.centimage.mapSceneToView( evt.scenePos() )
        except:
            return

        add = True
        remove = self.widgets['delete_item'].isChecked() or evt.button()==2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x,y,z) in enumerate(self.particleList):
            if sqrt( (x-pos.x())**2 + (y-pos.y())**2 + (z-self.slice)**2 ) < self.radius:
                add = False
                if remove:
                    self.remove_point(n-num_deleted_items,z)
                    self.subtomo_plots.delete_subplot([x,y,z])
                    num_deleted_items += 1

        if add and not remove:
            X, Y = pos.x(), pos.y()

            self.add_points(pos, X, Y, self.slice, self.radius,add=add)
            self.subtomo_plots.add_subplot(self.vol, self.particleList[-1])

    def add_points(self, pos, cx, cy, cz, radius, add=False):
        self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz))])

        pos.setX( cx-radius)
        pos.setY( cy-radius)
        self.circles_cent.append(circle(pos, size=(radius)*2))

        if abs(cz - self.slice) < self.radius:
            self.centimage.addItem(self.circles_cent[-1])

        pos.setX(cx - self.radius)
        pos.setY(cz - self.radius)
        self.circles_bottom.append(circle(pos, size=self.radius * 2))
        if add:
            self.bottomimage.addItem(self.circles_bottom[-1])

        pos.setX(cz - self.radius)
        pos.setY(cy - self.radius)
        self.circles_left.append(circle(pos, size=self.radius * 2))
        if add:
            self.leftimage.addItem(self.circles_left[-1])

    def remove_point(self, n, z, from_particleList=True):
        if from_particleList:
            self.particleList = self.particleList[:n] + self.particleList[n + 1:]

        if abs(z - self.slice) <= self.radius:
            self.centimage.removeItem(self.circles_cent[n])
        self.circles_cent = self.circles_cent[:n] + self.circles_cent[n + 1:]

        self.leftimage.removeItem(self.circles_left[n])
        self.circles_left = self.circles_left[:n] + self.circles_left[n + 1:]

        self.bottomimage.removeItem(self.circles_bottom[n])
        self.circles_bottom = self.circles_bottom[:n] + self.circles_bottom[n + 1:]

    def remove_all(self):
        for image in self.image_list:
            for child in image.allChildren():
                if type(child) == MyCircleOverlay:
                    image.removeItem(child)

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.particleList = []

    def load_image(self):
        if not self.title: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol  = read(self.title)
            self.vol = copy.deepcopy( vol2npy(vol) )
            self.vol = self.vol.T
        elif self.title.endswith('mrc'):
            self.vol = read_mrc(self.title)

        #self.vol[self.vol < -4.] = -4.
        self.backup = self.vol.copy()
        self.origin = self.vol.copy()

        self.vol[ self.vol < self.vol.min()] = self.vol.min()

        self.dim = self.vol.shape[0]
        self.slice = self.d = int(self.dim / 2)


        self.img1a = pg.ImageItem(self.vol.sum(axis=1))
        self.img1b = pg.ImageItem(self.vol.sum(axis=2))
        self.img1m = pg.ImageItem(self.vol[int(self.slice), :, :])

        self.leftimage.addItem(self.img1a)
        self.centimage.addItem(self.img1m)
        self.bottomimage.addItem(self.img1b)

        self.leftcanvas.setAspectLocked(True)
        self.bottomcanvas.setAspectLocked(True)

        self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.centcanvas.addItem(self.hist)

        self.replot_all()

    def update_circles(self, update=0):
        plist = copy.deepcopy(self.particleList)
        self.remove_all()

        for n, (cx, cy, cz) in enumerate(plist):
            radius = self.radius
            if abs(cz - self.slice) < self.radius:
                radius = sqrt(self.radius ** 2 - (cz - self.slice)**2)
            self.add_points(self.pos, cx, cy, cz, radius,add=True)
        self.slice += update

class GeneralSettings(QMainWindow, CommonFunctions):
    def __init__(self,parent):
        super(GeneralSettings, self).__init__(parent)
        self.stage='generalSettings_'
        self.pytompath = self.parent().pytompath


        #headers = ['Data Transfer', 'Tomographic Reconstruction', 'Particle Picking', 'Subtomogram Analysis', "Job Submission"]
        #subheaders  = [[],]*len(headers)
        #self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.setGeometry(50, 50, 300, 100)

        self.cwidget = QWidget()
        self.gridLayout = QGridLayout()
        self.setWindowModality(Qt.ApplicationModal)

        # self.gridLayout.setContentrsMargins(10, 10, 10, 10)

        self.setStyleSheet('background: #{};'.format(self.parent().mainc))
        self.cwidget.setLayout(self.gridLayout)
        self.setCentralWidget(self.cwidget)

        self.show()

class PlotterSubPlots(QMainWindow,CommonFunctions):
    def __init__(self, parent=None, width=400, size_subplot=80, size_subtomo=40, height=1000, offset_x=0, offset_y=0):
        super(PlotterSubPlots,self).__init__(parent)
        self.width = width
        self.height = height
        self.size_subplot = size_subplot
        self.size_subtomo = size_subtomo
        self.dx = offset_x
        self.dy = offset_y

        self.num_subtomo_per_row = round(self.width/self.size_subplot)

        self.size_policies()
        self.superscroll = QScrollArea(self)
        self.superscroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.superscroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.superscroll.setGeometry(self.dx,self.dy,self.width,self.height)

        self.scrollarea = QScrollArea()
        self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setGeometry(self.dx,self.dy,self.width,self.height)
        self.scrollarea.setSizePolicy(self.sizePolicyA)
        self.canvas = KeyPressGraphicsWindow(size=(self.width, self.height), border=None)
        self.scrollarea.setWidget(self.canvas)
        self.superscroll.setWidget(self.scrollarea)
        #self.setCentralWidget(self.scrollarea)
        self.setLayout(QHBoxLayout())
        self.setGeometry(self.dx,self.dy,self.width,self.height)
        self.init_variables()

    def init_variables(self):
        #layout.addWidget(self.canvas)
        self.vBoxList = []
        self.iItemList = []
        self.coordinates = []
        self.assigned = []
        self.index = [0,0]
        self.num_assigned = 0

    def keyPressEvent(self, evt):
        if evt.key() == Qt.Key_Escape:
            self.close()
            self.parent().close()

    def add_subplot(self, tomo, position):
        [x,y,z]  = position
        xmin = max(0, x - int(self.size_subtomo / 2))
        xmax = min(tomo.shape[1], x + int(self.size_subtomo / 2))
        ymin = max(0, y - int(self.size_subtomo / 2))
        ymax = min(tomo.shape[2], y + int(self.size_subtomo / 2))

        subtomo = zeros((int(self.size_subtomo),int(self.size_subtomo)),dtype=float)
        subtomo[:xmax-xmin,:ymax-ymin] = (tomo[position[2]].T)[xmin:xmax, ymin:ymax]

        self.position = position
        if self.index[1] == 0:

            self.newViewBox()
            self.iItemList.append( pg.ImageItem(subtomo))
            self.vBoxList[-1].addItem(self.iItemList[-1])
            self.vBoxList[-1].setRange(xRange=[0, self.size_subtomo], yRange=[0, self.size_subtomo], padding=0)
            self.vBoxList[-1].setAspectLocked(True)

            self.index = [self.index[0] + 1, 0]
        else:

            self.iItemList[self.index[0]].setImage(image=subtomo)
            self.assigned[self.index[0]] = 1
            self.num_assigned += 1
            self.coordinates[self.index[0]] = self.position + [self.index[0]]
            self.index = [self.num_assigned, 0]
            for n, a in enumerate(self.assigned):
                if a == 0:
                    self.index = [n,1]
                    break

    def delete_subplot(self, position):
        for x,y,z,index in self.coordinates:
            if x == position[0] and y == position[1] and z == position[2]:
                blank = zeros((int(self.size_subtomo), int(self.size_subtomo)), dtype=float)
                self.iItemList[index].setImage(image=blank)
                if index < self.index[0]:
                    self.index = [index,1]
                self.assigned[index] = 0
                self.num_assigned -= 1
                break

    def newViewBox(self):
        row = int(self.index[0] / self.num_subtomo_per_row)
        col = int(self.index[0] % self.num_subtomo_per_row)
        if row*self.size_subplot > self.height:
            self.scrollarea.setGeometry(self.dx,self.dy,self.width,row*self.size_subplot)
        self.vBoxList.append( self.canvas.addViewBox(row=row, col=col) )
        self.vBoxList[-1].setMenuEnabled(False)
        self.vBoxList[-1].setMouseEnabled(False, False)
        if self.index[0] == 0: self.vBoxList[-1].scene().sigMouseClicked.connect(self.mouseClicked)

        self.vBoxList[-1].setGeometry(0, 0, self.size_subplot, self.size_subplot)

        #
        self.assigned.append(1)
        self.coordinates.append(self.position+[self.index[0]])
        self.num_assigned += 1

    def mouseClicked(self,event):
        ID = -2
        try:
            event.pos()

            for id, vb in enumerate(self.vBoxList):
                P = vb.mapSceneToView(event.scenePos())
                if P.x() > -0.0001 and P.y() > 0.0001 and P.x() < self.size_subtomo and P.y() < self.size_subtomo:
                    ID = id
                    break
        except:
            return

        if ID > -1:
            self.parent().slice = self.coordinates[ID][2]
            self.parent().replot()
            self.parent().update_circles()

    def reset_display_subtomograms(self, particleList, volume ):
        for child in self.vBoxList:
            self.canvas.removeItem(child)

        self.init_variables()

        self.num_subtomo_per_row = int(self.width/self.size_subplot)

        for n, (x,y,z) in enumerate(particleList):
            self.add_subplot( volume, [x,y,z] )
        self.show()

if __name__ == "__main__":
    import sys
    import numpy
    import glob
    import os

    headers = ["mdoc file", "create", "input files", 'redo', 'delete']
    types = ['txt', 'checkbox', 'combobox', 'checkbox', 'checkbox']
    sizes = [QHeaderView.ResizeToContents, 20, QHeaderView.ResizeToContents,]


    tomodir = 'Juliette/03_Tomographic_Reconstruction'
    rawdir = 'Juliette/01_Raw_Nanographs'

    processed = sorted(glob.glob('{}/tomogram_*/sorted/*.mdoc'.format(tomodir)))
    unprocessed = sorted(glob.glob('{}/*.mdoc'.format(rawdir)))
    values = []

    for t in processed:
        values.append([t, False, [1, 2, 3], True, True])
    for t in unprocessed:
        values.append([t, True, [1, 2, 3], False, False])

    # values = [[1,1,[1,3]],[1,1,[1,2]],[1,1,[1,2,3]]]*10

    app = QApplication(sys.argv)
    mw = SimpleTable(headers, types, values)
    mw.show()
    sys.exit(app.exec_())
