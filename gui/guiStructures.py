import os
import copy
import pickle

import pyqtgraph as pg

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5 import QtWidgets, QtCore, QtGui

#from pytom.basic.functions import initSphere
from pytom.basic.files import pdb2em
from pytom.gui.guiStyleSheets import *
from pytom.gui.mrcOperations import *
from pytom.gui.guiFunctions import initSphere

import pytom.gui.guiFunctions as guiFunctions
import traceback
from numpy import zeros, meshgrid, arange, sqrt
import numpy as np

from scipy.ndimage.filters import gaussian_filter
from ftplib import FTP_TLS, FTP
import lxml.etree as et

class BrowseWindowRemote(QMainWindow):
    '''This class creates a new windows for browsing'''
    def __init__(self, parent=None, initdir='/',filter=[''],search='file',credentials=['','',''],outputline='',
                 validate=True):
        super(BrowseWindowRemote, self).__init__(parent)
        self.setGeometry(50, 50, 800, 300)
        self.pathdisplay = QLineEdit()
        self.pathdisplay.setStyleSheet(WHITE)
        self.pathdisplay.setEnabled(False)
        self.splitter0 = splitter0 = QSplitter(Qt.Vertical)
        self.topleft = QListWidget()
        self.topleft.setSelectionMode(QtGui.QAbstractItemView.ExtendedSelection)
        self.topleft.setStyleSheet(WHITE)
        self.topleft.itemDoubleClicked.connect(self.repopulate_folder_list)
        self.topright = QListWidget()
        self.topright.setStyleSheet(WHITE)
        if search == 'file': self.topright.itemDoubleClicked.connect(self.select_file)
        self.splitter1 = splitter1 = QSplitter(Qt.Horizontal)
        splitter1.addWidget(self.topleft)
        splitter1.addWidget(self.topright)
        splitter1.setSizes([300,500])

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
                 validate=False, run_upon_complete=print, title=''):
        self.finished = False
        super(SelectFiles,self).__init__(parent,initdir=initdir,filter=[''],search='file',credentials=['','',''],
                                         outputline=outputline, validate=validate)
        self.setWindowTitle(title)
        self.run_upon_complete = run_upon_complete
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
        self.run_upon_complete()

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
    def __init__(self, fn=None, args=[], sig=True):
        #super(ProcessRunnable,self).__init__(parent)
        super(Worker, self).__init__()
        self.signals = WorkerSignals()
        self.fn = fn
        self.args = tuple(list(args))

        if sig: self.args = tuple(list(args)+[self.signals])


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
        widget.setStyleSheet("QComboBox{background:white; selection-background-color: #1989ac;}")
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
                       maximum=0, minimum=0, stepsize=0, widgetType=QSpinBox, decimals=0):

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
        if decimals:
            widget.setDecimals(decimals)
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


        if logvar:
            try:
                widget.setPlainText(self.logbook[wname])
            except:
                self.logbook[wname] = ''

            widget.setObjectName(wname)
            #widget.textChanged.connect(lambda ignore, name=wname, text='textfield':
            #                           self.update_logbook(name=name, text=text))

        parent.addWidget(widget, self.row, self.column, rowspan, columnspan)
        self.widgets[wname] = widget
        self.items[self.row][self.column] = widget
        self.row += rstep
        self.column += cstep

    def insert_label_spinbox(self, parent, wname, text='', rowspan=1, columnspan=1, rstep=1, cstep=-1, width=0, height=0,
                             validator=None, tooltip='', value=None, att=False, enabled=True, maximum=0, minimum=0,
                             wtype=QSpinBox, stepsize=0, logvar=True, decimals=0):

        self.insert_label(parent, text=text, cstep=1, alignment=QtCore.Qt.AlignRight, tooltip=tooltip)
        self.insert_spinbox(parent, wname, validator=validator, width=width, enabled=enabled,
                            maximum=maximum, minimum=minimum, cstep=cstep, rstep=rstep, value=value,
                            widgetType=wtype, stepsize=stepsize, logvar=logvar, decimals=decimals)

    def insert_label_line_push(self, parent, textlabel, wname, tooltip='', text='', cstep=-2, rstep=1, validator=None,
                               mode='folder', remote=False, pushtext='Browse', width=150, filetype='',action='',
                               initdir='', enabled=False):

        if action == '':
            action = self.browse


        self.insert_label(parent, text=textlabel, cstep=1, alignment=QtCore.Qt.AlignRight, tooltip=tooltip)
        self.insert_lineedit(parent, wname, cstep=1, logvar=True, text='',validator=validator,width=width,enabled=enabled)

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

    def insert_label_checkbox(self, parent, wname, textlabel, tooltip='', cstep=-1, rstep=1, logvar=True):
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
                            paramsSbatch={}, xmlfilename='', exefilename='exe.sh', jobfield=False, id=''):

        self.insert_label_action_label(parent, 'Generate command', cstep=1, rstep=-1, sizepolicy=self.sizePolicyB,
                                       action=self.gen_action,
                                       params=[[mode + 'XMLText']+paramsXML,
                                               [mode+'CommandText'] + paramsCmd,
                                               paramsSbatch])


        self.insert_checkbox(parent,mode + 'queue',text='queue',cstep=-3,rstep=1,logvar=True,alignment=Qt.AlignLeft)

        if jobfield:
            self.insert_textfield(parent, mode + 'XMLText', columnspan=3, rstep=1, cstep=2, width=600,logvar=False)
            self.insert_label(parent, alignment=Qt.AlignRight, rstep=1, cstep=-2, sizepolicy=self.sizePolicyB)
        self.insert_textfield(parent, mode + 'CommandText', columnspan=3, rstep=1, cstep=2, width=600, logvar=False)
        self.insert_label_action_label(parent, 'Execute command', rstep=1, action=self.exe_action,
                                       params=[exefilename, mode+'CommandText', xmlfilename, mode+'XMLText', action,
                                               paramsAction])

    def exe_action(self, params):

        if params[4]: params[4](params[5])

        if 1:
            # Check if one needs to write pytom related XML file.
            if params[2]:
                if len(params[2][0]) > 1:
                    tempfilename = os.path.join(self.widgets[params[2][0]].text(), params[2][1])
                else:
                    tempfilename = params[2]
                jobfile = open(tempfilename,'w')

                jobfile.write(self.widgets[params[3]].toPlainText())
                jobfile.close()

            # Write executable file
            if len(params[0][0]) > 1:
                exefilename = os.path.join(self.widgets[params[0][0]].text(), params[0][1])
            else:
                exefilename = params[0]
            exefile = open(exefilename, 'w')

            exefile.write(self.widgets[params[1]].toPlainText())
            exefile.close()

            if len(self.widgets[params[1]].toPlainText().split('SBATCH') ) > 2:
                os.system('{} {}'.format(self.qcommand, exefilename))
            else:
                proc = Worker(fn=os.system, args=['sh {}'.format(exefilename)], sig=False)
                proc.start()
#                os.system('sh {}'.format(params[0]))
        else:
            print ('Please check your input parameters. They might be incomplete.')

    def gen_action(self, params):
        mode = params[1][0][:-len('CommandText')]
        for i in range(2):
            if params[i][0] in self.widgets.keys():
                text = params[i][-1]
                if params[i][1:-1]:
                    d = []
                    for a in params[i][1:-1]:
                        if type(a) == type([]):
                            d.append(os.path.join(self.widgets[a[0]].text(), a[1]))
                        elif a in self.widgets.keys():
                            datatype = type(self.widgets[a])
                            if datatype in (QSpinBox, QDoubleSpinBox, QLineEdit): d.append(self.widgets[a].text())
                            elif datatype == QComboBox: d.append(self.widgets[a].currentText())
                            elif datatype == QCheckBox: d.append(str(int(self.widgets[a].isChecked())))
                            else: print(a, datatype)
                        elif type(str(a)) == type(''):
                            d.append(a)
                        else: print(a)

                    text = text.format( d=d )
                if i==0: self.widgets[params[i][0]].setPlainText(text)
        # Check if user wants to submit to queue. If so, add queue header.
        if self.widgets[mode+'queue'].isChecked():
            d = params[2]
            folder = d['folder']
            num_jobs_per_node = d['num_jobs_per_node']
            time = d['time']
            partition = d['partition']
            suffix = d['suffix']
            num_nodes = d['num_nodes']
            id = d['id']
            print(self.qparams.keys())
            if id:
                partition, num_nodes, cores, time = self.qparams[id].values()

            if type(folder) == type([]):

                outfolder = ''

                for extra in folder:
                    try:
                        extra = self.widgets[extra].text()
                    except:
                        pass
                    outfolder = os.path.join( outfolder, extra)
                folder = outfolder

            if self.custom:
                text = self.genSettingsWidgets['v00_QParams_CustomHeaderTextField'].toPlainText() + '\n\n' + text
            else:
                text = guiFunctions.gen_queue_header(name=d['fname'], folder=folder, cmd=d['cmd'], modules=d['modules'],
                                                     qtype=self.qtype, partition=partition, time=time,suffix=suffix,
                                                     num_jobs_per_node=num_jobs_per_node, num_nodes=num_nodes) + text

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
            elif text == 'textfield':
                self.logbook[name] = self.widgets[name].toPlainText()

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

    def fill_tab(self, id, headers, types, values, sizes, tooltip=[],wname='v02_batch_aligntable_', connect=0, nn=False):
        try:
            self.tables[id].setParent(None)
            self.pbs[id].setParent(None)
            self.ends[id].setParent(None)
        except:

            pass

        self.tables[id] = SimpleTable(headers, types, values, sizes, tooltip=tooltip,connect=connect)
        self.widgets['{}{}'.format(wname,id)] = self.tables[id]

        if nn:
            num_nodes = QSpinBox()
            num_nodes.setValue(2)
            num_nodes.setRange(1,9)
            num_nodes.setPrefix('Num Nodes: ')
            
        
        try: 
            self.num_nodes[id] = num_nodes
        except:
            self.num_nodes ={}
            self.num_nodes[id] = 0

        self.pbs[id] = QPushButton('Run')
        self.pbs[id].setSizePolicy(self.sizePolicyC)
        self.ends[id] = QWidget()
        self.ends[id].setSizePolicy(self.sizePolicyA)

        for n, a in enumerate( (self.tables[id], self.num_nodes[id], self.pbs[id], self.ends[id]) ):
            if n==1 and nn == False: continue
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
        self.insert_label_line(parent,'Sixe X (px)',wname='size_template_x',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent,'Size Y (px)',wname='size_template_y',rstep=1,cstep=-1,
                               validator=QIntValidator())
        self.insert_label_line(parent,'Size Z (px)',wname='size_template_z',rstep=1,cstep=-1,
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
        sizeX = int(self.widgets['size_template_x'].text())
        sizeY = int(self.widgets['size_template_y'].text())
        sizeZ = int(self.widgets['size_template_z'].text())
        try:
            fname = os.path.join(self.parent().frmdir, 'FRM_mask.mrc')
        except:
            tomoname = os.path.basename(self.parent().widgets['v03_TemplateMatch_tomoFname'].text())
            filename, file_extension = os.path.splitext(tomoname)
            if not os.path.exists(os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename)):
                os.mkdir(os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename))
            fname = os.path.join(self.parent().templatematchfolder, 'cross_correlation', filename, 'TM_mask.mrc')

        maskfilename = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.mrc')[0])
        if maskfilename and not maskfilename.endswith('.mrc'): maskfilename += '.mrc'
        try:
            success = initSphere(sizeX, sizeY, sizeZ, radius=radius, smooth=smooth, filename=maskfilename)
            if success:
                self.parent().widgets[params[-1]].setText(maskfilename)
            else:
                self.popup_messagebox('Error', 'Mask Generation Failed',
                                      'Generation of the mask failed. Please select an existing mask, or generate a mask yourself.')

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
        pdb2em('{}/{}.pdb'.format(self.folder,pdbid), pixel_size, cube_size, chain=chain, fname=out_fname,
               densityNegative=True)

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

                    if values[v][i] == 'UNCHANGED':
                        widget.setEnabled(False)

                elif types[i] == 'spinbox':
                    widget = QSpinBox()
                    widget.setValue(values[v][i])

                elif types[i] == 'combobox':
                    print(values[v][i])
                    widget = QWidget()
                    cb = QComboBox()
                    l = QVBoxLayout(widget)
                    l.addWidget(cb)
                    cb.setContentsMargins(0, 0, 0, 0)
                    l.setContentsMargins(0, 0, 0, 0)
                    for t in values[v][i]:
                        cb.addItem(t.split('/')[-1])

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

                cb.setContentsMargins(10, 0, 10, 0)
                l.setContentsMargins(10, 0, 10, 0)

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
                    cb.addItem(value.split('/')[-1])
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

    def addTabs(self, headers, widget=QWidget, subheaders=[], offx=0,offy=0,dimx=900,dimy=721,soff=0, sizeX=900,sizeY=700):
        self.size_policies()
        self.scrollarea = QScrollArea(self)
        if soff: self.scrollarea.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        #self.scrollarea.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.scrollarea.setSizePolicy(self.sizePolicyA)
        self.scrollarea.setWidgetResizable(True)
        self.scrollarea.setContentsMargins(0,0,0,0)
        #self.scrollarea.setGeometry(offx, offy, dimx, dimy)

        self.scrollarea.setGeometry(0,0,sizeX,sizeY)

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


class CreateMaskTMOld(QMainWindow, CommonFunctions):
    def __init__(self,parent=None, fname=''):
        super(CreateMaskTM,self).__init__(parent)
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
        self.logbook = {}
        self.radius = 50
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0,0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'

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
        self.title = self.parent().parent().widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle( "Create Mask File for: {}".format(os.path.basename(self.title ) ))
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

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt,'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.',width=100)

        self.insert_label(prnt, text='Size Selection: ',cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ',cstep=-1,rstep=1, value=100, minimum=1, maximum=1000,width=100, stepsize =20)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt,text='Step Size',cstep=1,alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt,wname='step_size',cstep=-1, value=10, rstep=1,width=100,
                            minimum=1, maximum=int(self.vol.shape[0]/4))

        self.insert_label(prnt,text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt,'Number Particles','numSelected',validator=vInt,width=100)

        self.insert_checkbox(prnt, 'show_mask', cstep=1)
        self.insert_label(prnt, text='Show Mask', cstep=0, rstep=1, alignment=Qt.AlignLeft)
        self.widgets['show_mask'].stateChanged.connect(self.show_maskTM)
        #self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def show_maskTM(self):
        w = self.widgets['show_mask']
        width = self.widgets['show_mask'].text()

        if w.isChecked():
            Y,Z,X= numpy.meshgrid(numpy.arange(float(self.dim)),numpy.arange(float(self.dim)),numpy.arange(float(self.dim)))

            if len(self.particleList):
                try:
                    self.mask += 0
                except:
                    self.mask = numpy.zeros_like(self.vol)

                for x,y,z,r in self.particleList:

                    x = self.dim-x
                    y = self.dim-y
                    z = self.dim-z
                    X -= self.dim - x -0.5
                    Y -= self.dim - y -0.5
                    Z -= self.dim - z -0.5
                    R = numpy.sqrt(X**2+Y**2+Z**2)
                    self.mask[R<=r] = 1
                    X += self.dim - x -0.5
                    Y += self.dim - y -0.5
                    Z += self.dim - z -0.5

                self.vol = self.volg = self.vol*self.mask

            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

    def open_load(self, q):
        if q.text() == 'Save Project': self.save_particleList()
        elif q.text() == 'Load Project': self.load_particleList()

    def save_particleList(self):
        try: self.mask += 0
        except: 
            print('Please create a mask before saving.')
            return
        import mrcfile
        base, ext = os.path.splitext(self.title)
        fname = os.path.join(self.parent().parent().pickpartfolder, 'mask_' + os.path.basename(self.title.replace(ext, '')))
        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', fname, filter="MRC File (*.mrc);; EM File (*.em)")[0])
        if not fname.endswith('.mrc'): fname += '.mrc'
        if not self.mask.sum(): self.mask +=1
        mrcfile.new(fname, self.mask, overwrite=True)
        self.parent().widgets['filenameMask'].setText(fname)
    
    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "MRC File (*.mrc);; EM File (*.em)")[0])
        if not filename: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol  = read(self.title)
            self.mask = copy.deepcopy( vol2npy(vol) )
            self.mask = self.mask.T

        elif self.title.endswith('mrc'):
            self.mask = read_mrc(self.title)

        if not self.widgets['show_mask'].isChecked():
            self.widgets['show_mask'].setChecked(True)
        else:
            self.show_maskTM()

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        tree = et.parse(self.xmlfile)


        for particle in tree.xpath("Particle"):
            remove = True

            position = []

            for tag in ('X', 'Y', 'Z'):
                position.append(float(particle.xpath('PickPosition')[0].get(tag)))

            x, y, z = position

            for cx, cy, cz, score in self.particleList:
                if abs(x - cx) < 1 or abs(y - cy) < 1 or abs(z - cz) < 1:
                    remove = False
                    break
            if remove:
                self.remove_element(particle)

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
            except:
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

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
        
        self.radius = int(self.widgets['size_selection: '].text()) / 2
        

    def stateScoreChanged(self):
        pass

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

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
        remove = evt.button()==2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x,y,z,s) in enumerate(self.particleList):
            if 1:#sqrt( (x-pos.x())**2 + (y-pos.y())**2 + (z-self.slice)**2 ) < self.radius:
                add = False
                if remove:
                    self.remove_point(n-num_deleted_items,z)
                    self.subtomo_plots.delete_subplot([x,y,z])
                    num_deleted_items += 1

        if add and not remove:
            X, Y = pos.x(), pos.y()

            self.add_points(pos, X, Y, self.slice, self.radius, self.radius,add=add,score=self.radius)
            self.subtomo_plots.add_subplot(self.vol, self.particleList[-1])

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def remove_from_coords(self,coords):
        cx,cy,cz = coords[:3]
        for n, (x,y,z,s) in enumerate(self.particleList):
            if sqrt( (x-cx)**2 + (y-cy)**2 + (z-cz)**2 ) < self.radius:
                self.remove_point(n, z)
                self.subtomo_plots.delete_subplot([x, y, z])
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0.):
        print(cx,cy,cz,score)
        self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz)), score])
        if radius < 1: return
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

        self.volg = self.vol
        self.mask = numpy.zeros_like(self.vol)
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

        for n, (cx, cy, cz, cs) in enumerate(plist):
            
            radius = 0
            if abs(cz - self.slice) <= cs:
                radius = sqrt(cs ** 2 - (cz - self.slice)**2)

            if self.xmlfile and len(plist) > 100:
                add = False
            else: add = True
            self.add_points(self.pos, cx, cy, cz, cs, radius,add=add,score=cs)
        self.slice += update


class CreateMaskTM(QMainWindow, CommonFunctions):
    def __init__(self, parent=None, fname=''):
        super(CreateMaskTM, self).__init__(parent)
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
        self.logbook = {}
        self.radius = 50
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0, 0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(250, 750), border=True)
        self.leftimage = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(750, 750), border=True)
        self.centimage = w.addViewBox(row=0, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(750, 250), border=True)
        self.bottomimage = w2.addViewBox(row=0, col=0)
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]

        self.layout.addWidget(w1, 0, 0)
        self.layout.addWidget(w, 0, 1)
        self.layout.addWidget(w2, 1, 1)
        self.layout.addWidget(self.operationbox, 1, 0)
        self.title = self.parent().parent().widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle("Create Mask File for: {}".format(os.path.basename(self.title)))
        self.centcanvas.wheelEvent = self.wheelEvent
        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        # self.centcanvas.sigKeyPress.connect(self.keyPress)
        self.centcanvas.sigMouseReleased.connect(self.empty)
        self.load_image()
        self.leftimage.setXRange(0, self.vol.shape[0])

        self.add_controls(self.layout_operationbox)

        pg.QtGui.QApplication.processEvents()

    def wheelEvent(self, event):

        step = event.angleDelta().y() / 120
        increment = int(self.widgets['step_size'].text()) * step
        if self.slice + increment < self.vol.shape[0] and self.slice + increment > -1:
            self.update_circles(increment)
            self.replot()

    def keyPressEvent(self, evt):
        if Qt.Key_G == evt.key():
            w = self.widgets['apply_gaussian_filter']
            w.setChecked(w.isChecked() == False)

        if Qt.Key_Right == evt.key():
            if self.slice + int(self.widgets['step_size'].text()) < self.dim:
                update = int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if Qt.Key_Left == evt.key():
            if self.slice > int(self.widgets['step_size'].text()) - 1:
                update = -1 * int(self.widgets['step_size'].text())
                self.update_circles(update)
                self.replot()

        if evt.key() == Qt.Key_Escape:
            self.close()

    def add_controls(self, prnt):
        vDouble = QtGui.QDoubleValidator()
        vInt = QtGui.QIntValidator()
        self.widgets = {}
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, text='Size Selection: ', cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ', cstep=-1, rstep=1, value=100, minimum=1, maximum=1000,
                            width=100, stepsize=20)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt, text='Step Size', cstep=1, alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt, wname='step_size', cstep=-1, value=10, rstep=1, width=100,
                            minimum=1, maximum=int(self.vol.shape[0] / 4))

        self.insert_label(prnt, text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt, 'Number Particles', 'numSelected', validator=vInt, width=100)

        self.insert_checkbox(prnt, 'show_mask', cstep=1)
        self.insert_label(prnt, text='Show Mask', cstep=0, rstep=1, alignment=Qt.AlignLeft)
        self.widgets['show_mask'].stateChanged.connect(self.show_maskTM)
        # self.insert_lineedit(prnt, 'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.', width=100)

        self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def show_maskTM(self):
        w = self.widgets['show_mask']
        width = self.widgets['show_mask'].text()

        if w.isChecked():
            # Y,Z,X= numpy.meshgrid(numpy.arange(float(self.dim)),numpy.arange(float(self.dim)),numpy.arange(float(self.dim)))

            if len(self.particleList):
                try:
                    self.mask += 0
                except:
                    self.mask = numpy.zeros_like(self.vol)

                self.mask *= 0

                for x, y, z, r in self.particleList:
                    dim = r * 2. + 2.
                    Y, Z, X = numpy.meshgrid(numpy.arange(dim), numpy.arange(dim), numpy.arange(dim))

                    X -= dim / 2 - 0.5
                    Y -= dim / 2 - 0.5
                    Z -= dim / 2 - 0.5
                    R = numpy.sqrt(X ** 2 + Y ** 2 + Z ** 2)

                    startx, endx = int(max(0, x - dim // 2)), int(min(self.dim, x + dim // 2))
                    starty, endy = int(max(0, y - dim // 2)), int(min(self.dim, y + dim // 2))
                    startz, endz = int(max(0, z - dim // 2)), int(min(self.dim, z + dim // 2))

                    rsx, rex, rsy, rey, rsz, rez = 0, dim, 0, dim, 0, dim

                    if endx - startx < dim:
                        if x - dim // 2 < 0:
                            rsx = dim // 2 - x
                        else:
                            rex = dim // 2 - x + self.dim
                    if endy - starty < dim:
                        if y - dim // 2 < 0:
                            rsy = dim // 2 - y
                        else:
                            rey = dim // 2 - y + self.dim
                    if endz - startz < dim:
                        if z - dim // 2 < 0:
                            rsz = dim // 2 - z
                        else:
                            rez = dim // 2 - z + self.dim

                    rsx, rex, rsy, rey, rsz, rez = map(int, (rsx, rex, rsy, rey, rsz, rez))
                    R = R[rsz:rez, rsy:rey, rsx:rex]

                    self.mask[startz:endz, starty:endy, startx:endx][R <= r] = 1

                self.vol = self.volg = self.vol * (0.25 + self.mask)

            else:
                self.vol = self.volg
        else:
            self.vol = self.backup.copy()
        self.replot_all()

    def open_load(self, q):
        if q.text() == 'Save Project':
            self.save_particleList()
        elif q.text() == 'Load Project':
            self.load_particleList()

    def save_particleList(self):
        try:
            self.mask += 0
        except:
            print('Please create a mask before saving.')
            return
        import mrcfile
        base, ext = os.path.splitext(self.title)
        fname = os.path.join(self.parent().parent().pickpartfolder,
                             'mask_' + os.path.basename(self.title.replace(ext, '')))
        fname = str(
            QFileDialog.getSaveFileName(self, 'Save particle list.', fname, filter="MRC File (*.mrc);; EM File (*.em)")[
                0])
        if not fname.endswith('.mrc'): fname += '.mrc'
        if not self.mask.sum(): self.mask += 1
        mrcfile.new(fname, self.mask, overwrite=True)
        self.parent().widgets['filenameMask'].setText(fname)

    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "MRC File (*.mrc);; EM File (*.em)")[0])
        if not filename: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol = read(self.title)
            self.mask = copy.deepcopy(vol2npy(vol))
            self.mask = self.mask.T

        elif self.title.endswith('mrc'):
            self.mask = read_mrc(self.title)

        if not self.widgets['show_mask'].isChecked():
            self.widgets['show_mask'].setChecked(True)
        else:
            self.show_maskTM()

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        tree = et.parse(self.xmlfile)

        for particle in tree.xpath("Particle"):
            remove = True

            position = []

            for tag in ('X', 'Y', 'Z'):
                position.append(float(particle.xpath('PickPosition')[0].get(tag)))

            x, y, z = position

            for cx, cy, cz, score in self.particleList:
                if abs(x - cx) < 1 or abs(y - cy) < 1 or abs(z - cz) < 1:
                    remove = False
                    break
            if remove:
                self.remove_element(particle)

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
            except:
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

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

        #self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def sizeChanged(self):
        a = self.widgets['size_selection: '].text()
        if not a: return

        plist = copy.deepcopy(self.particleList)

        self.radius = int(self.widgets['size_selection: '].text()) / 2

    def stateScoreChanged(self):
        pass

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

    def replot(self):
        crop = self.vol[int(self.slice), :, :]
        self.img1m.setImage(image=crop.T)

        # self.centcanvas.removeItem(self.hist)
        # self.hist = pg.HistogramLUTItem()
        self.hist.setImageItem(self.img1m)
        self.hist.setLevels(numpy.median(crop) - crop.std() * 3, numpy.median(crop) + crop.std() * 3)
        # self.centcanvas.addItem(self.hist)
        del crop

    def mouseHasMoved(self, evt):

        try:
            pos = self.centimage.mapSceneToView(evt.scenePos())
        except:
            return

        add = True
        remove = evt.button() == 2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x, y, z, s) in enumerate(self.particleList):
            if sqrt((x - pos.x()) ** 2 + (y - pos.y()) ** 2 + (z - self.slice) ** 2) < self.radius:
                add = False
                if remove:
                    self.remove_point(n - num_deleted_items, z)

                    num_deleted_items += 1
            add=True

        if add and not remove:
            X, Y = pos.x(), pos.y()

            self.add_points(pos, X, Y, self.slice, self.radius, self.radius, add=add, score=self.radius)

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def remove_from_coords(self, coords):
        cx, cy, cz = coords[:3]
        for n, (x, y, z, s) in enumerate(self.particleList):
            if sqrt((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2) < self.radius:
                self.remove_point(n, z)
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0., new=True):
        print(cx, cy, cz, cs, radius)
        if new:
            self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz)), score])
            print(cx, cy, cz, score)

        if radius < 1: return

        pos.setX(cx - radius)
        pos.setY(cy - radius)
        self.circles_cent.append(circle(pos, size=(radius) * 2))

        if abs(cz - self.slice) < self.radius:
            self.centimage.addItem(self.circles_cent[-1])

        if new:
            pos.setX(cx - self.radius)
            pos.setY(cz - self.radius)
            self.circles_bottom.append(circle(pos, size=self.radius * 2))
            self.bottomimage.addItem(self.circles_bottom[-1])

            pos.setX(cz - self.radius)
            pos.setY(cy - self.radius)
            self.circles_left.append(circle(pos, size=self.radius * 2))
            self.leftimage.addItem(self.circles_left[-1])

    def remove_point(self, n, z, from_particleList=True):
        if from_particleList:
            self.particleList.pop(n)
            # self.particleList = self.particleList[:n] + self.particleList[n + 1:]

        if abs(z - self.slice) <= self.radius:
            self.centimage.removeItem(self.circles_cent[n])
        self.circles_cent.pop(n)  # = self.circles_cent[:n] + self.circles_cent[n + 1:]

        self.leftimage.removeItem(self.circles_left[n])
        self.circles_left.pop(n)  # = self.circles_left[:n] + self.circles_left[n + 1:]

        self.bottomimage.removeItem(self.circles_bottom[n])
        self.circles_bottom.pop(n)  # = self.circles_bottom[:n] + self.circles_bottom[n + 1:]

    def remove_all(self):
        for image in self.image_list:
            for child in image.allChildren():
                if type(child) == MyCircleOverlay:
                    image.removeItem(child)

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.particleList = []

    def remove_cent(self):
        for child in self.centimage.allChildren():
            if type(child) == MyCircleOverlay:
                self.centimage.removeItem(child)

        for i in range(len(self.circles_cent)):
            self.circles_cent.pop(0)

    def load_image(self):
        if not self.title: return

        if self.title.endswith('em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol = read(self.title)
            self.vol = copy.deepcopy(vol2npy(vol))
            self.vol = self.vol.T
        elif self.title.endswith('mrc'):
            self.vol = read_mrc(self.title)

        # self.vol[self.vol < -4.] = -4.
        self.backup = self.vol.copy()
        self.origin = self.vol.copy()

        self.vol[self.vol < self.vol.min()] = self.vol.min()

        self.volg = self.vol
        self.mask = numpy.zeros_like(self.vol)
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
        self.slice += update
        self.remove_cent()

        for n, (cx, cy, cz, cs) in enumerate(self.particleList):
            radius = 0
            if abs(cz - self.slice) <= cs:
                radius = sqrt(cs ** 2 - (cz - self.slice) ** 2)

            add = False
            self.add_points(self.pos, cx, cy, cz, cs, radius, add=add, score=cs, new=False)

        print(len(self.circles_bottom), len(self.circles_cent), len(self.circles_left), len(self.particleList))


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
        self.setGeometry(0, 0, 800, 800)
        self.operationbox = QWidget()
        self.layout_operationbox = prnt = QGridLayout()
        self.operationbox.setLayout(self.layout_operationbox)
        self.add_toolbar(self.open_load)
        self.logbook = {}
        self.radius = 7
        self.jump = 1
        self.current_width = 0.
        self.pos = QPoint(0,0)
        self.max_score = 1.
        self.min_score = 0.
        self.xmlfile = ''
        self.filetype = 'txt'

        self.circles_left = []
        self.circles_cent = []
        self.circles_bottom = []
        self.circles_list = [self.circles_left, self.circles_cent, self.circles_bottom]
        self.particleList = []

        self.leftcanvas = w1 = pg.GraphicsWindow(size=(200, 600), border=True)
        self.leftimage  = w1.addViewBox(row=0, col=0)
        self.leftimage.setMouseEnabled(False, False)

        self.centcanvas = w = KeyPressGraphicsWindow(size=(600, 600), border=True)
        self.centimage  = w.addViewBox(row=0, col=0, lockAspect=True)
        self.centimage.setMenuEnabled(False)
        self.target = w3 = pg.ImageView()

        self.bottomcanvas = w2 = pg.GraphicsWindow(size=(600, 200), border=True)
        self.bottomimage  = w2.addViewBox(row=0, col=0 )
        self.bottomimage.setMouseEnabled(False, False)

        self.image_list = [self.leftimage, self.centimage, self.bottomimage]



        self.layout.addWidget(w1,0,0)
        self.layout.addWidget(w,0,1)
        self.layout.addWidget(w2,1,1)
        self.layout.addWidget(self.operationbox,1,0)
        self.title = parent.widgets['v03_manualpp_tomogramFname'].text()
        if not self.title: self.title = 'Dummy Data'
        self.setWindowTitle( "Manual Particle Selection From: {}".format( os.path.basename(self.title)) )
        self.centcanvas.wheelEvent = self.wheelEvent

        self.centimage.scene().sigMouseClicked.connect(self.mouseHasMoved)
        self.leftimage.scene().sigMouseClicked.connect(self.mouseHasMovedLeft)
        self.bottomimage.scene().sigMouseClicked.connect(self.mouseHasMovedBottom)


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

        self.insert_checkbox(prnt, 'apply_gaussian_filter', cstep=1)
        self.insert_label(prnt, text='Gaussian Filter', cstep=1, alignment=Qt.AlignLeft)
        self.insert_lineedit(prnt,'width_gaussian_filter', validator=vDouble, rstep=1, cstep=-1, value='1.',width=100)

        self.insert_label(prnt, text='Size Selection: ',cstep=1)
        self.insert_spinbox(prnt, wname='size_selection: ',cstep=-1,rstep=1, value=14, minimum=1, maximum=1000,width=100)
        self.widgets['size_selection: '].valueChanged.connect(self.sizeChanged)

        self.insert_label(prnt,text='Step Size',cstep=1,alignment=Qt.AlignLeft)
        self.insert_spinbox(prnt,wname='step_size',cstep=-1, value=1, rstep=1,width=100,
                            minimum=1, maximum=int(self.vol.shape[0]/4))

        self.insert_label(prnt,text='', cstep=0)
        self.widgets['apply_gaussian_filter'].stateChanged.connect(self.stateGaussianChanged)

        self.insert_label_line(prnt,'Number Particles','numSelected',validator=vInt,width=100)
        self.insert_label_spinbox(prnt,'minScore','Minimal Score', value=0.,wtype=QDoubleSpinBox, minimum=0, maximum=1,
                                  stepsize=0.05, width=100, decimals=4)
        self.insert_label_spinbox(prnt, 'maxScore', 'Maximal Score', value=1., wtype=QDoubleSpinBox, minimum=0,
                                  maximum=1., stepsize=0.05, width=100, decimals=4, cstep=-2)

        self.widgets['minScore'].valueChanged.connect(self.stateScoreChanged)
        self.widgets['maxScore'].valueChanged.connect(self.stateScoreChanged)
        
        self.insert_checkbox(prnt, 'apply_mask', cstep=1)
        self.insert_label(prnt, text='Apply mask', cstep=1, alignment=Qt.AlignLeft)
        self.insert_pushbutton(prnt,'Create', rstep=1, cstep=-1, action=self.gen_mask,params=['filenameMask'])

        self.insert_lineedit(prnt,'filenameMask',  cstep=1, value='')
        params = ['file', self.widgets['filenameMask'], ['mrc', 'em'], False, self.parent().pickpartfolder]
        self.insert_pushbutton(prnt,'Browse', action=self.browse, params=params, width=100)
        self.widgets['filenameMask'].textChanged.connect(self.updateMask)




        #self.insert_label(prnt, sizepolicy=self.sizePolicyA)

    def updateMask(self):
        fname = self.widgets['filenameMask'].text()
        if fname.endswith('.em'):
            from pytom.basic.files import read
            from pytom_numpy import vol2npy
            vol  = read(fname)
            self.mask = copy.deepcopy( vol2npy(vol) )
            self.mask = self.mask.T
        elif self.title.endswith('mrc'):
            self.mask = read_mrc(fname)

        self.widgets['apply_mask'].setChecked(True)
        if self.filetype == 'xml': 
            self.load_xmlFile(self.FNAME)

    def gen_mask(self, params):
        self.genmask = CreateMaskTM(self)
        self.genmask.show()
        #self.widgets[params[0]].setText(fname)

    def open_load(self, q):
        if q.text() == 'Save Project': self.save_particleList()
        elif q.text() == 'Load Project': self.load_particleList()

    def save_particleList(self):
        if self.filetype == 'xml':
            self.save_particleList_xml()
        else:
            self.save_particleList_txt()

    def save_particleList_xml(self):
        if 1:
            self.remove_deselected_particles_from_XML()
        else:
            print('particleList not saved.')

    def save_particleList_txt(self):
        ext = '.'+os.path.basename(self.title).split('.')[-1]
        fname = os.path.join(self.parent().pickpartfolder,'coords_' + os.path.basename(self.title.replace(ext,'')))

        fname = str(QFileDialog.getSaveFileName( self, 'Save particle list.', fname, filter='*.txt')[0])
        if not fname.endswith('.txt'): fname += '.txt'

        if fname:

            folder = os.path.dirname(self.title)
            if 'INFR' in os.path.basename(self.title):
                key = 'INFR'
            else:
                key ='WBP'

            inputJobName = "cat {}/{}_reconstruction.sh | grep '{}' | {}"

            headerInfo = []
            for i in ('referenceMarkerIndex', 'projectionBinning'):
                if 1:
                    print(inputJobName.format(folder,key, i, "awk '{print $2}'"))
                    d = os.popen(inputJobName.format(folder,key, i, "awk '{print $2}'")).read()[:-1]
                else:
                    d = 1
                headerInfo.append(d)


            out = open(fname, 'w')
            out.write('#\n')
            out.write('#TOMONAME \t{}\n'.format(os.path.basename(self.title)))
            out.write('#MARKERINDEX \t{}\n'.format(headerInfo[0]))
            out.write('#BINNINGFACTOR \t{}\n'.format(headerInfo[1]))
            out.write('#\n')

            for x, y, z, s in self.particleList:
                outtext =  '{:8.0f} {:8.0f} {:8.0f}\n'.format(x,y,z)
                out.write(outtext)

    def load_particleList(self):
        filetype = 'txt'
        initdir = self.parent().pickpartfolder

        filename = str(QFileDialog.getOpenFileName(self, 'Open file', initdir, "Coordinate files (*.txt);; Particle List (*.xml)")[0])
        if not filename: return

        if filename.endswith('.txt'):
            self.filetype = 'txt'
            self.load_txtFile(filename)

        elif filename.endswith('.xml'):
            self.filetype = 'xml'
            self.load_xmlFile(filename)
            self.xmlfile = filename

    def remove_element(self, el):
        parent = el.getparent()
        if el.tail.strip():
            prev = el.getprevious()
            if prev:
                prev.tail = (prev.tail or '') + el.tail
            else:
                parent.text = (parent.text or '') + el.tail
        parent.remove(el)

    def remove_deselected_particles_from_XML(self):
        tree = et.parse(self.xmlfile)

        print(self.max_score, self.min_score)
        for particle in tree.xpath("Particle"):
            remove = True

            position = []

            for tag in ('X', 'Y', 'Z'):
                position.append(float(particle.xpath('PickPosition')[0].get(tag)))

            score = float(particle.xpath('Score')[0].get('Value'))
            x, y, z = position
            print(x, y, z, score)
            for cx, cy, cz, s in self.particleList:
                if abs(x - cx) < 1 and abs(y - cy) < 1 and abs(z - cz) < 1 and score <= self.max_score and score >= self.min_score:
                    remove = False
                    break
            if remove:
                self.remove_element(particle)

        fname = str(QFileDialog.getSaveFileName(self, 'Save particle list.', self.xmlfile, filter='*.xml')[0])
        if not fname.endswith('.xml'): fname += '.xml'
        if fname:
            try:
                tree.write(fname.replace('.xml', '_deselected.xml'), pretty_print=True)
            except:
                print('writing {} failed.'.format(fname))
        else:
            print('No file written.')

    def load_xmlFile(self,filename):
        from lxml import etree
        self.FNAME = filename
        xmlObj = etree.parse(filename)
        particles = xmlObj.xpath('Particle')

        self.remove_all()
        self.slice = int(self.vol.shape[0] / 2)

        for p in particles:
            try:    score = float( p.xpath('Score')[0].get('Value') )
            except: score = 0.5
            

            if score < self.max_score and score > self.min_score:
                
                x,y,z = map(float, (p.xpath('PickPosition')[0].get('X'), p.xpath('PickPosition')[0].get('Y'), p.xpath('PickPosition')[0].get('Z')))
                dx,dy,dz = self.vol.shape

                if self.mask[int(x),int(y),int(z)] or not self.widgets['apply_mask'].isChecked(): include =True
                else: include = False
                if not include: continue

                if abs(dx/2-x) > dx/2 or abs(dx/2-x) > dx/2 or abs(dx/2-x) > dx/2 :
                    print('particle not added: ', x,y,z)
                    continue
                self.add_points(self.pos, int(round(x)), int(round(y)), int(round(z)), score, self.radius, score=score)


        self.replot()
        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

        particlePos = [line.split() for line in open(filename).readlines() if 'PickPosition' in line]

    def load_txtFile(self, filename):
        particlePos = [map(float, line.split()) for line in open(filename).readlines() if line != '' and not '#' in line]

        self.remove_all()
        self.slice = int(self.vol.shape[0]/2)

        for x,y,z in particlePos:
            self.add_points(self.pos, int(round(x)), int(round(y)), int(round(z)), 0., self.radius, add=True)

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


        for nn, (x,y,z,s) in enumerate(plist):
            add= True
            if self.xmlfile and len(self.particleList) > 100: add=False
            self.add_points(QPoint(0,0), x, y, z, s, self.radius, add=add)

        self.subtomo_plots.size_subtomo = self.radius*2

        self.subtomo_plots.reset_display_subtomograms(self.particleList, self.vol)

    def stateScoreChanged(self):
        self.min_score = float(self.widgets['minScore'].value())
        self.max_score = float(self.widgets['maxScore'].value())
        if self.xmlfile and os.path.exists(self.xmlfile):
            self.load_xmlFile(self.xmlfile)

    def replot_all(self):
        self.replot()
        self.img1a.setImage(image=self.vol.sum(axis=2))
        self.img1b.setImage(image=self.vol.sum(axis=1).T)

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
        remove = evt.button()==2
        self.pos = pos

        if pos.x() < 0 or pos.y() < 0 or pos.x() >= self.vol.shape[1] or pos.y() >= self.vol.shape[2]:
            return

        num_deleted_items = 0
        for n, (x,y,z,s) in enumerate(self.particleList):
            if sqrt( (x-pos.x())**2 + (y-pos.y())**2 + (z-self.slice)**2 ) < self.radius:
                add = False
                if remove:
                    self.remove_point(n-num_deleted_items,z)
                    self.subtomo_plots.delete_subplot([x,y,z])
                    num_deleted_items += 1

        if add and not remove:
            X, Y = pos.x(), pos.y()

            self.add_points(pos, X, Y, self.slice, 0, self.radius,add=add)
            self.subtomo_plots.add_subplot(self.vol, self.particleList[-1])

        self.widgets['numSelected'].setText(str(len(self.particleList)))

    def mouseHasMovedBottom(self, evt):
        pos = self.bottomimage.mapSceneToView( evt.scenePos() )
        if pos.y() < 0 or pos.y() >= self.vol.shape[2]: return
        step = pos.y() - self.slice
        self.update_circles(step)
        self.replot()
        print('bottom: ', pos.x(), pos.y())

    def mouseHasMovedLeft(self, evt):
        pos = self.leftimage.mapSceneToView( evt.scenePos() )
        if pos.x() < 0 or pos.x() >= self.vol.shape[1]: return

        step = pos.x()-self.slice
        self.update_circles(step)
        self.replot()
        print('left: ', pos.x(), pos.y())

    def remove_from_coords(self,coords):
        cx,cy,cz = coords[:3]
        for n, (x,y,z,s) in enumerate(self.particleList):
            if sqrt( (x-cx)**2 + (y-cy)**2 + (z-cz)**2 ) < self.radius:
                self.remove_point(n, z)
                self.subtomo_plots.delete_subplot([x, y, z])
                break

    def add_points(self, pos, cx, cy, cz, cs, radius, add=False, score=0.):
        self.particleList.append([int(round(cx)), int(round(cy)), int(round(cz)), score])

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

        self.mask = numpy.ones_like(self.vol)
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

        for n, (cx, cy, cz, cs) in enumerate(plist):
            radius = self.radius
            if abs(cz - self.slice) < self.radius:
                radius = sqrt(self.radius ** 2 - (cz - self.slice)**2)
            if self.xmlfile and len(plist) > 100:
                add = False
            else: add = True
            self.add_points(self.pos, cx, cy, cz, cs, radius,add=add)
        self.slice += update

class QParams():
    def __init__(self, time=12, queue='defq', nodes=1, cores=20):
        self.time = time
        self.queue = queue
        self.nodes = nodes
        self.cores = cores

    def update(self, mode, parent):
        self.queue = parent.widgets[mode + 'queueName'].text()
        self.time  = parent.widgets[mode + 'maxTime'].value()
        self.nodes = parent.widgets[mode + 'numberOfNodes'].value()
        self.cores = parent.widgets[mode + 'numberOfCores'].value()
        with open(os.path.join(parent.projectname, 'qparams.pickle'), 'wb') as handle:
            pickle.dump(parent.qparams, handle, protocol=pickle.HIGHEST_PROTOCOL)

        parent.parent().qparams = parent.qparams
        if 1:
            for tab in (parent.parent().CD, parent.parent().TR, parent.parent().PP, parent.parent().SA):
                tab.qparams = parent.qparams
        else:
            pass

    def values(self):
        return [self.queue, self.nodes, self.cores, self.time]

class GeneralSettings(QMainWindow, GuiTabWidget, CommonFunctions):
    def __init__(self,parent):
        super(GeneralSettings, self).__init__(parent)
        self.stage='generalSettings_'
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname
        self.logbook = self.parent().logbook
        self.setGeometry(0, 0, 600, 500)
        self.qcommanddict = {'slurm': 'sbatch', 'sge': 'qsub', 'torque': 'qsub', 'none': 'none'}\

        headers = ['Queuing Parameters', 'Extra']
        subheaders = [[], ] * len(headers)

        self.addTabs(headers=headers, widget=GuiTabWidget, subheaders=subheaders, sizeX=600, sizeY=500)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}
        self.widgets = {}
        self.subprocesses = 10

        self.tabs = {'tab1': self.tab1,
                     'tab2': self.tab2,
                     }

        self.tab_actions = {'tab1': self.tab1UI,
                            'tab2': self.tab2UI,
                            }

        for i in range(len(headers)):
            t = 'tab{}'.format(i + 1)
            empty = 1 * (len(subheaders[i]) == 0)
            for j in range(len(subheaders[i]) + empty):
                tt = t + str(j + 1) * (1 - empty)
                if tt in ('tab1', 'tab2'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()

                if tt in ('tab1', 'tab2'):
                    self.tab_actions[tt]()

                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])

    def setQNames(self):
        self.qnames = ['defq', 'fastq']

    def updateJobName(self, mode='v00_QParams_'):
        jobname = self.widgets[mode + 'jobName'].currentText()
        self.currentJobName = jobname
        self.widgets[mode + 'queueName'].setText(self.qparams[self.currentJobName].queue)
        self.widgets[mode + 'maxTime'].setValue(self.qparams[self.currentJobName].time)
        self.widgets[mode + 'numberOfNodes'].setValue(self.qparams[self.currentJobName].nodes)
        self.widgets[mode + 'numberOfCores'].setValue(self.qparams[self.currentJobName].cores)

    def tab1UI(self):
        self.jobnames = ['SingleAlignment', 'BatchAlignment',
                         'ReconstructWBP', 'ReconstructINFR', 'BatchReconstruct',
                         'CTFDetermination', 'SingleCTFCorrection', 'BatchCTFCorrection',
                         'SingleTemplateMatch','SingleExtractCandidates','BatchTemplateMatch',
                         'SingleSubtomoReconstruct', 'BatchSubtomoReconstruct',
                         'FRMAlignment','GLocalAlignment',
                         'PairwiseCrossCorrelation', 'CPCA', 'AutoFocusClassification']
        self.setQNames()
        self.currentJobName = self.jobnames[0]

        if os.path.exists(os.path.join(self.projectname, 'qparams.pickle')):
            with open(os.path.join(self.projectname, 'qparams.pickle'), 'rb') as handle:
                self.qparams = pickle.load(handle)
        else:
            self.qparams = {}

        for jobname in self.jobnames:
            if not jobname in self.qparams.keys():
                try:
                    self.qparams[jobname] = QParams(queue=self.qnames[0])
                except:
                    self.qparams[jobname] = QParams()

        id = 'tab1'
        self.row, self.column = 0, 1
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_QParams_'

        w = 150
        last, reftilt = 10, 5
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_label_combobox(parent, 'Queuing System ', mode + 'qType', ['Slurm','Torque','SGE', 'None'],
                                    tooltip='Select a particleList which you want to plot.\n', logvar=True)
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_label_combobox(parent, 'Job Submission Parameters', mode + 'jobName', self.jobnames, logvar=True,
                                   tooltip='Select the jon for which you want to adjust the queuing parameters.\n')
        self.insert_label_line(parent, 'Submit to: ', mode+'queueName', value=self.qparams[self.currentJobName].queue)
        self.insert_label_spinbox(parent,mode+'numberOfNodes', 'Number of Nodes',
                                  value=self.qparams[self.currentJobName].nodes, minimum=1, stepsize=1,
                                  wtype=QSpinBox)
        self.insert_label_spinbox(parent,mode+'numberOfCores', 'Number of Cores',
                                  value=self.qparams[self.currentJobName].cores, minimum=1, maximum=10000,
                                  stepsize=1, wtype=QSpinBox)
        self.insert_label_spinbox(parent,mode+'maxTime', 'Maximum Time (hours)',
                                  value=self.qparams[self.currentJobName].time, minimum=1, stepsize=1,
                                  wtype=QSpinBox)
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w, columnspan=2)
        self.insert_checkbox(parent, mode + 'CustomHeader', 'use custom header for queue', cstep=0, rstep=1,
                            alignment=Qt.AlignLeft, logvar=True, columnspan=2)
        self.insert_textfield(parent,mode+'CustomHeaderTextField', columnspan=5,rstep=1,cstep=0, logvar=True)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

        self.widgets[mode + 'qType'].currentTextChanged.connect(lambda d, m=mode: self.updateQType(m))
        self.widgets[mode + 'CustomHeader'].stateChanged.connect(lambda d, m=mode: self.updateCustomHeader(m))
        self.widgets[mode + 'jobName'].currentTextChanged.connect(lambda d, m=mode: self.updateJobName(m))
        self.widgets[mode + 'queueName'].textChanged.connect(lambda d, m=mode:
                                                             self.qparams[self.currentJobName].update(mode,self))
        self.widgets[mode + 'numberOfNodes'].valueChanged.connect(lambda d, m=mode:
                                                                  self.qparams[self.currentJobName].update(m,self))
        self.widgets[mode + 'numberOfCores'].valueChanged.connect(lambda d, m=mode:
                                                                  self.qparams[self.currentJobName].update(m,self))
        self.widgets[mode + 'maxTime'].valueChanged.connect(lambda d, m=mode:
                                                            self.qparams[self.currentJobName].update(m,self))

        self.updateCustomHeader(mode)
        self.qparams[self.currentJobName].update(mode, self)

        for n, value in enumerate(self.qcommanddict.values()):
            if value != 'none':
                if os.path.exists( os.popen('which {}'.format(value)).read()[:-1] ):
                    self.widgets[mode + 'qType'].setCurrentIndex(n)
                    break
            else: self.widgets[mode + 'qType'].setCurrentIndex(n)

    def updateCustomHeader(self, mode):
        try:
            for tab in (self.parent().CD, self.parent().TR, self.parent().PP, self.parent().SA):
                tab.custom = self.widgets[mode + 'CustomHeader'].isChecked()
                tab.genSettingsWidgets = self.widgets
        except:
            pass

    def updateQType(self, mode):

        qtype = self.widgets[mode + 'qType'].currentText().lower()
        qcommand =  self.qcommanddict[qtype]

        self.parent().qtype = qtype
        self.parent().qcommand = qcommand
        active = (qtype != 'none')
        try:
            for tab in (self.parent().CD, self.parent().TR, self.parent().PP, self.parent().SA):
                tab.qtype = qtype
                tab.qcommand = qcommand
                for key in tab.widgets.keys():
                    if key.endswith('queue'):
                        if not active:
                            tab.widgets[key].setChecked(False)
                        tab.widgets[key].setEnabled(active)

        except:
            pass

    def showTMPlot(self, mode):
        from pytom.plotting.plottingFunctions import plotTMResults

        normal = self.widgets[mode + 'particleListNormal'].text()
        mirrored = self.widgets[mode + 'particleListMirrored'].text()

        plotTMResults([normal, mirrored], labels=['Normal', 'Mirrored'])

    def tab2UI(self):
        id = 'tab2'
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_PlotFSC'
        w = 150
        last, reftilt = 10, 5
        self.insert_label(parent, rstep=1, cstep=0, sizepolicy=self.sizePolicyB, width=w)
        self.insert_label_line_push(parent, 'FSC File (ascii)', mode + 'FSCFilename', mode='file', width=w,
                                    initdir=self.projectname,
                                    filetype='dat', tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_spinbox(parent, mode + 'BoxSize', text='Dimension of Image', tooltip='Box size of 3D object',
                                  value=64, minimum=1, stepsize=1, width=w)
        self.insert_label_spinbox(parent, mode + 'PixelSize', text='Pixel Size',
                                  tooltip='Pixel size of a voxel in teh object.',
                                  value=2.62, minimum=1, stepsize=1, wtype=QDoubleSpinBox, decimals=2, width=w)
        self.insert_label_spinbox(parent, mode + 'CutOff', text='Resolution Cutoff', value=0, minimum=0, stepsize=0.1,
                                  wtype=QDoubleSpinBox, decimals=3, width=w, cstep=0,
                                  tooltip='Cut-off used to determine the resolution of your object from the FSC curve. \nTypical values are 0.5 or 0.143')

        self.insert_pushbutton(parent, 'Plot!', action=self.showFSCPlot, params=mode, rstep=1, cstep=0)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

    def showFSCPlot(self, mode):
        from pytom.bin.plotFSC import plot_FSC
        filename = self.widgets[mode + 'FSCFilename'].text()
        pixel_size = self.widgets[mode + 'PixelSize'].value()
        box_size = self.widgets[mode + 'BoxSize'].value()
        cut_off = self.widgets[mode + 'CutOff'].value()
        show_image = True
        outFname = 'temp.png'
        if filename and outFname:
            plot_FSC(filename, pixel_size, boxsize=box_size, show_image=show_image, c=cut_off)


class PlotWindow(QMainWindow, GuiTabWidget, CommonFunctions):
    def __init__(self,parent):
        super(PlotWindow, self).__init__(parent)
        self.stage='generalSettings_'
        self.pytompath = self.parent().pytompath
        self.projectname = self.parent().projectname

        self.setGeometry(0,0,500,300)

        headers = ['Template Matching Results', 'FSC Curve']
        subheaders  = [[],]*len(headers)

        self.addTabs(headers=headers,widget=GuiTabWidget, subheaders=subheaders,sizeX=500,sizeY=300)

        self.table_layouts = {}
        self.tables = {}
        self.pbs = {}
        self.ends = {}
        self.checkbox = {}
        self.num_nodes = {}
        self.widgets={}
        self.subprocesses = 10

        self.tabs = {'tab1': self.tab1,
                     'tab2':  self.tab2,
                     }

        self.tab_actions = {'tab1':  self.tab1UI,
                            'tab2':  self.tab2UI,
                            }

        for i in range(len(headers)):
            t = 'tab{}'.format(i+1)
            empty = 1*(len(subheaders[i]) == 0)
            for j in range(len(subheaders[i])+empty):
                tt = t+str(j+1)*(1-empty)
                if tt in ('tab1', 'tab2'):
                    self.table_layouts[tt] = QGridLayout()
                else:
                    self.table_layouts[tt] = QVBoxLayout()

                if tt in ('tab1','tab2'):
                    self.tab_actions[tt]()


                tab = self.tabs[tt]
                tab.setLayout(self.table_layouts[tt])

    def tab1UI(self):

        id = 'tab1'
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_PlotTM'
        w = 150
        last, reftilt = 10, 5
        self.insert_label(parent, cstep=0, rstep=1, sizepolicy=self.sizePolicyB, width=w)
        self.insert_label_line_push(parent, 'particleList ', mode + 'particleListNormal', width=w,mode='file',
                                    filetype='xml',enabled=True,
                                    tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_line_push(parent, 'particleList Mirrored', mode + 'particleListMirrored', width=w,
                                    mode='file',filetype='xml', enabled=True,
                                    tooltip='Select a particleList which you want to plot.\n', cstep=-1)
        self.insert_pushbutton(parent,'Plot',action=self.showTMPlot, params=mode,rstep=1,cstep=0)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)


    def showTMPlot(self, mode):
        from pytom.plotting.plottingFunctions import plotTMResults

        normal = self.widgets[mode+'particleListNormal'].text()
        mirrored = self.widgets[mode+'particleListMirrored'].text()

        plotTMResults([normal, mirrored], labels=['Normal', 'Mirrored'])

    def tab2UI(self):
        id = 'tab2'
        self.row, self.column = 0, 0
        rows, columns = 20, 20
        self.items = [['', ] * columns, ] * rows
        parent = self.table_layouts[id]
        mode = 'v00_PlotFSC'
        w=150
        last, reftilt = 10, 5
        self.insert_label(parent, rstep=1, cstep=0, sizepolicy=self.sizePolicyB, width=w)
        self.insert_label_line_push(parent, 'FSC File (ascii)', mode + 'FSCFilename',mode='file',width=w, initdir=self.projectname,
                                    filetype='dat', tooltip='Select a particleList which you want to plot.\n')
        self.insert_label_spinbox(parent, mode+'BoxSize', text='Dimension of Image', tooltip='Box size of 3D object',
                                  value=64,minimum=1,stepsize=1, width=w)
        self.insert_label_spinbox(parent, mode +'PixelSize', text='Pixel Size', tooltip='Pixel size of a voxel in teh object.',
                                  value=2.62, minimum=1, stepsize=1, wtype=QDoubleSpinBox, decimals=2, width=w)
        self.insert_label_spinbox(parent, mode +'CutOff', text='Resolution Cutoff', value=0, minimum=0, stepsize=0.1,
                                  wtype=QDoubleSpinBox, decimals=3, width=w,cstep=0,
                                  tooltip='Cut-off used to determine the resolution of your object from the FSC curve. \nTypical values are 0.5 or 0.143')

        self.insert_pushbutton(parent,'Plot!',action=self.showFSCPlot, params=mode, rstep=1,cstep=0)
        self.insert_label(parent, cstep=1, rstep=1, sizepolicy=self.sizePolicyA)

    def showFSCPlot(self, mode):
        from pytom.bin.plotFSC import plot_FSC
        filename = self.widgets[mode+'FSCFilename'].text()
        pixel_size = self.widgets[mode+'PixelSize'].value()
        box_size = self.widgets[mode+'BoxSize'].value()
        cut_off = self.widgets[mode+'CutOff'].value()
        show_image=True
        outFname = 'temp.png'
        if filename and outFname:
            plot_FSC(filename, pixel_size, boxsize=box_size, show_image=show_image, c=cut_off )


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
        [x,y,z,s]  = position
        xmin = max(0, x - int(self.size_subtomo / 2))
        xmax = min(tomo.shape[1], x + int(self.size_subtomo / 2))
        ymin = max(0, y - int(self.size_subtomo / 2))
        ymax = min(tomo.shape[2], y + int(self.size_subtomo / 2))

        subtomo = zeros((int(self.size_subtomo),int(self.size_subtomo)),dtype=float)
        subtomo[:xmax-xmin,:ymax-ymin] = (tomo[position[2]-4:position[2]+4].sum(axis=0).T)[xmin:xmax, ymin:ymax]

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
        for x,y,z,s, index in self.coordinates:
            if x == position[0] and y == position[1] and z == position[2]:
                blank = zeros((int(self.size_subtomo), int(self.size_subtomo)), dtype=float)
                self.iItemList[index].setImage(image=blank)
                if index < self.index[0]:
                    self.index = [index,1]
                self.assigned[index] = 0
                self.num_assigned -= 1
                self.coordinates[index][:3] = [-1, -1, -1]
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

        if ID < 0: return

        if event.button() == 2:
            self.parent().remove_from_coords(self.coordinates[ID])

        elif self.coordinates[ID][2] > -1:
            print(self.coordinates[ID])

            #self.parent().pos.setX(self.coordinates[ID][0] - self.parent().radius)
            #self.parent().pos.setY(self.coordinates[ID][1] - self.parent().radius)
            #self.parent().centimage.addItem(circle(self.parent().pos, size=(self.parent().radius)*1.5, color=Qt.yellow))
            self.setWindowTitle("Error Score: {:6.3f}".format( self.coordinates[ID][3]))
            self.parent().slice = self.coordinates[ID][2]
            self.parent().replot()
            self.parent().update_circles()

            self.parent().pos.setX(self.coordinates[ID][0] - self.parent().radius)
            self.parent().pos.setY(self.coordinates[ID][1] - self.parent().radius)
            self.parent().centimage.addItem(circle(self.parent().pos, size=(self.parent().radius) * 2, color=Qt.red))
    def reset_display_subtomograms(self, particleList, volume ):
        for child in self.vBoxList:
            self.canvas.removeItem(child)

        self.init_variables()

        self.num_subtomo_per_row = int(self.width/self.size_subplot)

        for n, (x,y,z,s) in enumerate(particleList):
            self.add_subplot( volume, [x,y,z,s] )
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
