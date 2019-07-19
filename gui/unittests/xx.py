
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets
import os
import sys

from gui import PyTomGui


def create_project(qtbot,name,wtype=''):
    widget = PyTomGui()
    qtbot.addWidget(widget)
    widget.open_project(name)
    #widget.label.setText(name)
    wtypes = {'CD': widget.CD, 'TR': widget.TR, 'SA': widget.SA, 'PP': widget.PP}
    window = wtypes[wtype]
    return widget, window

def set_values(window, dd, mode):
    for key, value in dd:
        if type(window.widgets[mode + key]) in (QSpinBox, QDoubleSpinBox):
            window.widgets[mode + key].setValue(value)
        elif type(window.widgets[mode + key]) == QLineEdit:
            window.widgets[mode + key].setText(str(value))

def TR_starting_values(window):
    for key in window.widgets.keys():
        tempWidget = window.widgets[key]
        if type(tempWidget) == QLineEdit:
            assert tempWidget.text().replace(' ', '') == ''
            continue

def test_start_gui(qtbot):
    widget = PyTomGui()
    qtbot.addWidget(widget)
    assert widget.qcommand == "sbatch"

def test_Alignment_form(qtbot):
    mode = 'v02_Align_'
    wbpForm = open('alignment.sh', 'r').read()
    widget, window = create_project(qtbot, 'UnitTestProject', 'TR')
    qtbot.addWidget(widget)

    dd = (("LastAngle", 60), ("LastIndex", 36), ("RefMarkerIndex", 1), ("RefTiltIndex", 18), ("RotationTiltAxis", 0),
          ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", "UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    assert window.widgets[mode + 'CommandText'].toPlainText() == wbpForm

def test_INFR_form(qtbot):
    mode = 'v02_INFR_'
    infrForm = open('INFR_reconstruction.sh', 'r').read()
    widget, window = create_project(qtbot, 'UnitTestProject', 'TR')

    dd = (("LastAngle", 60), ("LastIndex", 37), ("RefMarkerIndex", 1), ("RefTiltIndex", 19), ("RotationTiltAxis", 0),
          ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", "UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    assert window.widgets[mode + 'CommandText'].toPlainText() == infrForm

def test_WBP_form(qtbot):
    mode = 'v02_WBP_'
    wbpForm = open('WBP_reconstruction.sh', 'r').read()
    widget, window = create_project(qtbot, 'UnitTestProject', 'TR')
    qtbot.addWidget(widget)

    dd = (("LastAngle", 60), ("LastIndex", 37), ("RefMarkerIndex", 1), ("RefTiltIndex", 18), ("RotationTiltAxis", 0),
          ("WeightingType", 1), ("queue", False),("FirstAngle", -60), ("FirstIndex", 0),
          ("FolderSorted", "UnitTestProject/03_Tomographic_Reconstruction/tomogram_001/sorted"))

    set_values(window, dd, mode)
    qtbot.mouseClick(window.widgets[mode + 'GeneratePushButton'], QtCore.Qt.LeftButton)

    assert window.widgets[mode + 'CommandText'].toPlainText() == wbpForm