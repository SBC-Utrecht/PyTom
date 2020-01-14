import sys
import os

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui, QtWidgets


class Dance(QMainWindow):
    def __init__(self, parent=None):
        super(Dance,self).__init__(parent)
        self.button = QPushButton('Run )')
        self.label = QLabel('Dummy')
        self.button.clicked.connect(self.update)

    def update(self):
        self.label.setText('Updated')


def test_hello(qtbot):
    '''
    pytompath = '/Users/gijs/Documents/PyTomPrivate'

    print(pytompath)
    if not pytompath:
        print('Pytom package is not available. Please load, or install Pytom.')
        sys.exit()

    update_vars = False
    for search in ('LD_LIBRARY_PATH','PATH','PYTHONPATH'):
        # Check if env vars include all paths set in paths.csh
        query_string = "cat {}/bin/paths.csh | grep 'setenv {}' | grep -v '${}'".format(pytompath, search,search)
        string = os.popen(query_string).read()[:-1].split()[-1]
        for new_lib in (string.split(':')):
            new_lib = new_lib.replace("'","")

            if not new_lib in os.environ[search].split(':'):
                os.environ[search] += ':'+new_lib
                update_vars = True
    #If any of the env vars are updated reopen this script.
    if update_vars:
        if len(sys.argv) < 2:
            pythonVersion = 'python{d[0]}.{d[1]}'.format( d=sys.version_info )
            path = os.popen('which {}'.format(pythonVersion)).read()[:-1]
            sys.argv = [path] + sys.argv

        os.execv(sys.argv[0],sys.argv)
    '''


    widget = Dance()
    qtbot.addWidget(widget)

    # click in the Greet button and make sure it updates the appropriate label
    #qtbot.mouseClick(widget.button, QtCore.Qt.LeftButton)
    
    assert widget.label.text() == "Dummy"

def test_updated(qtbot):
    widget = Dance()
    qtbot.mouseClick(widget.button, QtCore.Qt.LeftButton)
    assert widget.label.text() == "Updated"

