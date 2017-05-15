from PyQt5 import QtGui, QtCore  # (the example applies equally well to PySide)
import sys
import pyqtgraph as pg
import numpy as np

## Always start by initializing Qt (only once per application)
#app = QtGui.QApplication([])

app = QtCore.QCoreApplication.instance()
if app is None:
    app = QtGui.QApplication(sys.argv)

## Define a top-level widget to hold everything
w = QtGui.QWidget()
w.move(400,100)
#w.showFullScreen()

# change background color
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

## Create some widgets to be placed inside
btn = QtGui.QPushButton('press me')
text = QtGui.QLineEdit('enter text')
listw = QtGui.QListWidget()
pltW1 = pg.PlotWidget()
#pw2 = pg.PlotWidget()

x = np.arange(0,1000,1)
noise = np.random.normal(0,1,1000)/1
y = np.sin(x/10)+noise

plt1 = pltW1.plot(x,y, color='g')


## Create a grid layout to manage the widgets size and position
layout = QtGui.QGridLayout()
w.setLayout(layout)

## Add widgets to the layout in their proper positions
layout.addWidget(btn, 0, 0)   # button goes in upper-left
layout.addWidget(text, 1, 0)   # text edit goes in middle-left
layout.addWidget(listw, 2, 0)  # list widget goes in bottom-left
layout.addWidget(pltW1, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows
#layout.addWidget(pw2, 3, 1, 2, 1)

## Display the widget as a new window
w.show()

## Start the Qt event loop
app.exec_()