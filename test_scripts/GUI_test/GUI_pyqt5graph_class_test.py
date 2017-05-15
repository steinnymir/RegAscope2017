"""

GUI graph plotting test

author: S.Y.Agustsson
"""

from PyQt5 import QtGui, QtCore, QtWidgets  # (the example applies equally well to PySide)
import sys
import pyqtgraph as pg
import numpy as np


def rand(n):
    data = np.random.random(n)
    data[int(n*0.1):int(n*0.13)] += .5
    data[int(n*0.18)] += 2
    data[int(n*0.1):int(n*0.13)] *= 5
    data[int(n*0.18)] *= 20
    data *= 1e-12
    return data, np.arange(n, n+len(data)) / float(n)

#
#def updateData():
#    yd, xd = rand(10000)
#    p1.setData(y=yd, x=xd)


class PlottingGUI(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()




    def initUI(self):

        self.col = QtGui.QColor(0, 0, 0)
#        self.timer = QtCore.QTimer()
#        timer.timeout.connect(updateData)
#        timer.start(50)

        ## Define a top-level widget to hold everything
        #w = QtGui.QWidget()

        ## Create some widgets to be placed inside
        btn = QtGui.QPushButton('press me')
        btn.setCheckable(True)
        text = QtGui.QLineEdit('enter text')
        listw = QtGui.QListWidget()
        pw1 = pg.PlotWidget(name = 'plot1')
        #pw2 = pg.PlotWidget(name = 'plot2')

        ## Create a grid layout to manage the widgets size and position
        layout = QtGui.QGridLayout()
        self.setLayout(layout)

        ## Add widgets to the layout in their proper positions
        layout.addWidget(btn, 0, 0)   # button goes in upper-left
        layout.addWidget(text, 1, 0)   # text edit goes in middle-left
        layout.addWidget(listw, 2, 0)  # list widget goes in bottom-left
        layout.addWidget(pw1, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows
        #layout.addWidget(pw2, 3, 1, 2, 1)

        #data.updateData()
        #pw1.plotCurve(data)

        #p1 = pw1.plotXY()


    def plotXY(self):
        self.plot(x,y)


        self.show()




        #self.plot(np.random.normal(size=100), pen=(255,0,0), name="Red curve")

#    def updateData(self):
#        yd, xd = rand(10000)
#        self.setData(y=yd, x=xd)





if __name__ == '__main__':

    ## Always start by initializing Qt (only once per application)
    #app = QtGui.QApplication([])

    x = np.arange(0,1000,1)
    noise = np.random.normal(0,1,1000)/1
    y = np.sin(x/10)+noise

    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)
    plGUI = PlottingGUI()
    sys.exit(app.exec_())

    #pw1.plot(np.random.normal(size=100), pen=(255,0,0), name="Red curve")



#    ## Define a top-level widget to hold everything
#    w = QtGui.QWidget()
#
#    ## Create some widgets to be placed inside
#    btn = QtGui.QPushButton('press me')
#    text = QtGui.QLineEdit('enter text')
#    listw = QtGui.QListWidget()
#    plot = pg.PlotWidget()
#    otherplot = pg.PlotWidget()
#
#
#    ## Create a grid layout to manage the widgets size and position
#    layout = QtGui.QGridLayout()
#    w.setLayout(layout)
#
#    ## Add widgets to the layout in their proper positions
#    layout.addWidget(btn, 0, 0)   # button goes in upper-left
#    layout.addWidget(text, 1, 0)   # text edit goes in middle-left
#    layout.addWidget(listw, 2, 0)  # list widget goes in bottom-left
#    layout.addWidget(plot, 0, 1, 3, 1)  # plot goes on right side, spanning 3 rows
#    layout.addWidget(otherplot, 3, 1, 2, 1)

    ## Display the widget as a new window
    #w.show()

    ## Start the Qt event loop
