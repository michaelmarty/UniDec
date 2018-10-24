# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'iFAMS_V5_New.ui'
#
# Created by: PyQt5 UI code generator 5.10.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets
import os
from numpy import arange, sin, pi
import GuiTestFun
import numpy as np
from operator import itemgetter
import pyqtgraph as pg

pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')
pg.setConfigOption('leftButtonPan', False)


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1148, 839)
        MainWindow.setStyleSheet(
            "background-color: qlineargradient(spread:pad, x1:0, y1:0, x2:1, y2:0, stop:0 rgba(198, 198, 198, 255), stop:1 rgba(255, 255, 255, 255));R")
        MainWindow.setWindowState(QtCore.Qt.WindowMaximized)
        QtGui.QApplication.setStyle(QtGui.QStyleFactory.create("Plastique"))
        MainWindow.setWindowIcon(QtGui.QIcon('CapI.png'))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.frame = QtWidgets.QFrame(self.centralwidget)
        self.frame.setFrameShape(QtWidgets.QFrame.Box)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setLineWidth(4)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.widget = pg.PlotWidget(self.frame)
        self.widget.setObjectName("widget")
        self.widget.hide()
        self.gridLayout.addWidget(self.widget, 0, 0, 1, 1)
        self.widget_3 = pg.PlotWidget(self.frame)
        self.widget_3.hide()
        self.widget_3.setObjectName("widget_3")
        self.gridLayout.addWidget(self.widget_3, 0, 1, 1, 1)
        self.widget_2 = pg.PlotWidget(self.frame)
        self.widget_2.hide()
        self.widget_2.setObjectName("widget_2")
        self.gridLayout.addWidget(self.widget_2, 1, 0, 1, 1)
        self.widget_4 = pg.PlotWidget(self.frame)
        self.widget_4.setObjectName("widget_4")
        self.widget_4.hide()
        self.gridLayout.addWidget(self.widget_4, 1, 1, 1, 1)
        self.gridLayout_3.addWidget(self.frame, 0, 0, 1, 1)
        self.frame_2 = QtWidgets.QFrame(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame_2.sizePolicy().hasHeightForWidth())
        self.frame_2.setSizePolicy(sizePolicy)
        self.frame_2.setMaximumSize(QtCore.QSize(300, 16777215))
        self.frame_2.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setLineWidth(4)
        self.frame_2.setObjectName("frame_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.frame_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.toolBox = QtWidgets.QToolBox(self.frame_2)
        self.toolBox.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.toolBox.setObjectName("toolBox")
        self.toolBox.setMinimumWidth(270)
        self.page = QtWidgets.QWidget()
        self.page.setGeometry(QtCore.QRect(0, 0, 262, 604))
        self.page.setObjectName("page")
        self.widget1 = QtWidgets.QWidget(self.page)
        self.widget1.setGeometry(QtCore.QRect(1, 1, 261, 325))
        self.widget1.setObjectName("widget1")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.widget1)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.widget1)
        self.label.setMaximumSize(QtCore.QSize(16777215, 30))
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.lineEdit = QtWidgets.QLineEdit(self.widget1)
        self.lineEdit.setObjectName("lineEdit")
        self.verticalLayout.addWidget(self.lineEdit)
        self.label_2 = QtWidgets.QLabel(self.widget1)
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 30))
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.lineEdit_2 = QtWidgets.QLineEdit(self.widget1)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.verticalLayout.addWidget(self.lineEdit_2)
        self.label_3 = QtWidgets.QLabel(self.widget1)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.lineEdit_3 = QtWidgets.QLineEdit(self.widget1)
        self.lineEdit_3.setObjectName("lineEdit_3")
        self.verticalLayout.addWidget(self.lineEdit_3)
        self.pushButton = QtWidgets.QPushButton(self.widget1)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.Recalculate)
        self.verticalLayout.addWidget(self.pushButton)
        self.label_7 = QtWidgets.QLabel(self.widget1)
        self.label_7.setObjectName("label_7")
        self.verticalLayout.addWidget(self.label_7)
        self.lineEdit_7 = QtWidgets.QLineEdit(self.widget1)
        self.lineEdit_7.setObjectName("lineEdit_7")
        self.verticalLayout.addWidget(self.lineEdit_7)
        self.label_8 = QtWidgets.QLabel(self.widget1)
        self.label_8.setObjectName("label_8")
        self.verticalLayout.addWidget(self.label_8)
        self.lineEdit_8 = QtWidgets.QLineEdit(self.widget1)
        self.lineEdit_8.setObjectName("lineEdit_8")
        self.verticalLayout.addWidget(self.lineEdit_8)
        self.pushButton_3 = QtWidgets.QPushButton(self.widget1)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.FourierNoise)
        self.verticalLayout.addWidget(self.pushButton_3)
        self.toolBox.addItem(self.page, "")
        self.page_3 = QtWidgets.QWidget()
        self.page_3.setObjectName("page_3")
        self.widget2 = QtWidgets.QWidget(self.page_3)
        self.widget2.setGeometry(QtCore.QRect(0, 0, 261, 201))
        self.widget2.setObjectName("widget2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.widget2)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_4 = QtWidgets.QLabel(self.widget2)
        self.label_4.setObjectName("label_4")
        self.verticalLayout_2.addWidget(self.label_4)
        self.lineEdit_4 = QtWidgets.QLineEdit(self.widget2)
        self.lineEdit_4.setObjectName("lineEdit_4")
        self.verticalLayout_2.addWidget(self.lineEdit_4)
        self.label_5 = QtWidgets.QLabel(self.widget2)
        self.label_5.setObjectName("label_5")
        self.verticalLayout_2.addWidget(self.label_5)
        self.lineEdit_5 = QtWidgets.QLineEdit(self.widget2)
        self.lineEdit_5.setObjectName("lineEdit_5")
        self.verticalLayout_2.addWidget(self.lineEdit_5)
        self.label_6 = QtWidgets.QLabel(self.widget2)
        self.label_6.setObjectName("label_6")
        self.verticalLayout_2.addWidget(self.label_6)
        self.lineEdit_6 = QtWidgets.QLineEdit(self.widget2)
        self.lineEdit_6.setObjectName("lineEdit_6")
        self.verticalLayout_2.addWidget(self.lineEdit_6)
        self.pushButton_2 = QtWidgets.QPushButton(self.widget2)
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.clicked.connect(self.inputs)
        self.verticalLayout_2.addWidget(self.pushButton_2)
        self.toolBox.addItem(self.page_3, "")
        self.page_5 = QtWidgets.QWidget()
        self.page_5.setObjectName("page_5")
        self.widget3 = QtWidgets.QWidget(self.page_5)
        self.widget3.setGeometry(QtCore.QRect(0, 0, 261, 82))
        self.widget3.setObjectName("widget3")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.widget3)
        self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_9 = QtWidgets.QLabel(self.widget3)
        self.label_9.setObjectName("label_9")
        self.verticalLayout_3.addWidget(self.label_9)
        self.lineEdit_9 = QtWidgets.QLineEdit(self.widget3)
        self.lineEdit_9.setObjectName("lineEdit_9")
        self.verticalLayout_3.addWidget(self.lineEdit_9)
        self.pushButton_4 = QtWidgets.QPushButton(self.widget3)
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_4.clicked.connect(self.AverageHarm)
        self.verticalLayout_3.addWidget(self.pushButton_4)
        self.label_9.raise_()
        self.lineEdit_9.raise_()
        self.pushButton_4.raise_()
        self.pushButton_4.raise_()
        self.toolBox.addItem(self.page_5, "")
        self.page_4 = QtWidgets.QWidget()
        self.page_4.setObjectName("page_4")
        self.widget4 = QtWidgets.QWidget(self.page_4)
        self.widget4.setGeometry(QtCore.QRect(8, 3, 241, 134))
        self.widget4.setObjectName("widget4")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.widget4)
        self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.label_11 = QtWidgets.QLabel(self.widget4)
        self.label_11.setObjectName("label_11")
        self.verticalLayout_5.addWidget(self.label_11)
        self.lineEdit_11 = QtWidgets.QLineEdit(self.widget4)
        self.lineEdit_11.setObjectName("lineEdit_11")
        self.verticalLayout_5.addWidget(self.lineEdit_11)
        self.label_12 = QtWidgets.QLabel(self.widget4)
        self.label_12.setObjectName("label_12")
        self.verticalLayout_5.addWidget(self.label_12)
        self.lineEdit_12 = QtWidgets.QLineEdit(self.widget4)
        self.lineEdit_12.setObjectName("lineEdit_12")
        self.verticalLayout_5.addWidget(self.lineEdit_12)
        self.pushButton_6 = QtWidgets.QPushButton(self.widget4)
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_6.clicked.connect(self.Filter)
        self.verticalLayout_5.addWidget(self.pushButton_6)
        self.toolBox.addItem(self.page_4, "")
        self.gridLayout_2.addWidget(self.toolBox, 0, 1, 1, 1)
        self.gridLayout_3.addWidget(self.frame_2, 0, 1, 1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1148, 26))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuSave = QtWidgets.QMenu(self.menuFile)
        self.menuSave.setObjectName("menuSave")
        self.menuTools = QtWidgets.QMenu(self.menubar)
        self.menuTools.setObjectName("menuTools")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionLoad_Spectrum = QtWidgets.QAction(MainWindow)
        self.actionLoad_Spectrum.setObjectName("actionLoad_Spectrum")
        self.actionFourier_Spectrum = QtWidgets.QAction(MainWindow)
        self.actionFourier_Spectrum.setObjectName("actionFourier_Spectrum")
        self.actionEnvelope_Reconstructions = QtWidgets.QAction(MainWindow)
        self.actionEnvelope_Reconstructions.setObjectName("actionEnvelope_Reconstructions")
        self.actionZero_Charge_Spectrum = QtWidgets.QAction(MainWindow)
        self.actionZero_Charge_Spectrum.setObjectName("actionZero_Charge_Spectrum")
        self.actionFourier_Filtered_Spectrum = QtWidgets.QAction(MainWindow)
        self.actionFourier_Filtered_Spectrum.setObjectName("actionFourier_Filtered_Spectrum")
        self.actionQuit = QtWidgets.QAction(MainWindow)
        self.actionQuit.setObjectName("actionQuit")
        self.actionFourier_Transform = QtWidgets.QAction(MainWindow)
        self.actionFourier_Transform.setObjectName("actionFourier_Transform")
        self.actionCalculate_Charge_States_and_Subunit_Mass = QtWidgets.QAction(MainWindow)
        self.actionCalculate_Charge_States_and_Subunit_Mass.setObjectName(
            "actionCalculate_Charge_States_and_Subunit_Mass")
        self.actionReconstruct_Envelope_Functions = QtWidgets.QAction(MainWindow)
        self.actionReconstruct_Envelope_Functions.setObjectName("actionReconstruct_Envelope_Functions")
        self.actionModel_Spectrum_Generator = QtWidgets.QAction(MainWindow)
        self.actionModel_Spectrum_Generator.setObjectName("actionModel_Spectrum_Generator")
        self.actionSubunit_Statistics_Calculator = QtWidgets.QAction(MainWindow)
        self.actionSubunit_Statistics_Calculator.setObjectName("actionSubunit_Statistics_Calculator")
        self.actionRun_iFAMS_Analysis = QtWidgets.QAction(MainWindow)
        self.actionRun_iFAMS_Analysis.setObjectName("actionRun_iFAMS_Analysis")
        self.actionModel_Spectrum_Generator_2 = QtWidgets.QAction(MainWindow)
        self.actionModel_Spectrum_Generator_2.setObjectName("actionModel_Spectrum_Generator_2")
        self.actionSubunit_Statistic_Calculator = QtWidgets.QAction(MainWindow)
        self.actionSubunit_Statistic_Calculator.setObjectName("actionSubunit_Statistic_Calculator")
        self.menuSave.addAction('Fourier Spectrum', self.SaveFT)
        self.menuSave.addAction('Envelope_Reconstructions', self.SaveEnvel)
        self.menuSave.addAction('Zero_Charge_Spectrum', self.SaveZero)
        self.menuSave.addAction('Fourier_Filtered_Spectrum', self.SaveFilter)
        self.menuFile.addAction('Load Spectrum', self.MSplot)
        self.menuFile.addAction(self.menuSave.menuAction())
        self.menuFile.addSeparator()
        self.menuFile.addAction('&Quit', self.fileQuit,
                                QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuTools.addAction('Fourier Transform', self.FourierPlot)
        self.menuTools.addAction('Calculate Charge States and Subunit mass', self.CharAndSubCalc)
        self.menuTools.addAction('Run iFAMS Analysis', self.reconstruct)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuTools.menuAction())

        self.retranslateUi(MainWindow)
        self.toolBox.setCurrentIndex(0)
        self.toolBox.setStyleSheet(
            "QToolBox::tab{color: rgb(255, 255, 255);background-color: rgb(74, 74, 144)}QToolBox::tab:selected{color: rgb(255,255,255);background-color: rgb(5, 208, 70)}")
        self.toolBox.setItemText(0, "Fourier Transform")
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "iFAMS Version 5"))
        self.label.setText(_translate("MainWindow", "Minimum Frequency Used"))
        self.label_2.setText(_translate("MainWindow", "Minimum Peak Percentage"))
        self.label_3.setText(_translate("MainWindow", "Delta"))
        self.pushButton.setText(_translate("MainWindow", "Calculate"))
        self.label_7.setText(_translate("MainWindow", "Minimum Noise Frequency"))
        self.label_8.setText(_translate("MainWindow", "Maximum Noise Frequency"))
        self.pushButton_3.setText(_translate("MainWindow", "Calculate"))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page), _translate("MainWindow", " Fourier Spectrum"))
        self.label_4.setText(_translate("MainWindow", "Low Charge State"))
        self.label_5.setText(_translate("MainWindow", "High Charge State"))
        self.label_6.setText(_translate("MainWindow", "Subunit Mass"))
        self.pushButton_2.setText(_translate("MainWindow", "Calculate"))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_3),
                                 _translate("MainWindow", "Manually Enter Charge States"))
        self.label_9.setText(_translate("MainWindow", "Harmonics Used"))
        self.pushButton_4.setText(_translate("MainWindow", "Calculate"))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_5), _translate("MainWindow", "Harmonic Average"))
        self.label_11.setText(_translate("MainWindow", "Zero Frequency Data"))
        self.label_12.setText(_translate("MainWindow", "Number of Harmonics Used"))
        self.pushButton_6.setText(_translate("MainWindow", "Calculate"))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_4), _translate("MainWindow", "Fourier Filter"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuSave.setTitle(_translate("MainWindow", "Save"))
        self.menuTools.setTitle(_translate("MainWindow", "Analysis"))
        self.actionLoad_Spectrum.setText(_translate("MainWindow", "Load Spectrum"))
        self.actionFourier_Spectrum.setText(_translate("MainWindow", "Fourier Spectrum"))
        self.actionEnvelope_Reconstructions.setText(_translate("MainWindow", "Envelope Reconstructions"))
        self.actionZero_Charge_Spectrum.setText(_translate("MainWindow", "Zero Charge Spectrum"))
        self.actionFourier_Filtered_Spectrum.setText(_translate("MainWindow", "Fourier Filtered Spectrum"))
        self.actionQuit.setText(_translate("MainWindow", "Quit"))
        self.actionFourier_Transform.setText(_translate("MainWindow", "Fourier Transform"))
        self.actionCalculate_Charge_States_and_Subunit_Mass.setText(
            _translate("MainWindow", "Calculate Charge States and Subunit Mass"))
        self.actionReconstruct_Envelope_Functions.setText(_translate("MainWindow", "Run iFAMS analysis"))
        self.actionModel_Spectrum_Generator.setText(_translate("MainWindow", "Fourier Transform"))
        self.actionSubunit_Statistics_Calculator.setText(
            _translate("MainWindow", "Calculate Charge States and Subunit Mass"))
        self.actionRun_iFAMS_Analysis.setText(_translate("MainWindow", "Run iFAMS Analysis"))

    def fileQuit(self):
        MainWindow.close()

    def closeEvent(self, event):

        reply = QtGui.QMessageBox.question(self, 'Message',
                                           "Are you sure you want to quit?", QtGui.QMessageBox.Yes |
                                           QtGui.QMessageBox.No, QtGui.QMessageBox.No)

        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def MSplot(self):  # A function that plots the mass spectrum
        self.widget.clear()
        self.widget_4.hide()
        self.widget_3.hide()
        self.widget_2.hide()
        self.xdata, self.ydata, self.namebase = GuiTestFun.plot(self)
        self.yfull, self.expandedspan, self.xnew, self.ftspacing, self.paddedxnew, self.ynew = GuiTestFun.plot_function(
            self.xdata, self.ydata)
        self.widget.plot(self.xdata, self.ydata, pen=pg.mkPen(0, 0, 128), symbol=None)
        self.widget.setLabel('left', 'abundance')
        self.widget.setLabel('bottom', 'm/z')
        self.widget.show()

    def FourierPlot(
            self):  # A function that Fourier transfroms the mass spectrum, a plots the fourier transform next to th
        self.widget_3.clear()
        self.maxfreq = GuiTestFun.maxfreq(self, self.yfull, self.expandedspan)
        self.ftx, self.ABFT, self.FT = GuiTestFun.Fourier(self, self.maxfreq, self.yfull)
        self.refmaxtab = GuiTestFun.findmax(self.expandedspan, self.ABFT, self.ftx, 0.001, 5, 10)
        self.widget_3.plot(self.ftx, self.ABFT, pen='k', symbol=None)
        self.widget_3.plot(np.array(self.refmaxtab)[:, 0], np.array(self.refmaxtab)[:, 1],
                           pen=None, symbolPen='r', symbolBrush='r', symbol='o')
        self.widget_3.setLabel('left', 'amplitude')
        self.widget_3.setLabel('bottom', 'frequency')
        self.widget_3.show()

    def Recalculate(self):  # A Fuction to estimate peak maxima in your spectra
        lowend = float(str(self.lineEdit.text()))
        pctpkht = float(str(self.lineEdit_2.text()))
        delta = int(str(self.lineEdit_3.text()))
        self.refmaxtab = GuiTestFun.findmax(self.expandedspan, self.ABFT, self.ftx, lowend, delta, pctpkht)
        self.plotDots()

    def plotDots(self):  # This plots the dots in the spectrum
        self.widget_3.clear()
        self.widget_3.plot(self.ftx, self.ABFT, pen='k', symbol=None)
        self.widget_3.plot(np.array(self.refmaxtab)[:, 0], np.array(self.refmaxtab)[:, 1],
                           pen=None, symbolPen='b', symbolBrush='b', symbol='o')
        self.widget_3.setLabel('left', 'amplitude')
        self.widget_3.setLabel('bottom', 'frequency')
        self.widget_3.show()

    def CharAndSubCalc(self):  # calculates an initial estimate for the spacing of the group of evenly-spaced FFT
        # spectrum peaks that probably aren't overtone peaks
        self.colors = ['#434483', '#4DC29C', '#9F23A4', '#CC7647', '#C2201F', '#C49DD4', '#B8A2BC']

        self.newcalcX = []
        self.newcalcY = []
        self.refmaxtabCalc = np.array(self.refmaxtab)
        self.numchar = GuiTestFun.spacing(self, self.refmaxtab)
        self.omega = GuiTestFun.omega(self.refmaxtab, self.numchar)
        self.chargestates, self.chargestatesr = GuiTestFun.charge(self.refmaxtab, self.numchar, self.omega)
        self.chargestateints = [int(self.chargestatesr[i]) for i in range(0, len(self.chargestatesr))]
        for i in range(0, len(self.chargestatesr)):
            self.newcalcX.append(self.refmaxtabCalc[i, 0])
            self.newcalcY.append(self.refmaxtabCalc[i, 1])
        print(self.newcalcX)
        self.widget_3.clear()
        self.widget_3.plot(self.ftx, self.ABFT, pen='k', symbol=None)
        for i in range(len(self.newcalcY)):
            self.widget_3_scat = pg.ScatterPlotItem([self.newcalcX[i]], [self.newcalcY[i]], size=10,
                                                    pen=pg.mkPen(None))  # brush=pg.mkBrush(255, 255, 255, 120))
            color = QtGui.QColor(0, 0, 0)
            color.setNamedColor(self.colors[i])
            self.widget_3_scat.setBrush(color)
            self.widget_3.addItem(self.widget_3_scat)
        # self.widget_3.plot(self.newcalcX,self.newcalcY,pen=None, symbol="o", symbolpen='w',symbolbrush=pg.mkPen(self.colors[i] for i in range(0,len(self.newcalcX)),hues=360))
        self.submass, self.stdevmass = GuiTestFun.subunit(self.refmaxtab, self.numchar, self.omega, self.chargestates,
                                                          self.chargestatesr, self.ftx, self.ABFT)
        self.widget_3info = "The charge states are " + str(self.chargestateints) + " with a subunit mass of " + str(
            self.submass) + " +/- " + str(self.stdevmass)
        self.widget_3info2 = pg.TextItem()
        self.widget_3info2.setText(self.widget_3info, 'r')
        self.widget_3info2.setTextWidth(300)
        self.widget_3.addItem(self.widget_3info2)
        self.widget_3info2.setPos(self.newcalcX[-1], max(self.newcalcY))
        self.widget_3info2.updateTextPos()

    def FourierNoise(self):  # the calculates the white noise in the fourier spectrum

        minnoi = int(np.floor(float(str(self.lineEdit_7.text())) * 2 * self.expandedspan))
        maxnoi = int(np.floor(float(str(self.lineEdit_8.text())) * 2 * self.expandedspan))

        ##################### Real Noise Average #######################
        self.realnoiseavg = GuiTestFun.FTRealNoiAvg(minnoi, maxnoi, self.FT)
        realnoiseAvgtag = "average Re(noise) is "
        realnoiseavgFin = realnoiseAvgtag + str(self.realnoiseavg)

        ##################### Real Noise Average #######################

        ##################### Imag Noise Average #######################
        self.imagnoiseavg = GuiTestFun.FTImagNoiAvg(minnoi, maxnoi, self.FT)
        imagnoiavg = "average Im(noise) is "
        imagnoiavg2 = imagnoiavg + str(self.imagnoiseavg)

        ##################### Imag Noise Average #######################

        ##################### Real Noise Stdev #########################
        self.realnoisestdev = GuiTestFun.FTRealNoiStdev(minnoi, maxnoi, self.FT)
        realnoisvd = "stdev of Re(noise) is "
        realnoisvd2 = realnoisvd + str(self.realnoisestdev)

        ##################### Real Noise Stdev #########################

        ##################### Imag Noise Stdev #########################
        self.imagnoisestdev = GuiTestFun.FTImagNoiStdev(minnoi, maxnoi, self.FT)
        imagenoisvd = "stdev of Im(noise) is "
        imagenoisvd2 = imagenoisvd + str(self.imagnoisestdev)

        ##################### Imag Noise Stdev #########################

    def reconstruct(self):  # this function will window the fourier spectrum to reconstruct the peak

        self.ABIFT, self.ABIFTnoise, self.ABIFTintegral, self.msintegral = GuiTestFun.IFTEnvelope(self.chargestatesr,
                                                                                                  self.expandedspan,
                                                                                                  self.submass,
                                                                                                  self.ftx,
                                                                                                  self.ftspacing,
                                                                                                  self.FT,
                                                                                                  self.paddedxnew,
                                                                                                  self.xnew,
                                                                                                  self.ynew, self.ydata,
                                                                                                  self.realnoisestdev,
                                                                                                  self.imagnoisestdev)

        self.widget.clear()
        self.widget.plot(self.xdata, self.ydata, pen='k', symbol=None)
        self.widget.setLabel('left', 'abundance')
        self.widget.setLabel('bottom', 'm/z')
        self.widget.setTitle('Color corresponds to charge state Dot in Fourier Domain')
        self.widget.show()

        for i in range(0, len(self.chargestateints)):
            color = QtGui.QColor(0, 0, 0)
            color.setNamedColor(self.colors[i])
            self.widget.plot(self.xnew, self.ABIFT[i][0:int(len(self.xnew))],
                             pen=pg.mkPen(color=color, width=3), symbol=None)
            self.xnewfornoise = []
            self.ABIFTplusnoise = []
            self.ABIFTminusnoise = []
            for j in range(0, len(self.ABIFT[i])):
                if j % 50 == 0:
                    self.ABIFTnoise[i][int(j / 50)] = self.ABIFTnoise[i][int(j / 50)]
            for m in range(0, len(self.ABIFTnoise[i])):
                self.xnewfornoise.append(self.paddedxnew[int(50 * m)])
                self.ABIFTplusnoise.append(self.ABIFT[i][int(50 * m)] + self.ABIFTnoise[i][m])
                self.ABIFTminusnoise.append(self.ABIFT[i][int(50 * m)] - self.ABIFTnoise[i][m])
            self.widget.plot(self.xnewfornoise[0:int(len(self.xnew) / 50)],
                             self.ABIFTplusnoise[0:int(len(self.xnew) / 50)],
                             pen=pg.mkPen(color=color, symbol=None, width=3, style=QtCore.Qt.DashLine))
            self.widget.plot(self.xnewfornoise[0:int(len(self.xnew) / 50)],
                             self.ABIFTminusnoise[0:int(len(self.xnew) / 50)],
                             pen=pg.mkPen(color=color, symbol=None, width=3, style=QtCore.Qt.DashLine))
        self.zerochargeplot()

    def zerochargeplot(self):
        self.widget_2.clear()
        self.widget_2.addLegend()
        self.xrange, self.yfinal, self.yrangespec = GuiTestFun.zerocharge(self.ABIFT, self.xnew, self.chargestatesr)
        for i in range(0, len(self.yrangespec)):
            color = QtGui.QColor(0, 0, 0)
            color.setNamedColor(self.colors[i])
            self.widget_2.plot(self.xrange, self.yrangespec[i], pen=pg.mkPen(color=color, width=3), symbol=None)
        self.widget_2.plot(self.xrange, self.yfinal, pen=pg.mkPen(color='k', width=3, symbol=None),
                           name="zero charge spectrum")
        self.widget_2.setLabel('left', 'abundance')
        self.widget_2.setLabel('bottom', 'mass (Da)')
        self.widget_2.setTitle('Color corresponds to charge state Dot in Fourier Domain')
        self.widget_2.show()

    def Filter(self):
        self.widget_4.clear()
        self.widget_4.addLegend()
        ZFreqData = int(str(self.lineEdit_11.text()))
        OTnum = int(str(self.lineEdit_12.text()))
        self.reconstspec, self.reconstbaseline = GuiTestFun.FFTFilter(self.expandedspan, self.submass, ZFreqData,
                                                                      self.ftx,
                                                                      self.ftspacing, self.FT, self.chargestatesr,
                                                                      OTnum)
        FiltDataName = "Fourier Filtered Spectrum"
        BaselineName = "Fourier Baseline"
        BaselineSubName = "Baseline Subtracted Data"
        self.widget_4.plot(self.xdata, self.ydata, pen='k', symbol=None)
        self.widget_4.plot(self.xnew, self.reconstspec[0:int(len(self.xnew))],
                           pen=pg.mkPen(color=(0, 204, 204), symbol=None, width=2), name=FiltDataName, symbol=None)
        self.widget_4.plot(self.xnew, self.reconstbaseline[0:int(len(self.xnew))],
                           pen=pg.mkPen(color='r', symbol=None, width=3), name=BaselineName, symbol=None)
        self.widget_4.plot(self.xnew, self.ynew - self.reconstbaseline[0:int(len(self.xnew))],
                           pen=(0, 153, 0), name=BaselineSubName, symbol=None)
        self.widget_4.show()

    def AverageHarm(self):
        self.widget.clear()
        self.widget.plot(self.xdata, self.ydata, pen='k', symbol=None)
        self.widget.setLabel('left', 'abundance')
        self.widget.setLabel('bottom', 'm/z')

        self.chargestateints = [int(self.chargestatesr[i]) for i in range(0, len(self.chargestatesr))]

        ov = int(str(self.lineEdit_9.text()))
        self.ABIFT = GuiTestFun.AverageHarmFun(self.chargestatesr, self.expandedspan, self.submass, self.ftx,
                                               self.ftspacing,
                                               self.FT, self.paddedxnew, self.xnew, self.ynew, self.ydata, ov)
        for i in range(0, len(self.chargestateints)):
            color = QtGui.QColor(0, 0, 0)
            color.setNamedColor(self.colors[i])
            self.widget.plot(self.xnew, self.ABIFT[i][0:int(len(self.xnew))],
                             pen=pg.mkPen(color=color, width=3), symbol=None)

        self.widget.show()
        self.zerochargeplot()

    def SaveFT(self):  # saves the FT spectrum
        ftstring = "FT"
        ftfilename = self.namebase + ftstring + ".csv"
        ftarrayforexport = np.transpose([self.ftx, self.ABFT])
        np.savetxt(ftfilename, ftarrayforexport, fmt='%10.6f', delimiter=',')  # outputs the FFT data to its own file

    def SaveEnvel(self):  # saves the envelope function
        for i in range(0, len(self.chargestatesr)):
            ifftfilename = 0
            ifftfilenamenoise = 0
            iftstring = "IFFT"
            ifftfilename = self.namebase + iftstring + str(self.chargestateints[i]) + ".csv"
            ifftforexport = 0
            ifftforexport = np.transpose([self.xnew, self.ABIFT[i][0:int(len(self.xnew))]])
            np.savetxt(ifftfilename, ifftforexport, fmt='%10.6f', delimiter=',')  # outputs each charge-state-specific
            # spectrum to its own csv file
            self.xnewfornoise = []
            self.ABIFTplusnoise = []
            self.ABIFTminusnoise = []
            for j in range(0, len(self.ABIFT[i])):
                if j % 50 == 0:
                    self.ABIFTnoise[i][int(j / 50)] = self.ABIFTnoise[i][int(j / 50)]
            for m in range(0, len(self.ABIFTnoise[i])):
                self.xnewfornoise.append(self.paddedxnew[int(50 * m)])
                self.ABIFTplusnoise.append(self.ABIFT[i][int(50 * m)] + self.ABIFTnoise[i][m])
                self.ABIFTminusnoise.append(self.ABIFT[i][int(50 * m)] - self.ABIFTnoise[i][m])
            ifftfilenamenoise = self.namebase + iftstring + str(self.chargestateints[i]) + "noise.csv"
            ifftnoiseforexport = 0
            ifftnoiseforexport = np.transpose([self.xnewfornoise[0:int(len(self.xnew) / 50)],
                                               self.ABIFTminusnoise[0:int(len(self.xnew) / 50)],
                                               self.ABIFTplusnoise[0:int(len(self.xnew) / 50)]])
            np.savetxt(ifftfilenamenoise, ifftnoiseforexport,
                       fmt='%10.6f', header='m/z,abundance-noise,abundance+noise', delimiter=',')

    def SaveFilter(self):
        FIltString1 = "Fourier_Filtered_Spectrum"
        FIltString2 = "Fourier_Baseline"
        FIltString3 = "Baseline_Subtracted_Data"
        filenameFil1 = self.namebase + FIltString1 + ".csv"
        filenameFil2 = self.namebase + FIltString2 + ".csv"
        filenameFil3 = self.namebase + FIltString3 + ".csv"
        filforexport1 = np.transpose([self.xnew, self.reconstspec[0:int(len(self.xnew))]])
        filforexport2 = np.transpose([self.xnew, self.reconstbaseline[0:int(len(self.xnew))]])
        filforexport3 = np.transpose([self.xnew, self.ynew - self.reconstbaseline[0:int(len(self.xnew))]])
        np.savetxt(filenameFil1, filforexport1, fmt='%10.6f', delimiter=',')
        np.savetxt(filenameFil2, filforexport2, fmt='%10.6f', delimiter=',')
        np.savetxt(filenameFil3, filforexport3, fmt='%10.6f', delimiter=',')

    def SaveZero(self):
        for i in range(0, len(self.chargestateints)):
            ifftfilename = 0
            iftstring = "zerocharge"
            ifftfilename = self.namebase + iftstring + str(self.chargestateints[i]) + ".csv"
            ifftforexport = 0
            ifftforexport = np.transpose([self.xrange, self.yrangespec[i]])
            np.savetxt(ifftfilename, ifftforexport, fmt='%10.6f',
                       delimiter=',')  # outputs each charge-state-specific spectrum to
            # its own csv file
            print("file was exported")
        ifftfilename = 0
        iftstring = "zerocharge"
        ifftfilename = self.namebase + iftstring + ".csv"
        ifftforexport = 0
        ifftforexport = np.transpose([self.xrange, self.yfinal])
        np.savetxt(ifftfilename, ifftforexport, fmt='%10.6f', delimiter=',')

    def inputs(self):  # the input for the
        lowcharge = int(str(self.lineEdit_4.text()))
        highcharge = int(str(self.lineEdit_5.text()))
        self.submass = float(str(self.lineEdit_6.text()))
        self.chargestatesr = np.arange(lowcharge, highcharge + 1)
        self.refmaxtab = GuiTestFun.inputmax(self.chargestatesr, self.submass, self.ABFT, self.ftx)
        self.plotDots()
        self.CharAndSubCalc()


if __name__ == "__main__":
    import sys

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
