import sys
from PyQt5 import QtGui, QtWidgets, QtCore
import numpy as np
import pandas as pd
from PyQt5.QtWidgets import QApplication, QWidget
import classes as fsc
import charts
import prop
import fs
import os
import time

QApplication.setAttribute(
    QtCore.Qt.AA_EnableHighDpiScaling, True)  # enable high-dpi scaling
QApplication.setAttribute(
    QtCore.Qt.AA_UseHighDpiPixmaps, True)  # use high-dpi icons


class RC(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.dpath = prop.dpath
        self.tpath = prop.tpath
        self.fnd = 'fdata'
        self.dat_file = os.path.join(self.tpath, self.fnd + '.csv')
        self.fnl = 'floads'
        self.loads_file = os.path.join(self.tpath, self.fnl + '.csv')
        self.head = ['Mx', 'My', 'N']
        self.tloads = pd.read_csv(self.loads_file, usecols=self.head, sep=';')
        self.iF = 0
        self.msgBox = QtWidgets.QMessageBox()
        self.htbox = QtWidgets.QHBoxLayout()
        self.vbox = QtWidgets.QVBoxLayout()
        bar = QtWidgets.QMenuBar()
        file = bar.addMenu("File")
        fopen = file.addMenu("Open")
        fopen.addAction("Data")
        fopen.addAction("Loads")
        file.addAction("Generate")
        file.addAction("Quit")
        analysis = bar.addMenu("Analysis")
        analysis.addAction("All Loads")
        analysis.addAction("One Load")
        chrts = bar.addMenu("Plots")
        chprops = chrts.addMenu("Strain-Stress")
        chprops.addAction("Concrete")
        chprops.addAction("Rebar")
        chres = chrts.addMenu("Results")
        chres.addAction("Strain")
        chres.addAction("Stress")
        chres.addAction("Stress3D")
        help = bar.addMenu("Help")
        help.addAction("Help")
        help.addAction("About")
        file.triggered[QtWidgets.QAction].connect(self.processtrigger)
        analysis.triggered[QtWidgets.QAction].connect(self.processtrigger)
        chrts.triggered[QtWidgets.QAction].connect(self.processtrigger)
        help.triggered[QtWidgets.QAction].connect(self.processtrigger)
        self.vbox.addWidget(bar)
        self.hbox = QtWidgets.QHBoxLayout()
        self.cb1 = QtWidgets.QComboBox()
        self.cb1.setToolTip('Select member')
        self.cb2 = QtWidgets.QComboBox()
        self.cb2.setToolTip("Poisson's ratio for slab")
        self.cb3 = QtWidgets.QComboBox()
        self.cb3.setToolTip("Yb3")
        self.chb1 = QtWidgets.QCheckBox()
        self.chb2 = QtWidgets.QCheckBox()
        self.chb1.setToolTip('Whether do analysis second limit state')
        self.chb2.setToolTip('Whether do second stage analysis')
        self.qtwi = QtWidgets.QTableWidgetItem
        self.qtab = QtWidgets.QTableWidget()
        self.qtab.setToolTip('Select Load Case')
        self.button1 = QtWidgets.QPushButton("Run analysis")
        self.button2 = QtWidgets.QPushButton("View report")
        self.button3 = QtWidgets.QPushButton("Ok")
        self.button4 = QtWidgets.QPushButton("Chart")
        self.button4.setToolTip('Interaction Surface for Section')
        self.htbox.addWidget(self.cb1)
        self.htbox.addWidget(self.button4)
        self.button4.setVisible(False)
        self.htbox.addWidget(self.cb2)
        self.htbox.addWidget(self.chb1)
        self.htbox.addWidget(self.cb3)
        self.chb1.setText('2-st lim_st')
        self.chb2.setText('2-st stage')
        self.chb1.setChecked(False)
        self.hbox.addWidget(self.button1)
        self.hbox.addWidget(self.chb2)
        self.hbox.addWidget(self.button2)
        self.hbox.addWidget(self.button3)
        self.setLayout(self.vbox)
        self.vbox.addLayout(self.htbox)
        self.vbox.addWidget(self.qtab)
        self.vbox.addLayout(self.hbox)
        self.cb1.addItems(["frame", "slab"])
        self.S = self.cb1.currentText()
        self.cb2.addItems(["0.0", "0.2"])
        self.cb3.addItems(["1.0", "0.85"])
        self.cb2.setVisible(False)
        self.ltable(self.tloads.values, self.head)
        self.cb1.currentIndexChanged.connect(lambda: self.on_changed())
        self.cb2.currentIndexChanged.connect(self.on_sel_v)
        self.cb3.currentIndexChanged.connect(self.on_sel_gb3)
        self.button1.clicked.connect(self.on_clicked)
        self.button2.clicked.connect(self.on_clicked)
        self.button3.clicked.connect(self.on_clicked)
        self.button4.clicked.connect(self.on_clicked)
        self.qtab.cellClicked.connect(self.on_selected)
        self.qtab.itemSelectionChanged.connect(self.on_selected)
        self.qtab.selectRow(self.iF)
        self.chb1.toggled.connect(lambda: self.btnstate(self.chb1))
        self.chb2.toggled.connect(lambda: self.btnstate(self.chb2))
        self.chb2.setVisible(False)
        self.setWindowTitle("Reinforced Concrete")
        ico = QtGui.QIcon("python.png")
        self.setWindowIcon(ico)
        self.show()

    def openfile(self, title, type, pref2):
        os.system('cls')
        pref1 = self.S[0]
        mask = type + ' (' + pref1 + pref2 + '*.csv)'
        filename = QtWidgets.QFileDialog.getOpenFileName(
            self, title, self.dpath, mask)[0]
        fn = (filename.split('/')[-1])[:-4]
        if filename == '':
            fn = pref1 + pref2 + '.csv'
            filename = os.path.join(self.tpath, fn)
        self.stfs, self.prpr = np._NoValue, np._NoValue
        return filename, fn

    def processtrigger(self, q):
        sender = q.text()
        if sender == 'Data':
            self.dat_file, self.fnd = self.openfile(
                'Open data file', 'Data files', 'data')
        if sender == 'Loads':
            self.loads_file, self.fnl = self.openfile(
                'Open loads file', 'Load files', 'loads')
            self.changetable(self.chb1)
        if sender == 'Generate':
            os.system('cls')
            fn = 'GuiGen.py'
            if os.path.exists(fn):
                os.system(fn)
            else:
                print('File', fn, "doesn't exist !")
        if sender == 'Quit':
            os.system('cls')
            sys.exit(1)
        elif sender == 'All Loads':
            lst2 = self.chb1.checkState()
            gb3 = float(self.cb3.currentText())
            self.checkall(self.S, lst2, gb3)
        elif sender == 'One Load':
            self.analysis()
        elif sender == 'Concrete':
            try:
                charts.Plots.strainstress(
                    self.stfs, self.prpr, 'concrete')
            except:
                os.system('cls')
                print('First please run analysis')
        elif sender == 'Rebar':
            try:
                charts.Plots.strainstress(
                    self.stfs, self.prpr, 'rebar')
            except:
                os.system('cls')
                print('First please run analysis')
        elif sender == 'Strain':
            if self.all != False:
                charts.Plots.strains(self.result)
            elif self.S == 'frame':
                try:
                    charts.Plots.f_strain(self.rsl, self.stfs)
                except:
                    os.system('cls')
                    print('First please run analysis')
            else:
                try:
                    charts.Plots.sc_strain(self.rsl, self.stfs)
                except:
                    os.system('cls')
                    print('First please run analysis')
        elif sender == 'Stress':
            if self.S == 'frame':
                try:
                    charts.Plots.f_stress(self.rsl, self.stfs)
                except:
                    os.system('cls')
                    print('First please run analysis')
            else:
                try:
                    charts.Plots.sc_stress(self.rsl, self.stfs)
                except:
                    os.system('cls')
                    print('First please run analysis')
        elif sender == 'Stress3D':
            if self.S == 'frame':
                try:
                    charts.Plots.f_stress3D(self.rsl, self.stfs)
                except:
                    os.system('cls')
                    print('First please run analysis')
            else:
                print('Have to be developed')
        elif sender == 'Help':
            os.popen("RC.chm")
        elif sender == 'About':
            self.about()

    def ltable(self, l, h):
        r = l.shape[0]
        c = l.shape[1]
        self.qtab.setRowCount(r)
        self.qtab.setColumnCount(c)
        self.qtab.setHorizontalHeaderLabels((h))

        self.qtab.setVerticalHeaderLabels(
            ['L/C ' + str(n) + '  ' for n in np.arange(r)])
        for i in range(0, r):
            j = 0
            for j in range(0, c):
                cellinfo = self.qtwi(str(l[i, j]))
                self.qtab.setItem(i, j, cellinfo)
                j += 1
            i += 1

    def changetable(self, b):
        self.stfs, self.prpr = np._NoValue, np._NoValue
        if self.S == 'frame':
            self.cb2.setVisible(False)
            head1 = ['Mx', 'My', 'N']
            head2 = ['Mx2', 'My2', 'N2']
        else:
            self.cb2.setVisible(True)
            self.button4.setVisible(False)
            head1 = ['Nxx', 'Nyy', 'Nxy', 'Mxx', 'Myy', 'Mxy']
            head2 = ['Nxx2', 'Nyy2', 'Nxy2', 'Mxx2', 'Myy2', 'Mxy2']
        if b.isChecked() == False:
            self.head = head1
        else:
            self.head = head2
        self.tloads = pd.read_csv(self.loads_file, usecols=self.head, sep=';')
        self.ltable(self.tloads.values, self.head)
        self.qtab.selectRow(0)

    def on_changed(self):
        self.chb2.setChecked(False)
        self.chb2.setVisible(False)
        os.system('cls')
        self.all = False
        self.S = self.cb1.currentText()
        try:
            self.dat_file, self.fnd = self.openfile(
                'Open data file', 'Data files', 'data')
        except:
            self.dat_file = os.path.join(self.dpath, self.S[0] + 'data.csv')
        try:
            self.loads_file, self.fnl = self.openfile(
                'Open loads file', 'Load files', 'loads')
        except:
            self.loads_file = os.path.join(self.dpath, self.S[0] + 'loads.csv')
        if self.S == 'frame':
            self.cb2.setVisible(False)
        else:
            self.cb2.setVisible(True)
            self.button4.setVisible(False)
        self.changetable(self.chb1)

    def on_sel_gb3(self):
        self.chb2.setChecked(False)
        self.chb2.setVisible(False)
        os.system('cls')
        self.stfs, self.prpr = np._NoValue, np._NoValue

    def on_sel_v(self):
        self.chb2.setChecked(False)
        self.chb2.setVisible(False)
        os.system('cls')
        self.stfs, self.prpr = np._NoValue, np._NoValue

    def checkall(self, m, lst2, gb3):
        self.all = True
        lst2 = self.chb1.checkState()
        stg2 = self.chb2.checkState()
        gb3 = float(self.cb3.currentText())
        fss = fs.Start(prop)
        self.load = pd.read_csv(self.loads_file, sep=';')
        self.stfs, self.prpr, self.result = fss.strt(
            lst2, stg2, gb3, m, self.dat_file, self.load)

    def analysis(self):
        self.chb2.setVisible(False)
        self.chb2.setVisible(True)
        self.all = False
        pd.set_option('display.max_rows', None)
        lst2 = self.chb1.checkState()
        stg2 = self.chb2.checkState()
        gb3 = float(self.cb3.currentText())
        self.load = pd.read_csv(self.loads_file, sep=';')
        self.fname = self.fnd + '_' + self.fnl + '-' + str(self.iF)
        if self.S == 'frame':
            import fclasses
            os.system('cls')
            self.fc = fclasses.Start(fsc, prop, self.load, self.dat_file, self.fname)
            self.stfs, self.prpr, self.rsl = self.fc.strt(
                time, lst2, stg2, gb3, self.iF)
            self.button4.setVisible(True)
        else:
            import sclasses
            os.system('cls')
            v = float(self.cb2.currentText())
            self.sc = sclasses.Start(fsc, prop, self.load, self.dat_file, v, self.fname)
            self.stfs, self.prpr, self.rsl = self.sc.strt(
                time, lst2, stg2, gb3, self.iF)
        self.chb2.setChecked(False)

    def on_clicked(self):
        sender = self.sender()
        if sender.text() == 'Run analysis':
            self.analysis()
        elif sender.text() == 'View report':
            os.system('cls')
            repname = self.fname + '.txt'
            if os.path.exists(os.path.join(self.dpath, repname)):
                f = open(os.path.join(self.dpath, repname), 'r')
                print(f.read())
                f.close()
            else:
                print('First please run analysis')

        elif sender.text() == 'Chart':
            print('Please wait:')
            self.chart_inter()

        elif sender.text() == 'About':
            self.about()
        else:
            self.close()

    def about(self):
        self.msgBox.setIcon(self.msgBox.Information)
        self.msgBox.setText('Program =RC_v1.0= designed by G. Berezin in 2020')
        self.msgBox.setWindowTitle("About program")
        self.msgBox.setStandardButtons(self.msgBox.Ok)
        self.msgBox.exec()

    def on_selected(self):
        self.iF = self.qtab.currentRow()

    def btnstate(self, b):
        self.stfs, self.prpr = np._NoValue, np._NoValue
        if b.text() == "2-st lim_st":
            self.changetable(self.chb1)
            if b.isChecked() == True:
                self.chb2.setChecked(False)
        if b.text() == "2-st stage":
            if b.isChecked() == True:
                self.chb1.setChecked(False)

    def chart_inter(self):
        self.fc.prpr.lim_st = '1_st'
        self.fc.lds.ld()
        Fi = self.fc.intrn.interact()
        charts.Plots.chart_I(Fi, self.fc.lds.F[self.iF, :])

    def text(self):
        self.sender().text()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = RC()
    sys.exit(app.exec_())
