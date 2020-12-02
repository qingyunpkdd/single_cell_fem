import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5.QtCore import *
from .winform import Ui_Form
from PyQt5.QtCore import pyqtSignal, Qt
from collections import OrderedDict
from functools import partial
import threading
import time


class MyMainWindow(QMainWindow, Ui_Form):
    # project_para = pyqtSignal(dict)
    status_run = pyqtSignal(dict)

    def __init__(self, parent=None, para_que=None, status_que=None):
        super(MyMainWindow, self).__init__(parent)
        self.setupUi(self)
        self.ui_info = OrderedDict()
        self.para_que = para_que
        self.status_que = status_que
        self.status_worker = Worker(parent=parent, status__que=self.status_que)
        self.init_ui()

    def init_ui(self):
        self.pushButton_3.clicked.connect(partial(self.update_dir,
                                                  info="please input expression data dir",
                                                  key_="data_expr_dir"
                                                  ))
        self.pushButton_2.clicked.connect(partial(self.update_dir,
                                                  info="please input gene set data dir",
                                                  key_="gmt_dir"
                                                  ))
        self.pushButton.clicked.connect(partial(self.update_dir,
                                                info="please input output data dir",
                                                key_="output_dir"
                                                ))
        self.pushButton_4.clicked.connect(self.send_win_info)

        self.status_worker.sinOut.connect(self.update_status)
        # self.updat_thr = Update_status(self)

        self.status_worker.start()
        self.spinBox.setValue(5)
        self.progressBar.setValue(0)
        self.timeEdit.setDisplayFormat("HH:mm:ss")
        self.timeEdit.setTime(QTime.currentTime())
        self.setWindowTitle("Single-Cell FEM Calcuator!")

    def send_win_info(self):
        self.ui_info["core_num"] = self.spinBox.value()
        self.para_que.put(self.ui_info)
        self.para_que.put({"finish": "None"})
        self.pushButton_4.setEnabled(False)

    def update_status(self, cont_dic: dict):
        for k, v in cont_dic.items():
            if k == "print_out":
                self.textBrowser.append(v)
                self.textBrowser.moveCursor(self.textBrowser.textCursor().End)
            if k == "percent":
                self.progressBar.setValue(v)
            if k == "end":
                # self.timeEdit.setTime(self.timeEdit.time())
                self.pushButton_4.setEnabled(True)

    def update_dir(self, info: str, key_: str, init_dir="../../"):
        dir_ = QFileDialog.getExistingDirectory(self, info, init_dir)
        self.ui_info[key_] = dir_
        if key_ == "data_expr_dir":
            self.lineEdit.setText(dir_)
        if key_ == "gmt_dir":
            self.lineEdit_2.setText(dir_)
        if key_ == "output_dir":
            self.lineEdit_3.setText(dir_)


class Worker(QThread):
    sinOut = pyqtSignal(dict)

    def __init__(self, parent=None, status__que=None):
        super(Worker, self).__init__(parent)
        self.working = True
        self.status__que = status__que

    def __del__(self):
        self.working = False
        self.wait()

    def run(self):
        while self.working:
            cont = self.status__que.get()
            self.sinOut.emit(cont)


