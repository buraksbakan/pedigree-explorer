import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout, QWidget, QMessageBox, QTableWidget, QTableWidgetItem
)
from PyQt6.QtWidgets import QMainWindow, QHBoxLayout, QWidget
import pyqtgraph as pg
from PyQt6.QtCore import Qt
from PyQt6.QtCore import pyqtSignal
from PyQt6.QtGui import QPainter, QColor, QPen, QBrush


class ChromosomeWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.regions = []   # list of (start, end)
        self.chrom_length = 1

    def set_regions(self, regions):
        """regions = list of (start, end) from BED"""
        self.regions = regions
        if regions:
            self.chrom_length = max(end for _, end in regions)
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)

        w = self.width()
        h = self.height()

        # Draw chromosome backbone
        painter.setPen(QPen(Qt.GlobalColor.black, 3))
        painter.drawLine(50, h//2, w-50, h//2)

        # Draw regions
        for start, end in self.regions:
            x1 = 50 + (start / self.chrom_length) * (w - 100)
            x2 = 50 + (end   / self.chrom_length) * (w - 100)

            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QBrush(QColor(100, 180, 255)))
            painter.drawRect(int(x1), h//2 - 20, int(x2 - x1), 40)

class BEDViewer(QWidget):
    regions_loaded = pyqtSignal(list)   # emits list of (start, end)

    def __init__(self):
        super().__init__()
        self.setWindowTitle("BED File Viewer")
        self.resize(800, 600)

        # Table widget to display BED data
        self.table = QTableWidget()
        self.table.setColumnCount(0)
        self.table.setRowCount(0)

        # Button to open file
        self.open_button = QPushButton("Open BED File")
        self.open_button.clicked.connect(self.open_file)

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.open_button)
        layout.addWidget(self.table)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)
    
    def load_bed(self, file_path):
        with open(file_path, "r", encoding="utf-8") as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

        if not lines:
            QMessageBox.warning(self, "Empty File", "The BED file is empty or contains only comments.")
            return

        # Split lines into columns
        data = [line.split("\t") for line in lines]

        # Set table dimensions
        max_cols = max(len(row) for row in data)
        self.table.setColumnCount(max_cols)
        self.table.setRowCount(len(data))

        # Optional: Set headers for first 3 BED columns
        headers = ["chrom", "start", "end"] + [f"col{i+4}" for i in range(max_cols - 3)]
        self.table.setHorizontalHeaderLabels(headers)

        # Fill table
        for row_idx, row_data in enumerate(data):
            for col_idx, value in enumerate(row_data):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)  # Make read-only
                self.table.setItem(row_idx, col_idx, item)

        regions = []
        for row in data:
            if len(row) >= 3:
                try:
                    start = int(row[1])
                    end = int(row[2])
                    regions.append((start, end))
                except:
                    pass

        self.regions_loaded.emit(regions)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BED Chromosome Viewer")

        self.bed_viewer = BEDViewer()
        self.chrom_widget = ChromosomeWidget()

        # Connect BED viewer → chromosome widget
        self.bed_viewer.regions_loaded.connect(self.chrom_widget.set_regions)

        layout = QHBoxLayout()
        layout.addWidget(self.bed_viewer, 2)
        layout.addWidget(self.chrom_widget, 3)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
