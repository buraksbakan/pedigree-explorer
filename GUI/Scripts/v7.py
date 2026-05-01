import sys
import csv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout
)
from PyQt6.QtCore import Qt, pyqtSignal
import pyqtgraph as pg


# ============================================================
#  ChromosomePlot — Draws chromosome + BED regions using PyQtGraph
# ============================================================
class ChromosomePlot(pg.PlotWidget):
    def __init__(self):
        super().__init__()

        # Store BED regions as (start, end)
        self.regions = []
        self.chrom_length = 1

        # Configure PyQtGraph appearance
        self.setBackground('w')
        self.showGrid(x=True, y=False)
        self.setLabel('bottom', 'Genomic Position (bp)')
        self.setTitle("Chromosome Diagram")

        # Disable auto-scaling so we control the view
        self.enableAutoRange(False)

    def set_regions(self, regions):
        """
        Receive list of (start, end) tuples from BEDViewer.
        """
        self.regions = regions

        if regions:
            # Chromosome length = max end coordinate
            self.chrom_length = max(end for _, end in regions)

        self.update_plot()

    def update_plot(self):
        """
        Draw chromosome backbone + BED intervals.
        """
        self.clear()  # Remove previous drawings

        # Draw chromosome backbone as a horizontal line
        backbone = pg.PlotCurveItem(
            x=[0, self.chrom_length],
            y=[0, 0],
            pen=pg.mkPen(color='black', width=4)
        )
        self.addItem(backbone)

        # Draw each BED region as a colored rectangle
        for start, end in self.regions:
            width = end - start

            rect = pg.BarGraphItem(
                x=[start],
                height=[1],      # height of bar
                width=[width],   # width = region length
                y0=[-0.5],       # center the bar on the backbone
                brush=pg.mkBrush(100, 180, 255, 180)
            )
            self.addItem(rect)

        # Set visible range
        self.setXRange(0, self.chrom_length)
        self.setYRange(-1, 1)


# ============================================================
#  BEDViewer — Loads BED file and displays table
# ============================================================
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(list)

    def __init__(self):
        super().__init__()

        # Table to show BED contents
        self.table = QTableWidget()

        # Button to open BED file
        self.open_button = QPushButton("Open BED File")
        self.open_button.clicked.connect(self.open_file)

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.open_button)
        layout.addWidget(self.table)
        self.setLayout(layout)

    def open_file(self):
        """
        Open file dialog and load BED file.
        """
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open BED File",
            "",
            "BED Files (*.bed);;All Files (*)"
        )
        if file_path:
            try:
                self.load_bed(file_path)
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to read file:\n{e}")

    def load_bed(self, file_path):
        """
        Read BED file, populate table, and emit regions.
        """
        with open(file_path, "r", encoding="utf-8") as f:
            lines = [
                line.strip()
                for line in f
                if line.strip() and not line.startswith("#")
            ]

        if not lines:
            QMessageBox.warning(self, "Empty File", "The BED file is empty.")
            return

        # Split into columns
        data = [line.split("\t") for line in lines]

        # Configure table
        max_cols = max(len(row) for row in data)
        self.table.setColumnCount(max_cols)
        self.table.setRowCount(len(data))

        headers = ["chrom", "start", "end"] + [f"col{i+4}" for i in range(max_cols - 3)]
        self.table.setHorizontalHeaderLabels(headers)

        # Fill table
        for r, row in enumerate(data):
            for c, value in enumerate(row):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)

        # Extract BED regions
        regions = []
        for row in data:
            if len(row) >= 3:
                try:
                    start = int(row[1])
                    end = int(row[2])
                    regions.append((start, end))
                except:
                    pass

        # Emit regions to chromosome plot
        self.regions_loaded.emit(regions)


# ============================================================
#  MainWindow — Combines BEDViewer + ChromosomePlot
# ============================================================
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BED Chromosome Viewer")
        self.resize(1200, 600)

        # Create components
        self.bed_viewer = BEDViewer()
        self.chrom_plot = ChromosomePlot()

        # Connect BED viewer → chromosome plot
        self.bed_viewer.regions_loaded.connect(self.chrom_plot.set_regions)

        # Layout: left = BED table, right = chromosome plot
        layout = QHBoxLayout()
        layout.addWidget(self.bed_viewer, 2)
        layout.addWidget(self.chrom_plot, 3)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)


# ============================================================
#  Run Application
# ============================================================
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
