import sys
import csv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout, QComboBox
)
from PyQt6.QtCore import Qt, pyqtSignal
import pyqtgraph as pg


# plotting chromosomes
class ChromosomePlot(pg.PlotWidget):
    def __init__(self):
        super().__init__()

        self.single_chrom_mode = False

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
    
    def set_regions(self, chrom_dict):
        """
        Receive list of (start, end) tuples from BEDViewer.
        """
        self.chrom_dict = chrom_dict
        self.update_plot()

    def update_plot(self):
        
        self.clear()
        if not hasattr(self, "chrom_dict"):
            return

        if self.single_chrom_mode:
            self.plot_single_chromosome()
        else:
            self.plot_all_chromosomes()

    def show_all_chromosomes(self):
        self.single_chrom_mode = False
        self.update_plot()

    def show_single_chromosome(self, chrom_name):
        self.single_chrom_mode = True
        self.selected_chrom = chrom_name
        self.update_plot()

    def plot_all_chromosomes(self):
        """
        Draw chromosome backbone + BED intervals for all chromosomes.
        """
        # Compute the longest chromosome length
        max_length = max(
            max(end for _, end in regions)
            for regions in self.chrom_dict.values()
        )

        y_offset = 0
        spacing = 3  # horizontal spacing between chromosomes

        for chrom, regions in self.chrom_dict.items():

            chrom_length = max(end for _, end in regions)

            # Draw vertical chromosome backbone
            backbone = pg.PlotCurveItem(
                x=[y_offset, y_offset],
                y=[0, chrom_length],
                pen=pg.mkPen(color='black', width=3)
            )
            self.addItem(backbone)

            # Draw regions as vertical bars
            for start, end in regions:
                height = end - start
                rect = pg.BarGraphItem(
                    x=[y_offset - 0.5],   # horizontal position
                    width=[1],            # chromosome thickness
                    y0=[start],           # start coordinate
                    height=[height],      # region length
                    brush=pg.mkBrush(100, 180, 255, 180)
                )
                self.addItem(rect)

            # Add chromosome label
            label = pg.TextItem(chrom, anchor=(0.5, 1))
            label.setPos(y_offset, -200)
            self.addItem(label)

            # Move right for next chromosome
            y_offset += spacing

        # Set visible ranges
        self.setXRange(-2, y_offset + 2)
        self.setYRange(0, max_length)

    def plot_single_chromosome(self):
            # Make sure a chromosome was selected
        chrom_name = getattr(self, "selected_chrom", None)
        if chrom_name is None:
            return

        if chrom_name not in self.chrom_dict:
            return

        regions = self.chrom_dict[chrom_name]
        self.clear()

        chrom_length = max(end for _, end in regions)

        # Draw vertical backbone
        backbone = pg.PlotCurveItem(
            x=[0, 0],
            y=[0, chrom_length],
            pen=pg.mkPen(color='black', width=3)
        )
        self.addItem(backbone)

        # Draw regions
        for start, end in regions:
            height = end - start
            rect = pg.BarGraphItem(
                x=[-0.5],
                width=[1],
                y0=[start],
                height=[height],
                brush=pg.mkBrush(100, 180, 255, 180)
            )
            self.addItem(rect)

        # Label
        label = pg.TextItem(chrom_name, anchor=(0.5, 1))
        label.setPos(0, -200)
        self.addItem(label)

        # Adjust view
        self.setXRange(-2, 2)
        self.setYRange(0, chrom_length)





# loading BED file and presenting it in a table
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(dict)

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
        chrom_dict = {}   # { "chr1": [(start, end), ...], "chr2": [...], ... }

        for row in data:
            if len(row) >= 3:
                chrom = row[0]
                try:
                    start = int(row[1])
                    end = int(row[2])
                except:
                    continue

                if chrom not in chrom_dict:
                    chrom_dict[chrom] = []

                chrom_dict[chrom].append((start, end))


        # Emit regions to chromosome plot
        self.regions_loaded.emit(chrom_dict)



# using the BED file to create the visualisation
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BED Chromosome Viewer")
        self.resize(1200, 600)

        # Create components
        self.bed_viewer = BEDViewer()
        self.chrom_plot = ChromosomePlot()
        self.combobox = QComboBox()
        self.combobox.currentTextChanged.connect(self.chromosome_selected)
        self.combobox.addItem("All chromosomes")

        self.save_button = QPushButton("Export Image") # doesnt do anything atm

        # Connect BED viewer → chromosome plot
        self.bed_viewer.regions_loaded.connect(self.chrom_plot.set_regions)
        self.bed_viewer.regions_loaded.connect(self.populate_chrom_list)

        # # Layout: left = BED table, right = chromosome plot
        # layout = QHBoxLayout()
        # layout.addWidget(self.bed_viewer, 2)
        # layout.addWidget(self.chrom_plot, 3)
        # layout.addWidget(self.save_button)
        # layout.addWidget(self.combobox)

        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.chrom_plot)

        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.save_button)
        controls_layout.addWidget(self.combobox)

        plot_layout.addLayout(controls_layout)

        layout = QHBoxLayout()
        layout.addWidget(self.bed_viewer, 2)

        plot_container = QWidget()
        plot_container.setLayout(plot_layout)
        layout.addWidget(plot_container, 3)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def populate_chrom_list(self, chrom_dict):
        self.chrom_list = chrom_dict
        self.combobox.clear()
        self.combobox.addItem("All chromosomes")  # default option
        self.combobox.addItems(sorted(chrom_dict.keys()))

    def chromosome_selected(self, chrom_name):
        if chrom_name == "All chromosomes":
            self.chrom_plot.show_all_chromosomes()
        else:
            self.chrom_plot.show_single_chromosome(chrom_name)





#run application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
