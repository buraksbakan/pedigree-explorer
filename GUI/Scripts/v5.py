import sys
import csv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout, QWidget, QMessageBox, QTableWidget, QTableWidgetItem
)
import pyqtgraph as pg
from PyQt6.QtCore import Qt

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import BasicChromosome
from Bio import SeqIO

# Example: Create a chromosome diagram
def draw_chromosome(output_file="chromosome.pdf"):
    # Create a chromosome diagram
    diagram = BasicChromosome.Organism()
    diagram.page_size = (15 * cm, 4 * cm)  # Width x Height

    # Create a chromosome object
    chromosome = BasicChromosome.Chromosome("Chromosome 1")
    chromosome.scale_num = 1e-6  # Scale in megabases (Mb)

    # Add an "ideogram" (main chromosome body)
    # Here we simulate a chromosome length of 120 Mb
    length = 120_000_000
    centromere_pos = 60_000_000  # Example centromere position

    # Add the left arm
    left_arm = BasicChromosome.ChromosomeSegment(centromere_pos)
    left_arm.fill_color = colors.lightblue
    chromosome.add(left_arm)

    # Add the right arm
    right_arm = BasicChromosome.ChromosomeSegment(length - centromere_pos)
    right_arm.fill_color = colors.lightgreen
    chromosome.add(right_arm)

    # Add the chromosome to the diagram
    diagram.add(chromosome)

    # Draw to PDF
    diagram.draw(output_file, "Chromosome Example")

if __name__ == "__main__":
    try:
        draw_chromosome()
        print("Chromosome diagram saved as 'chromosome.pdf'")
    except Exception as e:
        print(f"Error: {e}")


class BEDViewer(QMainWindow):
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

    def open_file(self):
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


class FilePlotter(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("File Input & Graph Display - PyQt6")
        self.resize(800, 600)

        # Create plot widget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setBackground('w')  # White background
        self.plot_widget.showGrid(x=True, y=True)

        # Create button to load file
        self.load_button = QPushButton("Load Data File")
        self.load_button.clicked.connect(self.load_file)

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.load_button)
        layout.addWidget(self.plot_widget)

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def load_file(self):
        """Open file dialog and load data from CSV/TXT."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open Data File",
            "",
            "Data Files (*.csv *.txt);;All Files (*)"
        )

        if not file_path:
            return  # User cancelled

        try:
            x_data, y_data = self.read_data(file_path)
            if not x_data or not y_data:
                raise ValueError("No valid numeric data found.")

            self.plot_widget.clear()
            self.plot_widget.plot(x_data, y_data, pen=pg.mkPen(color='b', width=2), symbol='o')
            self.plot_widget.setLabel('left', 'Y Axis')
            self.plot_widget.setLabel('bottom', 'X Axis')
            self.plot_widget.setTitle(f"Plot from {file_path.split('/')[-1]}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load file:\n{str(e)}")

    def read_data(self, file_path):
        """Read numeric data from CSV/TXT file."""
        x_data, y_data = [], []
        with open(file_path, 'r', newline='') as file:
            reader = csv.reader(file)
            for row in reader:
                # Skip empty or invalid rows
                if len(row) < 2:
                    continue
                try:
                    x_val = float(row[0].strip())
                    y_val = float(row[1].strip())
                    x_data.append(x_val)
                    y_data.append(y_val)
                except ValueError:
                    continue  # Skip non-numeric rows
        return x_data, y_data


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = FilePlotter()
    window.show()
    viewer = BEDViewer()
    viewer.show()
    sys.exit(app.exec())
