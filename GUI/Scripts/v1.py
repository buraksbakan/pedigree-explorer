import sys
import csv
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout, QWidget, QMessageBox
)
import pyqtgraph as pg


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
    sys.exit(app.exec())
