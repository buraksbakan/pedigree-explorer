import sys
from PyQt6.QtWidgets import QApplication, QLabel, QPushButton, QVBoxLayout, QWidget, QFileDialog
from PyQt6.QtGui import QPixmap, QImage
from PyQt6.QtCore import Qt


class ImageExporter(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PyQt6 Image Export Example")
        self.resize(400, 300)

        # QLabel to display the image
        self.image_label = QLabel(alignment=Qt.AlignmentFlag.AlignCenter)

        # Buttons
        self.load_button = QPushButton("Load Image")
        self.save_button = QPushButton("Export Image")
        self.save_button.setEnabled(False)  # Disabled until an image is loaded

        # Layout
        layout = QVBoxLayout()
        layout.addWidget(self.image_label)
        layout.addWidget(self.load_button)
        layout.addWidget(self.save_button)
        self.setLayout(layout)

        # Connect signals
        self.load_button.clicked.connect(self.load_image)
        self.save_button.clicked.connect(self.export_image)

        self.current_pixmap = None

    def load_image(self):
        """Load an image from disk into QLabel."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Image", "", "Images (*.png *.jpg *.jpeg *.bmp *.gif)"
        )
        if file_path:
            pixmap = QPixmap(file_path)
            if not pixmap.isNull():
                self.current_pixmap = pixmap
                self.image_label.setPixmap(pixmap.scaled(
                    self.image_label.size(), Qt.AspectRatioMode.KeepAspectRatio
                ))
                self.save_button.setEnabled(True)
            else:
                print("Failed to load image.")

    def export_image(self):
        """Export the current image to a chosen file."""
        if self.current_pixmap:
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save Image", "", "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg)"
            )
            if file_path:
                # Save as QPixmap
                if not self.current_pixmap.save(file_path):
                    print("Failed to save image.")
                else:
                    print(f"Image saved to {file_path}")


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ImageExporter()
    window.show()
    sys.exit(app.exec())
