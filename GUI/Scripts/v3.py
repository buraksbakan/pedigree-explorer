import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QGraphicsView, QGraphicsScene, QGraphicsRectItem, QGraphicsTextItem
)
from PyQt6.QtGui import QBrush, QPen, QColor
from PyQt6.QtCore import Qt, QRectF


class ChromosomeItem(QGraphicsRectItem):
    """Custom QGraphicsRectItem representing a chromosome."""
    def __init__(self, x, y, width, height, label, color=QColor("lightblue")):
        super().__init__(x, y, width, height)
        self.setBrush(QBrush(color))
        self.setPen(QPen(Qt.GlobalColor.black, 1))
        self.setFlag(QGraphicsRectItem.GraphicsItemFlag.ItemIsMovable, True)  # Allow dragging

        # Add label
        text_item = QGraphicsTextItem(label, self)
        text_item.setDefaultTextColor(Qt.GlobalColor.black)
        text_item.setPos(width + 5, height / 4)  # Position label to the right


class ChromosomeViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Chromosome Viewer - PyQt6")
        self.setGeometry(200, 200, 800, 600)

        # Create a QGraphicsScene
        self.scene = QGraphicsScene()
        self.scene.setSceneRect(0, 0, 1000, 800)

        # Create a QGraphicsView
        self.view = QGraphicsView(self.scene, self)
        self.setCentralWidget(self.view)

        # Draw chromosomes
        self.draw_chromosomes()

    def draw_chromosomes(self):
        """Draws a set of chromosomes in the scene."""
        chromosome_data = [
            ("Chr 1", 200),
            ("Chr 2", 180),
            ("Chr 3", 160),
            ("Chr 4", 150),
            ("Chr 5", 140),
        ]

        x_offset = 50
        y_offset = 50
        spacing = 50

        for i, (label, length) in enumerate(chromosome_data):
            # Draw chromosome as a vertical rectangle
            chrom = ChromosomeItem(
                x_offset + i * spacing,
                y_offset,
                20,  # width
                length,  # height
                label
            )
            self.scene.addItem(chrom)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    viewer = ChromosomeViewer()
    viewer.show()
    sys.exit(app.exec())

