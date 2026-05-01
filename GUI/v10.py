import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout, QComboBox, 
    QColorDialog, QToolBar, QWidgetAction
)
from PyQt6.QtCore import Qt, pyqtSignal
import pyqtgraph as pg
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go

#from plotly import graph_objs as go

# plotting chromosomes
class ChromosomePlot(QWidget):
    cleared = pyqtSignal()

    def __init__(self):
        super().__init__()

        self.chrom_dict = {}
        
        self.view = QWebEngineView()
        layout = QVBoxLayout()
        layout.addWidget(self.view)
        self.setLayout(layout)

        # COLOUR!

        self.current_highlight_color = "red"
        self.used_highlight_colors = set()

        self.highlight_colors = [
            ("red", "darkred"),
            ("orange", "darkorange"),
            ("green", "darkgreen"),
            ("purple", "indigo"),
            ("cyan", "darkcyan")
        ]
        
        self.active_annotation_filter = None   # None = no colouring

        #self.use_annotation_colors = False
        self.annotation_colors = {}
        self.color_index = 0

        self.highlight_index = 0
        self.highlight_shapes = []  # store shape IDs so we can clear them later

    def save_image(self, file_path):
        if hasattr(self, "fig") and file_path:
            self.fig.write_image(file_path)

    def set_regions(self, chrom_dict):
        """
        Receive list of (start, end) tuples from BEDViewer.
        """
        self.chrom_dict = chrom_dict
        self.update_plot()

    def show_all_chromosomes(self):
        self.single_chrom_mode = False
        self.update_plot()

    def update_plot(self):
        if not self.chrom_dict:
            return

        fig = self._plot_all_chromosomes()
        self.fig = fig

        html = fig.to_html(include_plotlyjs="cdn", full_html=False)
        self.view.setHtml(html)


    def _plot_all_chromosomes(self):
        fig = go.Figure()

        # vertical layout: each chromosome at a different x
        x_offset = 0
        spacing = 3

        max_length = max(
            max(end for start, end, *_ in regions)
            for regions in self.chrom_dict.values()
        )

        for chrom, regions in sorted(self.chrom_dict.items()):
            chrom_length = max(end for start, end, *_ in regions)

            # backbone
            fig.add_shape(
                type="rect",
                x0=x_offset - 0.4, x1=x_offset + 0.4,
                y0=0, y1=chrom_length,
                line=dict(color="black", width=2),
                fillcolor="lightgray"
            )

            # regions
            for start, end, annotation in regions:

                # Decide whether to colour this region
                if self.active_annotation_filter is None:
                    color = "royalblue"  # default
                elif self.active_annotation_filter == "ALL":
                    color = self.get_color_for_annotation(annotation)
                elif annotation == self.active_annotation_filter:
                    color = self.get_color_for_annotation(annotation)
                else:
                    color = "lightgray"  # dim unselected annotations

                fig.add_shape(
                    type="rect",
                    x0=x_offset - 0.4, x1=x_offset + 0.4,
                    y0=start, y1=end,
                    fillcolor=color,
                    opacity=0.6,
                    line=dict(width=0)
                )

            # label as annotation
            fig.add_annotation(
                x=x_offset,
                y=-max_length * 0.05,
                text=chrom,
                showarrow=False,
                yanchor="top"
            )

            x_offset += spacing

        # legend 
        legend_entries = []

        # Annotation legend
        if self.active_annotation_filter == "ALL":
            for annotation, color in self.annotation_colors.items():
                legend_entries.append((f"Annotation: {annotation}", color))

        elif isinstance(self.active_annotation_filter, str):
            color = self.annotation_colors.get(self.active_annotation_filter)
            if color:
                legend_entries.append((f"Annotation: {self.active_annotation_filter}", color))

        # Add annotation legend traces
        for name, color in legend_entries:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=12, color=color),
                name=name
            ))

        # Add highlight legend traces (only used colours)
        for fill in self.used_highlight_colors:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],
                mode="markers",
                marker=dict(size=12, color=fill),
                name=f"Highlight: {fill}"
            ))


        fig.update_xaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[-spacing, x_offset]
        )

        fig.update_yaxes(
            title="Genomic position (bp)",
            range=[-max_length * 0.1, max_length]
        )

        fig.update_layout(
            height=600,
            margin=dict(l=40, r=20, t=40, b=40),
            showlegend=True,
            title="All chromosomes",
            legend=dict(
                title="Colour Legend",
                orientation="v",
                x=1.02,
                y=1,
                bordercolor="black",
                borderwidth=1
            )
        )
        return fig

    
    def highlight_region(self, chrom, start, end):
        # Determine x-offset for this chromosome
        chroms = sorted(self.chrom_dict.keys())
        index = chroms.index(chrom)
        x_offset = index * 3  # same spacing used in _plot_all_chromosomes()

        # Use the currently selected highlight colour
        fill = self.current_highlight_color
        line = "black"

        # Track that this colour is used (for legend)
        self.used_highlight_colors.add(fill)

        # Add highlight rectangle
        shape = dict(
            type="rect",
            x0=x_offset - 0.4, x1=x_offset + 0.4,
            y0=start, y1=end,
            fillcolor=fill,
            opacity=0.6,
            line=dict(color=line, width=2)
        )

        self.fig.add_shape(shape)
        self.highlight_shapes.append(shape)

        # Re-render
        html = self.fig.to_html(include_plotlyjs="cdn", full_html=False)
        self.view.setHtml(html)
    
    def get_color_for_annotation(self, annotation):
        if annotation not in self.annotation_colors:
            fill, line = self.highlight_colors[self.color_index % len(self.highlight_colors)]
            self.annotation_colors[annotation] = fill
            self.color_index += 1
        return self.annotation_colors[annotation]
    
    def clear_highlights(self):
        self.highlight_shapes = []
        self.used_highlight_colors = set()
        self.current_highlight_color = "red"

        self.annotation_colors = {}
        self.color_index = 0
        self.active_annotation_filter = None

        self.update_plot()
        self.cleared.emit() # notify main window the highlight was cleared




# loading BED file and presenting it in a table
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(dict)
    region_selected = pyqtSignal(str, int, int)

    def __init__(self):
        super().__init__()

        # Table to show BED contents
        self.table = QTableWidget()

        # Button to open BED file
        self.open_button = QPushButton("Open BED File")
        self.open_button.clicked.connect(self.open_file)

        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.cellClicked.connect(self.row_clicked)
        
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

                annotation = row[3] if len(row) > 3 else None
                chrom_dict[chrom].append((start, end, annotation))


        # Emit regions to chromosome plot
        self.regions_loaded.emit(chrom_dict)
    
    def row_clicked(self, row, col):
        chrom = self.table.item(row, 0).text()
        start = int(self.table.item(row, 1).text())
        end   = int(self.table.item(row, 2).text())
        self.region_selected.emit(chrom, start, end)


# using the BED file to create the visualisation
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BED Chromosome Viewer")
        self.resize(1200, 600)

        # table view (left) and plot (right)
        self.bed_viewer = BEDViewer()
        self.chrom_plot = ChromosomePlot()

        # connect bed view with plot
            # send regions to the plot
        self.bed_viewer.regions_loaded.connect(self.chrom_plot.set_regions)
            # populate annotation dropdown
        self.bed_viewer.regions_loaded.connect(self.populate_chrom_list)
            # highlight region when clicked on the table row
        self.bed_viewer.region_selected.connect(self.chrom_plot.highlight_region)


        # menu bar
        menu = self.menuBar()
        file_menu = menu.addMenu("File")
        view_menu = menu.addMenu("View")

            # export
        file_menu.addAction("Export Image", self.export_image)
        file_menu.addAction("Quit", self.close)

            ## when highlights clear the annotation dropdown is reset
        view_menu.addAction("Clear Highlights", self.chrom_plot.clear_highlights)
            ## make click highlights + choose colour
        view_menu.addAction("Pick Highlight Colour", self.pick_highlight_color)
        
        # # export button 
        # self.save_button = QPushButton("Export Image")
        # self.save_button.clicked.connect(self.export_image)

        # # make click highlights + choose colour
        # self.pick_color_button = QPushButton("Pick Highlight Colour")
        # self.pick_color_button.clicked.connect(self.pick_highlight_color)
        
        # # clear both annotation and click highlights
        # self.clear_button = QPushButton("Clear Highlights")
        # self.clear_button.clicked.connect(self.chrom_plot.clear_highlights)
        #     # when highlights clear the annotation dropdown is reset
        # self.chrom_plot.cleared.connect(self.reset_annotation_dropdown)

        # annotations to colour bed regions
        self.annotation_dropdown = QComboBox()
        self.annotation_dropdown.addItem("No annotation colouring")
        self.annotation_dropdown.addItem("All annotations")
        self.annotation_dropdown.currentTextChanged.connect(self.annotation_selection_changed)
        

        
        
        # layout
            # plot layout is vertical so plot on top buttons on the bottom
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.chrom_plot)

        # add the buttons
        controls_layout = QHBoxLayout()
        #controls_layout.addWidget(self.save_button)
        #controls_layout.addWidget(self.clear_button)
        controls_layout.addWidget(self.annotation_dropdown)
        #controls_layout.addWidget(self.pick_color_button)
        
        # combine the layouts
        plot_layout.addLayout(controls_layout)

        # main layout
        layout = QHBoxLayout()
        layout.addWidget(self.bed_viewer, 2) # add table

        plot_container = QWidget()
        plot_container.setLayout(plot_layout)
        layout.addWidget(plot_container, 3) # add plot and buttons

        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    
    def pick_highlight_color(self):
        # highlight colour picker
        color = QColorDialog.getColor()
        if color.isValid():
            self.chrom_plot.current_highlight_color = color.name()

    def reset_annotation_dropdown(self):
        # reset dropdown when highlights are cleared
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.setCurrentText("No annotation colouring")
        self.annotation_dropdown.blockSignals(False)

    def annotation_selection_changed(self, text):
        # annotation dropdown changes
        # clear all highlights and legend
        self.chrom_plot.clear_highlights()

        # set the new annotation mode
        if text == "No annotation colouring":
            self.chrom_plot.active_annotation_filter = None
        elif text == "All annotations":
            self.chrom_plot.active_annotation_filter = "ALL"
        else:
            self.chrom_plot.active_annotation_filter = text

        # redraw with the new mode
        self.chrom_plot.update_plot()

    def populate_chrom_list(self, chrom_dict):
        # populate the dropdown box with the stuff in annotation column
        self.chrom_list = chrom_dict

        # extract unique annotations
        annotations = set()
        for regions in chrom_dict.values():
            for _, _, annotation in regions:
                if annotation:
                    annotations.add(annotation)

        # fill dropdown
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("No annotation colouring")
        self.annotation_dropdown.addItem("All annotations")
        for a in sorted(annotations):
            self.annotation_dropdown.addItem(a)


    def chromosome_selected(self, chrom_name):
        if chrom_name == "All chromosomes":
            self.chrom_plot.show_all_chromosomes()
        else:
            return

    def export_image(self):   
        # save image
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Chromosome Plot",
            "",
            "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg);;SVG Files (*.svg)"
        )

        if not file_path:
            return  # user cancelled

        try:
            self.chrom_plot.save_image(file_path)
            print(f"Image saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")



#run application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
