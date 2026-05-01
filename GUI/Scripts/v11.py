import sys
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout, QComboBox, 
    QColorDialog, QStatusBar, QInputDialog
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QAction
import pyqtgraph as pg
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go

#from plotly import graph_objs as go

# plotting chromosomes
class ChromosomePlot(QWidget):
    cleared = pyqtSignal()

    def __init__(self):
        super().__init__()

        #self.chrom_dict = {}
        self.fig = go.Figure()
        
        self.view = QWebEngineView()
        layout = QVBoxLayout()
        layout.addWidget(self.view)
        self.setLayout(layout)
        self.all_regions = []   # list of (chrom_dict, headers) for each loaded BED file

        self.colour_by_column = None
        self.chrom_lengths = {
            "chr1": 248956422,
            "chr2": 242193529,
            "chr3": 198295559,
            "chr4": 190214555,
            "chr5": 181538259,
            "chr6": 170805979,
            "chr7": 159345973,
            "chr8": 145138636,
            "chr9": 138394717,
            "chr10": 133797422,
            "chr11": 135086622,
            "chr12": 133275309,
            "chr13": 114364328,
            "chr14": 107043718,
            "chr15": 101991189,
            "chr16": 90338345,
            "chr17": 83257441,
            "chr18": 80373285,
            "chr19": 58617616,
            "chr20": 64444167,
            "chr21": 46709983,
            "chr22": 50818468,
            "chrX": 156040895,
            "chrY": 57227415,
        }


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

    def set_regions(self, chrom_dict, headers):
        # Store this BED file as one track
        self.all_regions.append((chrom_dict, headers))
        self.update_plot()

    def show_all_chromosomes(self):
        self.single_chrom_mode = False
        self.update_plot()

    def update_plot(self):
        if not self.all_regions:
            return

        try:
            fig = self._plot_all_chromosomes()
            self.fig = fig
        except Exception as e:
            print("Plot error:", e)
            return

        html = self.fig.to_html(include_plotlyjs="cdn", full_html=False)
        self.view.setHtml(html)



    def _plot_all_chromosomes(self):
        fig = go.Figure()

        # vertical layout: each chromosome at a different x
        spacing = 3
        chroms = self.get_all_chromosomes()

        if not chroms:
            return fig
        
        # determine max chrom length for y axis
        max_length = max(
            self.chrom_lengths.get(chrom, 0)
            for chrom in chroms
        )

        # draw chromosome backbones
        for i, chrom in enumerate(chroms):
            x_offset = i * spacing
            chrom_length = self.chrom_lengths.get(chrom, 0)

            fig.add_shape(
                type="rect",
                x0=x_offset - 0.4, x1=x_offset + 0.4,
                y0=0, y1=chrom_length,
                line=dict(color="black", width=2),
                fillcolor="lightgray"
            )
            # Label chromosome
            fig.add_annotation(
                x=x_offset,
                y=-max_length * 0.05,
                text=chrom,
                showarrow=False,
                yanchor="top"
            )

        # draw regions form all loaded bed files
        for chrom_dict, headers in self.all_regions:
            for chrom in chroms:

                if chrom not in chrom_dict:
                    continue

                x_offset = chroms.index(chrom) * spacing
                regions = chrom_dict[chrom]

                for region in regions:
                    start = region["start"]
                    end = region["end"]
                    row = region["columns"]

                    # Determine annotation value for colouring
                    if self.colour_by_column is None:
                        annotation_value = None
                    else:
                        col_index = headers.index(self.colour_by_column)
                        annotation_value = row[col_index]

                    # Annotation filtering
                    if self.active_annotation_filter is None:
                        passes_filter = True
                    elif self.active_annotation_filter == "ALL":
                        passes_filter = True
                    else:
                        passes_filter = (annotation_value == self.active_annotation_filter)

                    # Choose colour
                    if not passes_filter:
                        color = "lightgray"
                    else:
                        if annotation_value is None:
                            color = "royalblue"
                        else:
                            color = self.get_color_for_annotation(annotation_value)

                    # Draw region
                    fig.add_shape(
                        type="rect",
                        x0=x_offset - 0.4, x1=x_offset + 0.4,
                        y0=start, y1=end,
                        fillcolor=color,
                        opacity=0.6,
                        line=dict(width=0)
                    )

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


        # axes and layout
        fig.update_xaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            range=[-spacing, (len(chroms) - 1) * spacing + spacing]
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
        chroms = self.get_all_chromosomes()
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

    def get_all_chromosomes(self):
        chroms = set()
        for chrom_dict, headers in self.all_regions:
            chroms.update(chrom_dict.keys())
        return sorted(chroms)
    
    def clear_only_highlights(self):
        self.highlight_shapes = []
        self.used_highlight_colors = set()
        self.current_highlight_color = "red"
        self.update_plot()
        self.cleared.emit()







# loading BED file and presenting it in a table
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(dict, list) # sends chrom: start, end, annot
    region_selected = pyqtSignal(str, int, int) # when row is clicked sends chrom, start, end

    def __init__(self):
        super().__init__()

        # Table to show BED contents
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.cellClicked.connect(self.row_clicked) 
        
        # Layout of table
        layout = QVBoxLayout()
        layout.addWidget(self.table)
        self.setLayout(layout)

    def load_bed(self, file_path):
        # read BED file, populate table, and emit regions.
        
        with open(file_path, "r", encoding="utf-8") as f:
            lines = [
                line.strip()
                for line in f
                if line.strip() and not line.startswith("#")
            ]

        if not lines:
            QMessageBox.warning(self, "Empty File", "The BED file is empty.")
            return

        # split into columns by tab
        data = [line.split("\t") for line in lines]

        # Configure table
        max_cols = max(len(row) for row in data)
        self.table.setColumnCount(max_cols)
        self.table.setRowCount(len(data))

        headers = ["chrom", "start", "end"] + [f"col{i+4}" for i in range(max_cols - 3)]
        self.table.setHorizontalHeaderLabels(headers)
        self.headers = headers


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
                except ValueError:
                    continue

                chrom_dict.setdefault(chrom, []).append({
                    "start": start,
                    "end": end,
                    "columns": row
                })
 
                # if chrom not in chrom_dict:
                #     chrom_dict[chrom] = []

                # annotation = row[3] if len(row) > 3 else None
                # chrom_dict[chrom].append((start, end, annotation))


        # Emit regions to chromosome plot
        self.regions_loaded.emit(chrom_dict, self.headers)
    
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

        # core widgets
        # table view (left) and plot (right)
        self.bed_viewer = BEDViewer()
        self.chrom_plot = ChromosomePlot()

        # annotations to colour bed regions
        self.annotation_dropdown = QComboBox()
        self.annotation_dropdown.addItem("No annotation colouring")
        self.annotation_dropdown.addItem("All annotations")
        self.annotation_dropdown.currentTextChanged.connect(self.annotation_selection_changed)
        
        self.colour_by_dropdown = QComboBox()
        self.colour_by_dropdown.addItem("No colouring")
        self.colour_by_dropdown.currentTextChanged.connect(self.colour_by_changed)

        
        # connect bed view with plot
        self.bed_viewer.regions_loaded.connect(self.handle_regions_loaded)
            # highlight region when clicked on the table row
        self.bed_viewer.region_selected.connect(self.chrom_plot.highlight_region)
            # clear highlight resets dropdown
        self.chrom_plot.cleared.connect(self.reset_annotation_dropdown)

        self.all_headers = set()
        self.all_column_values = {}   # {column_name: set(values)}

        # Build UI
        self.create_layout()
        self.create_menu()

        # Status bar
        self.setStatusBar(QStatusBar(self))

        
    def create_layout(self):
        # builds main window layout

        # Right Side:
    
            # plot layout is vertical so plot on top buttons on the bottom
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.chrom_plot)
            ## add the buttons
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.annotation_dropdown)
            ## combine the layouts
        plot_layout.addLayout(controls_layout)
            ## create plot container
        plot_container = QWidget()
        plot_container.setLayout(plot_layout)

        # Left Side:
        
        self.left_panel = QVBoxLayout()
        self.left_panel.addWidget(self.bed_viewer)

        left_container = QWidget()
        left_container.setLayout(self.left_panel)

        # Main horizontal layout
        main_layout = QHBoxLayout()
        main_layout.addWidget(left_container, 2)
        main_layout.addWidget(plot_container, 3)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def create_menu(self):
        # creates the menu bar and actions

        menu = self.menuBar()

        # file menu
        file_menu = menu.addMenu("File")
        file_menu.addAction("Import BED File", self.load_primary_bed)
        file_menu.addAction("Add Secondary BED File", self.load_secondary_bed)
        file_menu.addAction("Export Image", self.export_image)
        file_menu.addAction("Quit", self.close)

        # view menu
        view_menu = menu.addMenu("View")
        view_menu.addAction("Clear Highlights", self.chrom_plot.clear_only_highlights)
        view_menu.addAction("Pick Highlight Colour", self.pick_highlight_color)
        view_menu.addAction("Choose Colour By Column", self.open_colour_by_dialog)

###### MENU ACTIONS

    def load_primary_bed(self):
        """Loads the main BED file into the first table."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open BED File", "", "BED Files (*.bed);;All Files (*)"
        )
        if file_path:
            self.bed_viewer.load_bed(file_path)

    def load_secondary_bed(self):
        """Loads a second BED file into a new table and links it to the plot."""
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Secondary BED File", "", "BED Files (*.bed);;All Files (*)"
        )
        if not file_path:
            return

        secondary_viewer = BEDViewer()
        secondary_viewer.regions_loaded.connect(self.handle_regions_loaded)
        secondary_viewer.region_selected.connect(self.chrom_plot.highlight_region)
        secondary_viewer.load_bed(file_path)

        # Add to left panel
        self.left_panel.addWidget(secondary_viewer)

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

######### UI HELPERS
    
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
       
        self.chrom_plot.clear_only_highlights()

        self.chrom_plot.highlight_shapes = []
        self.chrom_plot.used_highlight_colors = set()

        # set the new annotation mode
        if text == "No annotation colouring":
            self.chrom_plot.active_annotation_filter = None
            self.chrom_plot.colour_by_column = None
            self.chrom_plot.annotation_colors = {}
            self.chrom_plot.color_index = 0
            self.chrom_plot.update_plot()
            return
        elif text == "All annotations":
            self.chrom_plot.active_annotation_filter = "ALL"
            self.chrom_plot.annotation_colors = {}
            self.chrom_plot.color_index = 0
        else:
            self.chrom_plot.active_annotation_filter = text

        # redraw with the new mode
        self.chrom_plot.update_plot()

    def update_annotation_dropdown(self):
        col = self.colour_by_dropdown.currentText()

        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()

        self.annotation_dropdown.addItem("No annotation colouring")

        if col != "No colouring":
            values = sorted(self.all_column_values.get(col, []))
            self.annotation_dropdown.addItem("All annotations")
            for v in values:
                self.annotation_dropdown.addItem(v)

        self.annotation_dropdown.blockSignals(False)
        

    def colour_by_changed(self, text):
        if text == "No colouring":
            self.chrom_plot.colour_by_column = None
        else:
            self.chrom_plot.colour_by_column = text

        # reset annot columns
        self.chrom_plot.annotation_colors = {}
        self.chrom_plot.color_index = 0

        self.update_annotation_dropdown()
        self.chrom_plot.update_plot()
        
    def update_colour_by_dropdown(self):
        self.colour_by_dropdown.blockSignals(True)
        self.colour_by_dropdown.clear()
        self.colour_by_dropdown.addItem("No colouring")
        for h in sorted(self.all_headers):
            self.colour_by_dropdown.addItem(h)
        self.colour_by_dropdown.blockSignals(False)

    def open_colour_by_dialog(self):
        col, ok = QInputDialog.getItem(
            self,
            "Select Column for Colouring",
            "Choose a column:",
            sorted(self.all_headers),
            editable=False
        )

        if ok:
            self.colour_by_dropdown.setCurrentText(col)
            self.update_annotation_dropdown()
            self.chrom_plot.update_plot()


    def handle_regions_loaded(self, chrom_dict, headers):
        # Add headers to global list
        for h in headers:
            self.all_headers.add(h)

        # Extract values for each column
        for chrom, regions in chrom_dict.items():
            for region in regions:
                row = region["columns"]
                for i, value in enumerate(row):
                    col_name = headers[i]
                    self.all_column_values.setdefault(col_name, set()).add(value)

        # Update dropdowns
        self.update_colour_by_dropdown()
        self.update_annotation_dropdown()

        # Send regions to plot
        self.chrom_plot.set_regions(chrom_dict, headers)








#run application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
