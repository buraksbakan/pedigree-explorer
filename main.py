import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout, QComboBox, 
    QColorDialog, QStatusBar, QDialog, QCheckBox, QLabel, QGroupBox, 
    QSizePolicy, QAbstractItemView
)
from PyQt6.QtCore import Qt, pyqtSignal, QTimer
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
from collections import defaultdict
import webbrowser
import re

# Ideogram imports
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Rectangle
import pyideogram
from itertools import combinations


def resource_path(relative_path):
    if hasattr(sys, "_MEIPASS"):
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)


# GRCh38 chromosome lengths (bp)
_CHROM_LENGTHS_BP = {
    "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
    "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
    "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
    "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
    "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
    "chr16": 90338345,  "chr17": 83257441,  "chr18": 80373285,
    "chr19": 58617616,  "chr20": 64444167,  "chr21": 46709983,
    "chr22": 50818468,  "chrX": 156040895,  "chrY": 57227415,
}
_MAX_CHROM_LEN = max(_CHROM_LENGTHS_BP.values())

# Plain colors: no banding, centromere pinch only
_PLAIN_BAND_colorS = {
    "gneg": (0.92, 0.92, 0.92), "gpos25": (0.92, 0.92, 0.92),
    "gpos50": (0.92, 0.92, 0.92), "gpos75": (0.92, 0.92, 0.92),
    "gpos100": (0.92, 0.92, 0.92), "acen": (0.4, 0.4, 0.4),
    "gvar": (0.92, 0.92, 0.92), "stalk": (0.92, 0.92, 0.92),
}

# Annotation color palette (matches ChromosomePlot's highlight_colors order)
_IDEO_PALETTE = [
    "#E74C3C", "#FF8C00", "#27AE60", "#8E44AD", "#00BCD4",
    "#F39C12", "#2E86C1", "#D35400", "#1ABC9C", "#C0392B",
]


class IdeogramPlot(QWidget):
    """
    pyideogram-based chromosome ideogram.  View-only — no click interaction.
    Reads color state from the parent MainWindow via set_color_state().
    """
    cleared = pyqtSignal()

    def __init__(self):
        super().__init__()
        #print("IDEO INIT", id(self))

        self.all_regions = [] 

        self.highlights = []  # list of (chrom, start, end, annotation_value)
        self.highlight_colors = {}  # annotation_value -> color
        self.current_highlight_column = None
        self.current_annotation_filter = None
        self.default_highlight_color = "#4A90D9"  # neutral blue
        self.current_highlight_color = "#D97A4A"

        self.headers = []
        self.color_by_column = None
        self.annotation_colors = {}
        self.color_index = 0

        self.gene_index = {}   # chrom → list of (start, end, name)
        self.show_all_genes = False
        self.show_all_cytobands = True
        self.selected_genes = set()


        # Matplotlib figure & canvas
        self.figure = Figure(figsize=(16, 8), dpi=100)
        self.figure.patch.set_facecolor("#FAFAFA")
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)



    # ── public API (called from MainWindow) ────────────────────────────
    def set_regions(self, chrom_dict, headers):
        """Replaces the region list with a new dataset and refreshes the plot"""
        #print("SET_REGIONS CALLED")

        self.all_regions = [(chrom_dict, headers)]
        self.headers = headers
        self.update_plot()

    def update_plot(self):
        """Clears and redraws figure"""
        #print("UPDATE PLOT JUST RAN")

        self.figure.clear()
        if not self.all_regions:
            self.canvas.draw()
            return
        self.draw_overview()
        self.canvas.draw()


    def add_highlight(self, chrom, start, end, annotation_value=None):
        """Adds new highlighted genomic interval when table clicked to highlight list and triggers a plot redraw"""
        #print("ADD HIGHLIGHT", chrom, start, end, annotation_value)

        self.highlights.append((chrom, start, end, annotation_value))
        self.update_plot()

    def save_image(self, file_path):
        """Saves image as a PNG using MatPlotLib"""
        if not file_path:
            return

        # Ensure the figure is freshly drawn
        self.figure.clear()
        self.draw_overview()
        
        #testing
        #print("Saving figure with axes:", len(self.figure.axes))

        # Save using Matplotlib
        self.figure.savefig(file_path, dpi=300, bbox_inches="tight")

    def draw_overview(self):
        """Draws multi-chromosomal ideogram view with cytobands, highlighted BED regions and annotation legend"""
        #print("DRAW OVERVIEW HAPPENS NOW")

        chroms = set()
        for chrom_dict, _ in self.all_regions:
            chroms.update(chrom_dict.keys())

        # numeric sort for chr1–chr22
        def chrom_key(c):
            return int(c.replace("chr", ""))

        chroms = sorted(chroms, key=chrom_key)


        if not chroms:
            return

        n = len(chroms)
        gs = self.figure.add_gridspec(
            nrows=2, ncols=n,
            height_ratios=[20, 1],
            wspace=0.08, hspace=0.02,
            left=0.02, right=0.82, top=0.92, bottom=0.04,
        )

        margin = _MAX_CHROM_LEN * 0.01

        # Reset legend entries each draw
        legend_entries = {}


        for i, chrom in enumerate(chroms):
            #print("THIS IS THE FOR I CHROM IN ENUMERATE BIT")
            ax = self.figure.add_subplot(gs[0, i])


            # Draw ideogram backbone
            # Determine chromosome length once
            chrom_length = _CHROM_LENGTHS_BP.get(chrom, 50_000_000)

            # Draw ideogram backbone
            if chrom in _CHROM_LENGTHS_BP:
                # Draw backbone with or without chromosomes
                if self.show_all_cytobands:
                    pyideogram.ideogramv(chrom, ax=ax, color=_PLAIN_BAND_colorS)
                else:
                    pyideogram.ideogramv(chrom, ax=ax, color=_PLAIN_BAND_colorS, edgecolor = None)
            else:
                # Cytobands OFF or chromosome unknown → draw plain backbone
                ax.add_patch(Rectangle(
                    (-0.4, 0), 0.8, chrom_length,
                    facecolor=(0.92, 0.92, 0.92),
                    edgecolor="black"
                ))

            ax.set_xlim(-0.6, 0.6)
            ax.set_ylim(-margin, _MAX_CHROM_LEN + margin)

            # Draw BED regions
            for chrom_dict, headers in self.all_regions:
                #print("THIS IS THE DRAW BED REGIONS BIT")
                if chrom not in chrom_dict:
                    continue

                for region in chrom_dict[chrom]:
                    start = region["start"]
                    end = region["end"]
                    row = region["columns"]

                    # Determine annotation value
                    annotation_value = None
                    if self.color_by_column and self.color_by_column in headers:
                        col_idx = headers.index(self.color_by_column)
                        annotation_value = row[col_idx] if col_idx < len(row) else None

                    # CASE 1 — No coloring selected → neutral blue
                    if self.color_by_column is None:
                        #print("CASE 1")
                        region_color = "#4A90D9"

                    # CASE 2 — color-by active
                    else:
                        #print("CASE 2")
                        if annotation_value is not None:
                            # Assign palette color if new
                            if annotation_value not in self.annotation_colors:
                                self.annotation_colors[annotation_value] = _IDEO_PALETTE[
                                    self.color_index % len(_IDEO_PALETTE)
                                ]
                                self.color_index += 1

                            region_color = self.annotation_colors[annotation_value]
                            legend_entries[annotation_value] = region_color
                        else:
                            # Unannotated regions stay neutral blue
                            region_color = "#4A90D9"

                        # Apply annotation filter: hide non-matching
                        if self.current_annotation_filter not in (None, "ALL", annotation_value):
                            continue

                    # Draw region
                    #print("DRAW REGION")
                    ax.add_patch(Rectangle(
                        (-0.48, start),
                        0.96,
                        end - start,
                        facecolor=region_color,
                        alpha=0.55,
                        edgecolor="none",
                        zorder=5,
                    ))

            # Push ideogram + regions behind highlights
            for artist in ax.get_children():
                try:
                    artist.set_zorder(0)
                except Exception:
                    pass

            # Draw highlights
            for hchrom, hstart, hend, annot in self.highlights:
                #print("DRAW HIGHLIGHT")
                if hchrom != chrom:
                    continue

                # Highlight color (independent of region colors)
                if annot is None:
                    hcolor = self.current_highlight_color
                else:
                    if annot not in self.highlight_colors:
                        self.highlight_colors[annot] = _IDEO_PALETTE[
                            len(self.highlight_colors) % len(_IDEO_PALETTE)
                        ]
                    hcolor = self.highlight_colors[annot]

                # Apply annotation filter to highlights too
                if self.current_annotation_filter not in (None, "ALL", annot):
                    continue


                ax.add_patch(Rectangle(
                    (-0.55, hstart),
                    1.10,
                    hend - hstart,
                    facecolor="none",
                    edgecolor=hcolor,
                    linewidth=6.0,
                    zorder=9999,
                ))

            # Chromosome label
            ax_label = self.figure.add_subplot(gs[1, i])
            label_text = chrom.replace("chr", "")
            ax_label.text(
                0.5, 0.5, label_text,
                ha="center", va="center",
                fontsize=7, fontweight="bold"
            )
            ax_label.set_xlim(0, 1)
            ax_label.set_ylim(0, 1)
            ax_label.axis("off")

            ax.set_ylabel("")
            ax.set_yticks([])
            ax.set_xticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)

            

        # legend_entries already built during region drawing
        # but now remove entries that are filtered out
        if self.current_annotation_filter not in (None, "ALL"):
            legend_entries = {
                k: v for k, v in legend_entries.items()
                if k == self.current_annotation_filter
            }

        # Legend for annotation colors
        if legend_entries:
            patches = [
                Patch(facecolor=c, alpha=0.6, label=v, edgecolor="gray", linewidth=0.5)
                for v, c in legend_entries.items()
            ]

            # Create a dedicated legend axis on the right side
            ax_leg = self.figure.add_axes([0.85, 0.15, 0.07, 0.70])
            ax_leg.axis("off")

            ax_leg.legend(
                handles=patches,
                loc="upper left",
                fontsize=7,
                framealpha=0.9,
                title="Comparisons",
                title_fontsize=8,
            )



        self.figure.suptitle(
            "Ideogram View  (pyideogram / GRCh38)",
            fontsize=12,
            fontweight="bold",
            y=0.97,
        )






class RegionPopup(QDialog):
    def __init__(self, chrom, start, end, genes, intervals=None, parent=None):
        """Initalizes the popup window, builds the gene table, 
        interval controls, and renders the plot"""
        #print("Genes passed to popup:", genes)

        super().__init__(parent)
        self.setWindowTitle(f"{chrom}:{start}-{end}")
        self.intervals = intervals or []
        self.interval_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"][:len(self.intervals)]

        self.chrom = chrom
        self.start = start
        self.end = end
        # Clean gene list: ensure start/end are valid integers
        clean_genes = []
        for entry in genes:
            # Detect format:
            # If first element is a string and NOT an Ensembl ID → it's the gene name
            if isinstance(entry[0], str) and not entry[0].startswith("ENSG"):
                # Format: (name, start, end, id)
                n, s, e, i = entry
            else:
                # Format: (start, end, name, id)
                s, e, n, i = entry

            try:
                s = int(s)
                e = int(e)
                clean_genes.append((s, e, n, i))
            except Exception:
                continue

        self.genes = clean_genes



        self.region_color = "royalblue"

        layout = QHBoxLayout(self)

        # TABLE
        self.table = QTableWidget()
        self.table.setColumnCount(4)
        self.table.setHorizontalHeaderLabels(["Gene", "Start", "End", "Ensembl ID"])
        self.table.setRowCount(len(self.genes))

        for i, (gstart, gend, gname, stable_id) in enumerate(self.genes):
            display_name = gname if gname else "-"
            self.table.setItem(i, 0, QTableWidgetItem(display_name))
            self.table.setItem(i, 1, QTableWidgetItem(str(gstart)))
            self.table.setItem(i, 2, QTableWidgetItem(str(gend)))
            self.table.setItem(i, 3, QTableWidgetItem(stable_id))



        # PLOT VIEW
        self.view = QWebEngineView()

        # color BUTTON
        self.color_button = QPushButton("Change Region Colour")
        self.color_button.clicked.connect(self.pick_region_color)
        self.interval_color_button = QPushButton("Colour of three-way comparisons")
        self.interval_color_button.clicked.connect(self.pick_interval_colors)

        # SAVE BUTTON
        self.export_gene_table_button = QPushButton("Save Gene Table")
        self.export_gene_table_button.clicked.connect(self.export_gene_table)

        # SHOW/UNSHOW GENE LABELS
        self.show_genes_checkbox = QCheckBox("Show all genes")
        self.show_genes_checkbox.setChecked(True)
        self.show_genes_checkbox.stateChanged.connect(self.on_show_genes_checkbox)
        self.show_all_genes = True
        self.selected_genes = set()
        self.table.cellClicked.connect(self.on_gene_table_clicked)


        # BUTTONS
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.color_button)
        controls_layout.addWidget(self.interval_color_button)
        controls_layout.addWidget(self.export_gene_table_button)

        # RIGHT SIDE LAYOUT
        right_side = QVBoxLayout()
        right_side.addLayout(controls_layout)
        right_side.addWidget(self.show_genes_checkbox)
        right_side.addWidget(self.view)
        
        right_side.addStretch()

        layout.addWidget(self.table, stretch=2)
        layout.addLayout(right_side, stretch=3)
        
        self.update_plot()

        

    def pick_region_color(self):
        """Opens a colour picker to reassign the region's 
        highlight colour on the plot"""
        color = QColorDialog.getColor()
        if color.isValid():
            self.region_color = color.name()
            self.update_plot()

    def pick_interval_colors(self):
        """Opens a colour picker to reassign the 
        interval colours on the plot"""
        if not self.intervals:
            return
        
        new_colors = []
        for label, istart, iend in self.intervals:
            color = QColorDialog.getColor()
            if color.isValid():
                new_colors.append(color.name())
            else:
                new_colors.append("#1f77b4")
        self.interval_colors = new_colors
        self.update_plot()

    def update_plot(self):
        """Rebuilds plot of the selected genomic region, 
        including the region box, intervals, ticks, and gene labels."""
        fig = go.Figure()
        start = self.start
        end = self.end
        genes = self.genes

        # Add padding above/below region
        pad = int((end - start) * 0.15)
        ymin = max(0, start - pad)
        ymax = end + pad
        
        # Invisible trace to force axis scaling
        fig.add_trace(go.Scatter(
            x=[0.5, 0.5],
            y=[ymin, ymax],
            mode="markers",
            marker=dict(opacity=0),
            showlegend=False
        ))
        
        # Grey fade above region
        fig.add_shape(
            type="rect",
            x0=0.2, x1=0.8,
            y0=ymin, y1=start,
            fillcolor="lightgrey",
            opacity=0.3,
            line=dict(width=0)
        )
        # Grey fade below region
        fig.add_shape(
            type="rect",
            x0=0.2, x1=0.8,
            y0=end, y1=ymax,
            fillcolor="lightgrey",
            opacity=0.3,
            line=dict(width=0)
        )

        # Draw vertical region
        fig.add_shape(
            type="rect",
            x0=0.2, x1=0.8,
            y0=start, y1=end,
            fillcolor=self.region_color,
            opacity=0.5,
            line=dict(width=2)
        )

        # Distinct colors for intervals
        interval_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]
        lane_width = 0.15
        lane_spacing = 0.05

        for idx, item in enumerate(self.intervals):
            # Accept both (label, start, end) and (label, (start, end))
            if len(item) == 3 and not isinstance(item[1], tuple):
                label, istart, iend = item
            else:
                label, (istart, iend) = item

            color = interval_colors[idx % len(interval_colors)]

            x0 = 0.2 + idx * (lane_width + lane_spacing)
            x1 = x0 + lane_width

            fig.add_shape(
                type="rect",
                x0=x0, x1=x1,
                y0=istart, y1=iend,
                fillcolor=color,
                opacity=0.4,
                line=dict(width=1, color=color)
            )

            if istart < start:
                fig.add_shape(
                    type="rect",
                    x0=x0, x1=x1,
                    y0=ymin, y1=start,
                    fillcolor=color,
                    opacity=0.15,
                    line=dict(width=0)
                )

            if iend > end:
                fig.add_shape(
                    type="rect",
                    x0=x0, x1=x1,
                    y0=end, y1=ymax,
                    fillcolor=color,
                    opacity=0.15,
                    line=dict(width=0)
                )


        legend_items = []
        for idx, item in enumerate(self.intervals):
            if len(item) == 3 and not isinstance(item[1], tuple):
                label, istart, iend = item
            else:
                label, (istart, iend) = item

            color = interval_colors[idx % len(interval_colors)]
            legend_items.append(
                dict(name=label, marker=dict(color=color), mode="lines")
            )


        for item in legend_items:
            fig.add_trace(go.Scatter(
                x=[None], y=[None],  # invisible
                mode="markers",
                marker=dict(size=10, color=item["marker"]["color"]),
                name=item["name"]
            ))


        # Tick marks
        for pos in [start, end]:
            fig.add_shape(type="line", x0=0.15, x1=0.80, y0=pos , y1=pos,
                        line=dict(color="black", width=2))
            fig.add_annotation(x=0.95, y=pos, text=f"{pos:,}",
                            showarrow=False, font=dict(size=10))
            


        # GENE LABELS
        if self.show_all_genes:
            visible_genes = self.genes
        else:
            visible_genes = [g for g in self.genes if g[2] in self.selected_genes]


        for gstart, gend, gname, gid in visible_genes:
            mid = (gstart + gend) / 2
            label = gname if gname else gid

            fig.add_shape(
                type="line",
                x0=0.82, x1=0.92,
                y0=mid, y1=mid,
                line=dict(color="black", width=2)
            )

            fig.add_annotation(
                x=0.95,
                y=mid,
                text=label,
                showarrow=False,
                font=dict(size=10),
                xanchor="left"
            )


        # Vertical axis = genomic coordinates
        fig.update_yaxes(
            range=[ymin, ymax],  
            autorange=False,
            title="Position (bp)"
        )
        fig.update_xaxes(visible=False)
        
        fig.update_layout(
            autosize=False,
            width=450,
            height=1000,
            margin=dict(l=10, r=10, t=10, b=10)
        )

        html = fig.to_html(include_plotlyjs="cdn", full_html=False)
        self.view.setHtml(html)
        self.fig = fig


    def on_show_genes_checkbox(self, state):
        """Toggles whether all or select genes are shown on the plot"""
        self.show_all_genes = bool(state)
        self.update_plot()

    def on_gene_table_clicked(self, row, column):
        """Adds or removes a clicked gene from the selected genes and updates the plot"""
        gene_name = self.table.item(row, 0).text()
        if not gene_name:
            gene_name = self.table.item(row, 3).text()


        if gene_name in self.selected_genes:
            self.selected_genes.remove(gene_name)
        else:
            self.selected_genes.add(gene_name)

        # If checkbox is OFF, only selected genes will show
        # If checkbox is ON, all genes show anyway
        self.update_plot()

    def export_gene_table(self):
        """Exports the gene table for the region as a CSV file"""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Gene Table",
            "",
            "CSV Files (*.csv)"
        )

        if not file_path:
            return

        try:
            with open(file_path, "w", encoding="utf-8") as f:
                # Write header
                headers = [self.table.horizontalHeaderItem(i).text()
                        for i in range(self.table.columnCount())]
                f.write(",".join(headers) + "\n")

                # Write rows
                for row in range(self.table.rowCount()):
                    values = []
                    for col in range(self.table.columnCount()):
                        item = self.table.item(row, col)
                        values.append(item.text() if item else "")
                    f.write(",".join(values) + "\n")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save table:\n{e}")






# loading BED file and presenting it in a table
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(dict, list, list, str)
    region_selected = pyqtSignal(str, int, int)
    three_way_overlaps_found = pyqtSignal(list)


    def __init__(self, parent = None):
        """Initalizes the BED viewer table and creates layout"""
        super().__init__(parent)

        # Table to show BED contents
        self.table = QTableWidget()
        self.table.setMinimumSize(0, 0)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.cellClicked.connect(self.row_clicked)
        
        # Layout
        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.table)
        self.setLayout(layout)


    def open_file(self):
        """Open file dialog and load BED file."""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open File",
            "",
            "BED Files (*.bed);;All Files (*)"
        )
        if file_path:
            self.current_filename = os.path.basename(file_path)
            try:
                self.load_bed(file_path)

            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to read file:\n{e}")

    def load_bed(self, file_path):
        """Read BED file, compute three-way overlaps, populate table, and emit regions."""

        filename = os.path.basename(file_path)

        with open(file_path, "r", encoding="utf-8") as f:
            lines = [
                line.strip()
                for line in f
                if line.strip() and not line.startswith("#")
            ]

        if not lines:
            QMessageBox.warning(self, "Empty File", "The file is empty.")
            return

        # Split into columns
        data = [line.split("\t") for line in lines]

        # Determine number of columns
        max_cols = max(len(row) for row in data)

        # Build headers
        headers = ["chrom", "start", "end", "comparisons"] + [
            f"field{i}" for i in range(5, max_cols + 1)
        ]
        self.headers = headers

        # Pad rows
        for row in data:
            while len(row) < len(headers):
                row.append(None)

        # Fill table
        self.table.setColumnCount(len(headers))
        self.table.setRowCount(len(data))
        self.table.setHorizontalHeaderLabels(headers)

        for col in range(4, len(headers)):
            self.table.setColumnHidden(col, True)

        for r, row in enumerate(data):
            for c, value in enumerate(row):
                item = QTableWidgetItem("" if value is None else str(value))
                item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)

        # Build chrom_dict for plotting
        chrom_dict = {}

        for row in data:
            # chrom = row[0].strip().lower().replace("chr", "")
            # chrom = f"chr{chrom}"
            chrom = self.normalize_chrom(row[0])

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

        # Pre-populate autosomes
        for i in range(1, 23):
            chrom_dict.setdefault(f"chr{i}", [])

        # Group regions by chromosome
        from collections import defaultdict
        by_chrom = defaultdict(list)

        for row in data:
            if len(row) >= 4:
                # chrom = row[0].strip().lower().replace("chr", "")
                # chrom = f"chr{chrom}"
                chrom = self.normalize_chrom(row[0])

                try:
                    start = int(row[1])
                    end = int(row[2])
                except ValueError:
                    continue

                annot = row[3]
                by_chrom[chrom].append((start, end, annot))

        # -------------------------------
        # ROBUST 3-WAY OVERLAP DETECTION
        # -------------------------------

        from itertools import combinations

        def merge_intervals(intervals):
            if not intervals:
                return []
            intervals = sorted(intervals)
            merged = [intervals[0]]
            for s, e in intervals[1:]:
                ls, le = merged[-1]
                if s <= le:
                    merged[-1] = (ls, max(le, e))
                else:
                    merged.append((s, e))
            return merged

        three_way_overlaps = []

        for chrom, regions in by_chrom.items():
            # Group by annotation
            annot_groups = defaultdict(list)
            for start, end, annot in regions:
                annot_groups[annot].append((start, end))

            # Merge intervals within each annotation
            merged = {
                annot: merge_intervals(iv_list)
                for annot, iv_list in annot_groups.items()
            }

            # All combinations of 3 annotation labels
            for a1, a2, a3 in combinations(merged.keys(), 3):
                iv1 = merged[a1]
                iv2 = merged[a2]
                iv3 = merged[a3]

                # Brute-force over merged intervals
                for s1, e1 in iv1:
                    for s2, e2 in iv2:
                        for s3, e3 in iv3:
                            overlap_start = max(s1, s2, s3)
                            overlap_end   = min(e1, e2, e3)
                            if overlap_start < overlap_end:
                                three_way_overlaps.append(
                                    (chrom, overlap_start, overlap_end,
                                    [(a1, s1, e1), (a2, s2, e2), (a3, s3, e3)])
                                )

        # Deduplicate by (chrom, start, end)
        unique = {}
        for chrom, s, e, intervals in three_way_overlaps:
            key = (chrom, s, e)
            if key not in unique:
                unique[key] = intervals

        self.three_way_overlaps = [
            (chrom, s, e, unique[(chrom, s, e)])
            for (chrom, s, e) in unique
        ]

        # Build all_column_values
        self.all_column_values = defaultdict(set)
        for row in data:
            for col_index, value in enumerate(row):
                header = self.headers[col_index]
                if value is not None:
                    self.all_column_values[header].add(value)
        
        for entry in self.three_way_overlaps:
            print(entry)


        # Emit to GUI
        self.regions_loaded.emit(chrom_dict, self.headers, data, filename)
        self.three_way_overlaps_found.emit(self.three_way_overlaps)


    def normalize_chrom(self, s):
        s = s.strip()
        s = s.lower()
        s = re.sub(r"^chr", "", s)
        s = re.sub(r"[^\w]", "", s)   # remove hidden unicode
        return f"chr{s}"

            
    
    def row_clicked(self, row, col):
        """"""
        #print("A ROW WAS CLICKED")
        chrom = self.table.item(row, 0).text()
        chrom = chrom.strip().lower().replace("chr", "")
        chrom = f"chr{chrom}"

        start = int(self.table.item(row, 1).text())
        end   = int(self.table.item(row, 2).text())

        self.region_selected.emit(chrom, start, end)





# using the BED file to create the visualisation
class MainWindow(QMainWindow):
    def __init__(self):
        """Initalizes the main window, creates widgets, and connects actions to functions"""
        super().__init__()
        self.setWindowTitle("IBDogram Viewer")
        self.resize(1200, 600)

        self.current_filename = "-"

        # core widgets
        # table view (left) and plot (right)
        self.bed_viewer = BEDViewer(self)
        self.ideogram_plot = IdeogramPlot()

        # annotations to color bed regions
        self.annotation_dropdown = QComboBox()
        self.annotation_dropdown.addItem("Total IBD segments")
        self.annotation_dropdown.addItem("All pairwise comparisons")
        self.annotation_dropdown.currentTextChanged.connect(self.annotation_selection_changed)

        # cytobands
        self.show_cytobands_checkbox = QCheckBox("Show Cytobands")
        self.show_cytobands_checkbox.setChecked(True)
        self.show_cytobands_checkbox.stateChanged.connect(self.on_show_cytobands_checkbox)

        # # Table for 3-way overlaps
        self.overlap_table = QTableWidget()
        self.overlap_table.setColumnCount(3)
        self.overlap_table.setHorizontalHeaderLabels(["chrom", "start", "end"])
        self.overlap_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.overlap_table.cellClicked.connect(self.overlap_clicked)
        self.overlap_table.hide()
        self.bed_viewer.three_way_overlaps_found.connect(self.update_overlap_visibility)



        # labels
        self.pairwise_label = QLabel("Pairwise IBD Segments (Click to zoom)")
        self.overlap_label = QLabel("Three-way IBD Overlaps (Click to zoom)")
        self.overlap_label.hide()

        self.summary_box = QGroupBox("File Summary")
        self.summary_layout = QVBoxLayout()
        self.summary_box.setLayout(self.summary_layout)

        self.summary_filename = QLabel("File: -")
        self.summary_rows = QLabel("Number of IBD segments: -")
        self.summary_total_comparisons = QLabel("Comparison(s): -")
        self.summary_each_comparison = QLabel("")

        # Connect signal from BEDViewer
        # overlap table
        self.bed_viewer.three_way_overlaps_found.connect(self.update_overlap_table)
        # pop up for zoom
        self.bed_viewer.region_selected.connect(self.open_region_popup)

        # connect bed view with plot
        self.bed_viewer.regions_loaded.connect(self.handle_regions_loaded)
            # highlight region when clicked on the table row
        self.bed_viewer.region_selected.connect(self.region_clicked)
            # clear highlight resets dropdown
        self.ideogram_plot.cleared.connect(self.reset_annotation_dropdown)

        self.all_headers = set()
        self.all_column_values = {}   # {column_name: set(values)}

        # Build UI
        self.create_layout()
        self.create_menu()

        # Status bar
        self.setStatusBar(QStatusBar(self))

        #load gene files
        self.gene_index = {}  # chrom → list of (start, end, name)
        self.ideogram_plot.gene_index = self.gene_index
        self.load_gene_annotations()
        

    def create_layout(self):
        """Builds the main window layout"""
        # builds main window layout

        self.summary_layout.addWidget(self.summary_filename)
        self.summary_layout.addWidget(self.summary_rows)
        self.summary_layout.addWidget(self.summary_total_comparisons)
        self.summary_layout.addWidget(self.summary_each_comparison)

        # Right Side:
    
        # plot layout is vertical so plot on top buttons on the bottom
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.ideogram_plot)
            ## add the buttons
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.annotation_dropdown)
        controls_layout.addWidget(self.show_cytobands_checkbox)
        

        ## combine the layouts
        plot_layout.addLayout(controls_layout)
            ## create plot container
        plot_container = QWidget()
        plot_container.setLayout(plot_layout)


        # Left Side:
        
        self.left_panel = QVBoxLayout()
        self.left_panel.addWidget(self.pairwise_label)
        self.left_panel.addWidget(self.bed_viewer, stretch = 2)
        self.left_panel.addWidget(self.overlap_label)
        self.left_panel.addWidget(self.overlap_table, stretch = 3) 

        self.overlap_table.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding
        )

        left_container = QWidget()
        left_container.setLayout(self.left_panel)
        left_container.setSizePolicy(QSizePolicy.Policy.Expanding,
                             QSizePolicy.Policy.Expanding)
        self.bed_viewer.setSizePolicy(QSizePolicy.Policy.Expanding,
                              QSizePolicy.Policy.Expanding)


        # Main horizontal layout (left + right)
        content_layout = QHBoxLayout()
        content_layout.addWidget(left_container, 2)
        content_layout.addWidget(plot_container, 3)

        # Wrap everything in a vertical layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.summary_box)
        main_layout.addLayout(content_layout)

        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

    def create_menu(self):
        # creates the menu bar and actions

        menu = self.menuBar()

        # file menu
        file_menu = menu.addMenu("File")
        file_menu.addAction("Open File", self.bed_viewer.open_file)
        file_menu.addAction("Save Image", self.export_image)
        file_menu.addAction("Quit", self.close)

        # help menu
        help_menu = menu.addMenu("Help")
        help_menu.addAction("GitHub", self.open_webbrowser)



    ######################
    # MENU ACTIONS

    # Source - https://stackoverflow.com/a/42845714
    # Posted by Achayan
    # Retrieved 2026-04-20, License - CC BY-SA 3.0

    def open_webbrowser(self):
        """Opens the project's GitHub repository in the users default browser"""
        # if you want please remember to check its valid url or not
        url = "https://github.com/amarshall1312/pedigree-explorer"
        webbrowser.open(url, new=0, autoraise=True) 


    def export_image(self):  
        """Opens a save file dialog and exports the ideogram plot to PNG""" 
        # save image
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save Ideogram",
            "",
            "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg);;SVG Files (*.svg)"
        )
        

        if not file_path:
            return  # user cancelled

        try:
            self.ideogram_plot.save_image(file_path)
            #print(f"Image saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")

    def handle_regions_loaded(self, chrom_dict, headers, full_rows, filename):
        """Processes a loaded BED file, building value maps, updating dropdowns, 
        updating the summary field and sending regions to the ideogram to plot"""
        #print("HANDLE REGIONS LOADED RAN, filename =", self.current_filename)

        self.current_filename = filename
        self.headers = headers[:]   # store ordered list
        
        # Store headers
        for h in headers:
            self.all_headers.add(h)

        # Build column → values map
        self.all_column_values = {h: set() for h in headers}

        for row in full_rows:
            for i, value in enumerate(row):
                self.all_column_values[headers[i]].add(value)

        # Update dropdowns
        # self.update_color_by_dropdown()
        self.update_annotation_dropdown(reset=False)

        # Send regions to plot
        self.ideogram_plot.set_regions(chrom_dict, headers)
        
        # color by comparisons
        self.apply_color_by_column("comparisons")

        # Update summary box
        self.summary_filename.setText(f"File: {self.current_filename}")

        # Number of rows
        self.summary_rows.setText(f"Number of IBD segments: {len(full_rows)}")

        # Comparisons column
        comparisons = sorted(self.all_column_values["comparisons"])

        self.summary_total_comparisons.setText(
            f"Comparison(s): {len(comparisons)} ({', '.join(comparisons)})"
        )

        # Count each comparison
        counts = []
        for comp in comparisons:
            count = sum(1 for row in full_rows if row[3] == comp)
            counts.append(f"{comp}: {count}")

        self.summary_each_comparison.setText("   ".join(counts))
        QTimer.singleShot(0, lambda: self.resize_table_columns(self.bed_viewer.table))


    ######## COMPARISON AND ANNOTATIONS ###########
    def region_clicked(self, chrom, start, end):
        """Highlights the clicked BED region on the ideogram"""
        chrom = chrom.lower().replace("chr", "")
        #print("REGION CLICKED", chrom, start, end)  # TEMP
        self.ideogram_plot.add_highlight(chrom, start, end, None)

    def reset_annotation_dropdown(self):
        """Resets comparison filter dropdown"""
        #print("RESET ANNOTATION DROPDOWN JUST RAN")
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("Total IBD segments")
        self.annotation_dropdown.blockSignals(False)


    def annotation_selection_changed(self, text):
        """Updates comparison filter dropdown and colours by user selection"""
        #print("ANNOTATION SELECTION CHANGED JUST RAN")

        if text == "Total IBD segments":
            # Turn OFF color-by mode completely
            self.ideogram_plot.color_by_column = None
            self.ideogram_plot.current_annotation_filter = None

            # Reset palette so next color-by starts fresh
            self.ideogram_plot.annotation_colors = {}
            self.ideogram_plot.color_index = 0

        elif text == "All pairwise comparisons":
            # Turn color-by mode BACK ON using the last selected column
            if self.ideogram_plot.color_by_column is None:
                # Restore the last chosen column
                self.ideogram_plot.color_by_column = self.last_color_by_column

            self.ideogram_plot.current_annotation_filter = "ALL"

        else:
            # Specific annotation selected
            if self.ideogram_plot.color_by_column is not None:
                self.ideogram_plot.current_annotation_filter = text

        self.ideogram_plot.update_plot()


    def update_annotation_dropdown(self, reset=True):
        """Updates the comparison filter dropdown with all comparison values upon loadout"""
        #print("UPDATE ANNOTATION DROPDOWN JUST RAN")

        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("Total IBD segments")
        self.annotation_dropdown.addItem("All pairwise comparisons")

        col = self.ideogram_plot.color_by_column

        if col and col in self.all_column_values:
            for val in sorted(self.all_column_values[col]):
                self.annotation_dropdown.addItem(val)

        if reset:
            self.annotation_dropdown.setCurrentText("Total IBD segments")

        self.annotation_dropdown.blockSignals(False)

    ######### COLOUR ############

    # previously open color by dialog, but dialog isnt needed now
    def apply_color_by_column(self, col):
        """Colours by comparison value, rebuilds comparison dropdown,
        redraws the ideogram with the highlight """

        self.ideogram_plot.color_by_column = col
        self.last_color_by_column = col

        # Reset annotation colors so they rebuild
        self.ideogram_plot.annotation_colors = {}
        self.ideogram_plot.color_index = 0

        self.update_annotation_dropdown(reset=False)

        # Populate annotation dropdown
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("Total IBD segments")
        self.annotation_dropdown.addItem("All pairwise comparisons")

        for val in sorted(self.all_column_values[col]):
            self.annotation_dropdown.addItem(val)

        # Set filter before drawing
        self.ideogram_plot.current_annotation_filter = "ALL"

        # Draw first
        self.ideogram_plot.update_plot()

        # Now force dropdown selection AFTER plot redraw

        self.annotation_dropdown.setCurrentText("All pairwise comparisons")
        self.annotation_dropdown.blockSignals(False)



    ####### CYTOBANDS #########

    def on_show_cytobands_checkbox(self, state):
        """Updates the plot to show or hide cytobands"""
        self.ideogram_plot.show_all_cytobands = bool(state)        
        self.ideogram_plot.update_plot()

    ####### OVERLAPS ########
    
    def update_overlap_table(self, overlaps):
        """Populates the three-way overlap table"""
        normalized = []
        self.overlap_table.setRowCount(len(overlaps))

        for r, item in enumerate(overlaps):
            if len(item) == 3:
                chrom, start, end = item
                interval_list = []   # fallback
            else:
                chrom, start, end, interval_list = item

            # store normalized 4‑tuple
            normalized.append((chrom, start, end, interval_list))

            # fill table
            self.overlap_table.setItem(r, 0, QTableWidgetItem(chrom))
            self.overlap_table.setItem(r, 1, QTableWidgetItem(str(start)))
            self.overlap_table.setItem(r, 2, QTableWidgetItem(str(end)))

        self.resize_table_columns(self.overlap_table)
        #self.resize_table_columns(self.bed_viewer.table)
        self.overlaps = normalized

    def update_overlap_visibility(self, overlaps):
        """Show or hide the overlap table based on whether overlaps exist."""
        if overlaps:
            self.overlap_label.show()
            self.overlap_table.show()
        else:
            self.overlap_label.hide()
            self.overlap_table.hide()

    def overlap_clicked(self, row, col):
        """Opens a region popup for the selected overlap"""
        #print("AN OVERLAP WAS CLICKED")

        # When user clicks an overlap row, show in popup.
        chrom, start, end, interval_list = self.overlaps[row]
        
        genes = self.get_genes_in_region(chrom, start, end)

        # determine annotation value from highlight column
        annot = interval_list[0][0] if interval_list else None
        self.ideogram_plot.add_highlight(chrom, start, end, annot)

        popup = RegionPopup(chrom, start, end, genes=genes, intervals=interval_list, parent=self)
        popup.exec()


##################### REGION POPUP ##################

    def open_region_popup(self, chrom, start, end):
        """Creates and displays a popup showing genes and intervals for the selected region"""
        # Get genes in this region
        genes = self.get_genes_in_region(chrom, start, end)

        # Create popup
        popup = RegionPopup(
            chrom,
            start,
            end,
            genes,
            intervals=None,
            parent=self
        )
        popup.exec()

    def get_genes_in_region(self, chrom, start, end):
        """Returns all genes whose coordiantes overlap the region"""
        #print("Get genes in region just ran")

        results = []
        for gstart, gend, gname, gid in self.gene_index.get(chrom, []):
            if not (gend < start or gstart > end):
                results.append((gname, gstart, gend, gid))
        return results

    def load_gene_annotations(self):
        """Loads per-chromosome gene BED files to build a searchable gene index"""
        #print("load_gene_annotations just ran")

        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

        for chrom in chroms:
            genes = []

            path = resource_path(f"resources/ensembl_bed/{chrom}.bed")
            #print("Loading gene file:", path, "Exists:", os.path.exists(path))

            if not os.path.exists(path):
                self.gene_index[chrom] = []
                continue

            try:
                with open(path, "r", encoding="utf-8") as f:
                    for line in f:
                        parts = line.strip().split("\t")

                        if len(parts) < 5:
                            continue
                        if not parts[1].isdigit() or not parts[2].isdigit():
                            continue

                        start = int(parts[1])
                        end = int(parts[2])
                        gene_name = parts[3]
                        stable_id = parts[4]

                        genes.append((start, end, gene_name, stable_id))

            except Exception as e:
                print(f"Error reading {chrom}: {e}")
                genes = [] # fallback

            self.gene_index[chrom] = genes


    def resizeEvent(self, event):
        """Resizes table columns and redraws ideogram when the main window is resized"""
        super().resizeEvent(event)

        # BED table
        if self.bed_viewer.table.columnCount() > 0:
            self.resize_table_columns(self.bed_viewer.table)

        # overlap table
        if self.overlap_table.isVisible() and self.overlap_table.columnCount() > 0:
            self.resize_table_columns(self.overlap_table)

        self.ideogram_plot.draw_overview()


    def resize_table_columns(self, table):
        """Evenly distributes table width across all columns"""
        #print("BED viewport width:", self.bed_viewer.table.viewport().width())

        total_width = table.viewport().width()
        
        # count only visible columns
        visible_cols = [i for i in range(table.columnCount()) if not table.isColumnHidden(i)]
        num_cols = len(visible_cols)

        if num_cols == 0:
            return

        col_width = total_width // num_cols
        for i in range(num_cols):
            table.setColumnWidth(i, col_width)

#run application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
