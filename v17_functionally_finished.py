import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QPushButton, QFileDialog, QVBoxLayout,
    QWidget, QMessageBox, QTableWidget, QTableWidgetItem, QHBoxLayout, QComboBox, 
    QColorDialog, QStatusBar, QInputDialog, QDialog, QCheckBox
)
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
from collections import defaultdict
from PyQt6.QtWidgets import QAbstractItemView
import pkgutil

# Ideogram imports
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Rectangle
import pyideogram

data = pkgutil.get_data("IBD_Viewer", "resources/chromosome_bed/chr1.bed")

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    if hasattr(sys, '_MEIPASS'):
        # Running inside PyInstaller bundle
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

#from plotly import graph_objs as go

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

# Plain colours: no banding, centromere pinch only
_PLAIN_BAND_COLORS = {
    "gneg": (0.92, 0.92, 0.92), "gpos25": (0.92, 0.92, 0.92),
    "gpos50": (0.92, 0.92, 0.92), "gpos75": (0.92, 0.92, 0.92),
    "gpos100": (0.92, 0.92, 0.92), "acen": (0.75, 0.25, 0.25),
    "gvar": (0.92, 0.92, 0.92), "stalk": (0.92, 0.92, 0.92),
}

# Annotation colour palette (matches ChromosomePlot's highlight_colors order)
_IDEO_PALETTE = [
    "#E74C3C", "#FF8C00", "#27AE60", "#8E44AD", "#00BCD4",
    "#F39C12", "#2E86C1", "#D35400", "#1ABC9C", "#C0392B",
]


class IdeogramPlot(QWidget):
    """
    pyideogram-based chromosome ideogram.  View-only — no click interaction.
    Reads colour state from the parent MainWindow via set_colour_state().
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
        self.colour_by_column = None
        self.annotation_colors = {}
        self.color_index = 0

        self.gene_index = {}   # chrom → list of (start, end, name)
        self.show_all_genes = False
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
        self.all_regions = [(chrom_dict, headers)]
        self.headers = headers
        self.update_plot()

    def update_plot(self):
        #print("UPDATE PLOT JUST RAN")
        self.figure.clear()
        if not self.all_regions:
            self.canvas.draw()
            return
        self._draw_overview()
        self.canvas.draw()


    def add_highlight(self, chrom, start, end, annotation_value=None):
        #print("ADD HIGHLIGHT on", id(self))
        #print("ADD HIGHLIGHT", chrom, start, end, annotation_value)  # TEMP
        self.highlights.append((chrom, start, end, annotation_value))
        self.update_plot()


    def clear_highlights(self):
        self.highlights = []
        self.highlight_colors = {}
        self.current_highlight_column = None
        self.current_annotation_filter = None
        self.update_plot()
        self.cleared.emit()

    def save_image(self, file_path):
        if hasattr(self, "fig") and file_path:
            self.fig.write_image(file_path)

    def _get_color_for_annotation(self, value):
        if value not in self.annotation_colors:
            self.annotation_colors[value] = _IDEO_PALETTE[
                self.color_index % len(_IDEO_PALETTE)
            ]
            self.color_index += 1
        return self.annotation_colors[value]

    def _draw_overview(self):
        # Collect chromosomes from data
        #print("DRAW OVERVIEW HAPPENS NOW")
        chroms = set()
        for chrom_dict, _ in self.all_regions:
            chroms.update(chrom_dict.keys())
        chroms = sorted(chroms)

        if not chroms:
            return

        n = len(chroms)
        gs = self.figure.add_gridspec(
            nrows=2, ncols=n,
            height_ratios=[20, 1],
            wspace=0.08, hspace=0.02,
            left=0.02, right=0.90, top=0.92, bottom=0.04,
        )

        margin = _MAX_CHROM_LEN * 0.01

        # Reset legend entries each draw
        legend_entries = {}


        for i, chrom in enumerate(chroms):
            #print("THIS IS THE FOR I CHROM IN ENUMERATE BIT")
            ax = self.figure.add_subplot(gs[0, i])

            # Draw ideogram backbone
            if chrom in _CHROM_LENGTHS_BP:
                chrom_length = _CHROM_LENGTHS_BP[chrom]
                pyideogram.ideogramv(chrom, ax=ax, color=_PLAIN_BAND_COLORS)
            else:
                chrom_length = 50_000_000
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
                    if self.colour_by_column and self.colour_by_column in headers:
                        col_idx = headers.index(self.colour_by_column)
                        annotation_value = row[col_idx] if col_idx < len(row) else None

                    # CASE 1 — No colouring selected → neutral blue
                    if self.colour_by_column is None:
                        #print("CASE 1")
                        region_color = "#4A90D9"

                    # CASE 2 — Colour-by active
                    else:
                        #print("CASE 2")
                        if annotation_value is not None:
                            # Assign palette colour if new
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

            genes = self.get_genes_in_region(chrom, 0, chrom_length)
            self.draw_gene_labels(ax, chrom, 0, chrom_length, genes)

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

                # Highlight colour (independent of region colours)
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

        # Legend for annotation colours
        if legend_entries:
            #print("LEGEND ENTRIES")
            ax_leg = self.figure.add_axes([0.91, 0.40, 0.08, 0.50])
            ax_leg.axis("off")
            patches = [
                Patch(facecolor=c, alpha=0.6, label=v, edgecolor="gray", linewidth=0.5)
                for v, c in legend_entries.items()
            ]
            ax_leg.legend(
                handles=patches,
                loc="upper left",
                fontsize=7,
                framealpha=0.9,
                title="Annotations",
                title_fontsize=8,
            )

        self.figure.suptitle(
            "Ideogram View  (pyideogram / GRCh38)",
            fontsize=12,
            fontweight="bold",
            y=0.97,
        )
        
    def get_genes_in_region(self, chrom, start, end):
        results = []
        for gstart, gend, gname in self.gene_index.get(chrom, []):
            if not (gend < start or gstart > end):
                results.append((gname, gstart, gend))
        return results
    
    def draw_gene_labels(self, ax, chrom, start, end, genes):
        # Decide which genes to draw
        if self.show_all_genes:
            visible_genes = genes
        else:
            visible_genes = [g for g in genes if g[0] in self.selected_genes]

        if not visible_genes:
            return
                
        # Compute midpoints
        gene_positions = []
        for gname, gstart, gend in genes:
            mid = (gstart + gend) / 2
            gene_positions.append([mid, gname, gstart, gend])

        # Sort by midpoint (top to bottom)
        gene_positions.sort(key=lambda x: x[0])

        # Apply vertical staggering
        min_gap = (end - start) * 0.02   # 2% of region height
        last_y = None

        for i in range(len(gene_positions)):
            y, gname, gstart, gend = gene_positions[i]

            if last_y is not None and abs(y - last_y) < min_gap:
                y = last_y + min_gap

            gene_positions[i][0] = y
            last_y = y


        if self.show_genes_checkbox.isChecked():
            for y, gname, gstart, gend in gene_positions:
                ax.text(
                    0.85, y,
                    gname,
                    fontsize=9,
                    va="center",
                    ha="left",
                    color="black",
                    zorder=9999,
                )
                # Gene bar
                self.fig.add_shape(
                    type="line",
                    x0=0.3, x1=0.7,
                    y0=y, y1=y,
                    line=dict(color="black", width=2)
                )






class RegionPopup(QDialog):
    def __init__(self, chrom, start, end, genes, intervals=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"{chrom}:{start}-{end}")
        self.intervals = intervals or []
        self.interval_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"][:len(self.intervals)]

        self.chrom = chrom
        self.start = start
        self.end = end
        self.genes = genes
        self.region_colour = "royalblue"

        layout = QHBoxLayout(self)

        # TABLE
        self.table = QTableWidget()
        self.table.setColumnCount(3)
        self.table.setHorizontalHeaderLabels(["Gene", "Start", "End"])
        self.table.setRowCount(len(genes))

        for i, (gname, gstart, gend) in enumerate(genes):
            self.table.setItem(i, 0, QTableWidgetItem(gname))
            self.table.setItem(i, 1, QTableWidgetItem(str(gstart)))
            self.table.setItem(i, 2, QTableWidgetItem(str(gend)))

        # PLOT VIEW
        self.view = QWebEngineView()

        # COLOUR BUTTON
        self.color_button = QPushButton("Change Region Colour")
        self.color_button.clicked.connect(self.pick_region_color)
        self.interval_color_button = QPushButton("Change Interval Colours")
        self.interval_color_button.clicked.connect(self.pick_interval_colors)

        # SAVE BUTTON
        self.save_button = QPushButton("Export Image")
        self.save_button.clicked.connect(self.export_image)

        # SHOW/UNSHOW GENE LABELS
        self.show_genes_checkbox = QCheckBox("Show gene midpoints")
        self.show_genes_checkbox.setChecked(True)
        self.show_genes_checkbox.stateChanged.connect(self.on_show_genes_checkbox)
        self.show_all_genes = True
        self.selected_genes = set()
        self.table.cellClicked.connect(self.on_gene_table_clicked)


        # BUTTONS
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(self.color_button)
        controls_layout.addWidget(self.save_button)
        controls_layout.addWidget(self.interval_color_button)

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
        color = QColorDialog.getColor()
        if color.isValid():
            self.region_colour = color.name()
            self.update_plot()

    def pick_interval_colors(self):
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
            fillcolor=self.region_colour,
            opacity=0.5,
            line=dict(width=2)
        )

        # Distinct colours for intervals
        interval_colors = ["#1f77b4", "#2ca02c", "#ff7f0e", "#d62728", "#9467bd"]

        lane_width = 0.15
        lane_spacing = 0.05

        for idx, (label, istart, iend) in enumerate(self.intervals):
            color = self.interval_colors[idx]

            # lane positioning
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

            fig.add_annotation(
                x= x1 + 0.02,
                y=(istart + iend) / 2,
                text=label,
                showarrow=False,
                font=dict(size=10, color=color),
                xanchor="left"
            )

            # fade above if interval extends upward
            if istart < start:
                fig.add_shape(
                    type="rect",
                    x0=x0, x1=x1,
                    y0=ymin, y1=start,
                    fillcolor=color,
                    opacity=0.15,
                    line=dict(width=0)
                )

            # fade below if interval extends downward
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
        for idx, (label, istart, iend) in enumerate(self.intervals):
            color = self.interval_colors[idx]
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
            visible_genes = [g for g in self.genes if g[0] in self.selected_genes]

        for gname, gstart, gend in visible_genes:
            mid = (gstart + gend) / 2

            fig.add_shape(
                type="line",
                x0=0.82, x1=0.92,
                y0=mid, y1=mid,
                line=dict(color="black", width=2)
            )

            fig.add_annotation(
                x=0.95,
                y=mid,
                text=gname,
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
            width=300,
            height=1000,
            margin=dict(l=20, r=20, t=20, b=20)
        )

        html = fig.to_html(include_plotlyjs="cdn", full_html=False)
        self.view.setHtml(html)
        self.fig = fig


    def on_show_genes_checkbox(self, state):
        self.show_all_genes = bool(state)
        self.update_plot()

    def on_gene_table_clicked(self, row, column):
        gene_name = self.table.item(row, 0).text()

        if gene_name in self.selected_genes:
            self.selected_genes.remove(gene_name)
        else:
            self.selected_genes.add(gene_name)

        # If checkbox is OFF, only selected genes will show
        # If checkbox is ON, all genes show anyway
        self.update_plot()


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
            self.fig.save_image(file_path)
            #print(f"Image saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")




# loading BED file and presenting it in a table
class BEDViewer(QWidget):
    # Signal emitted when BED regions are loaded
    regions_loaded = pyqtSignal(dict, list, list)
    region_selected = pyqtSignal(str, int, int)
    three_way_overlaps_found = pyqtSignal(list)


    def __init__(self):
        super().__init__()

        # Table to show BED contents
        self.table = QTableWidget()
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.cellClicked.connect(self.row_clicked)
        
        # Layout
        layout = QVBoxLayout()
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
        #print("BED FILE LOADED")
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

        # Determine number of columns
        max_cols = max(len(row) for row in data)

        # Build headers (generic but consistent)
        
        headers = ["chrom", "start", "end"] + [f"field{i}" for i in range(4, max_cols+1)]
        self.headers = headers

        # Pad all rows so they match header length
        for row in data:
            while len(row) < len(headers):
                row.append(None)

        # Fill table
        self.table.setColumnCount(len(headers))
        self.table.setRowCount(len(data))
        self.table.setHorizontalHeaderLabels(headers)

        for r, row in enumerate(data):
            for c, value in enumerate(row):
                item = QTableWidgetItem("" if value is None else str(value))
                item.setFlags(item.flags() ^ Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)

        # Build chrom_dict for plotting
        # Extract BED regions
        chrom_dict = {}   # { "chr1": [ { "start":..., "end":..., "columns":[...] }, ... ], ... }


        for row in data:
            chrom = row[0].strip().lower().replace("chr", "")
            chrom = f"chr{chrom}"

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

        # Group regions by chromosome
        by_chrom = defaultdict(list)  # chrom -> list of (start, end, annot)

        for row in data:
            if len(row) >= 4:
                chrom = row[0].strip().lower().replace("chr", "")
                chrom = f"chr{chrom}"

                try:
                    start = int(row[1])
                    end = int(row[2])
                except ValueError:
                    continue
                annot = row[3]
                by_chrom[chrom].append((start, end, annot))

        three_way_overlaps = []

        for chrom, regions in by_chrom.items():
            n = len(regions)
            # brute-force triple combinations: fine for typical BED sizes
            for i in range(n):
                for j in range(i + 1, n):
                    for k in range(j + 1, n):
                        (s1, e1, a1) = regions[i]
                        (s2, e2, a2) = regions[j]
                        (s3, e3, a3) = regions[k]

                        # need three distinct annotations (e.g. a-b, b-c, a-c)
                        if len({a1, a2, a3}) < 3:
                            continue

                        overlap_start = max(s1, s2, s3)
                        overlap_end   = min(e1, e2, e3)

                        if overlap_start < overlap_end:
                            interval_list = [
                                (a1, s1, e1),
                                (a2, s2, e2),
                                (a3, s3, e3)
                            ]

                            three_way_overlaps.append(
                                (chrom, overlap_start, overlap_end, interval_list)
                            )


        # deduplicate
        # dedupe by (chrom, start, end) but keep interval lists
        unique = {}
        for chrom, s, e, intervals in three_way_overlaps:
            key = (chrom, s, e)
            if key not in unique:
                unique[key] = intervals

        self.three_way_overlaps = [
            (chrom, s, e, unique[(chrom, s, e)])
            for (chrom, s, e) in unique
        ]


        # Emit regions to chromosome plot
        self.regions_loaded.emit(chrom_dict, self.headers, data)

        # emit overlaps to mainwindow
        self.three_way_overlaps_found.emit(self.three_way_overlaps)

    
    def row_clicked(self, row, col):
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
        super().__init__()
        self.setWindowTitle("BED Chromosome Viewer")
        self.resize(1200, 600)

        # core widgets
        # table view (left) and plot (right)
        self.bed_viewer = BEDViewer()
        self.ideogram_plot = IdeogramPlot()

        # annotations to colour bed regions
        self.annotation_dropdown = QComboBox()
        self.annotation_dropdown.addItem("No annotation colouring")
        self.annotation_dropdown.addItem("All annotations")
        self.annotation_dropdown.currentTextChanged.connect(self.annotation_selection_changed)

        
        self.colour_by_dropdown = QComboBox()
        self.colour_by_dropdown.addItem("No colouring")
        self.colour_by_dropdown.currentTextChanged.connect(self.colour_by_changed)

        # # Table for 3-way overlaps
        self.overlap_table = QTableWidget()
        self.overlap_table.setColumnCount(3)
        self.overlap_table.setHorizontalHeaderLabels(["chrom", "start", "end"])
        self.overlap_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.overlap_table.cellClicked.connect(self.overlap_clicked)
        self.overlap_table.hide()


        # Connect signal from BEDViewer
        self.bed_viewer.three_way_overlaps_found.connect(self.update_overlap_table)

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
        # builds main window layout

        # Right Side:
    
            # plot layout is vertical so plot on top buttons on the bottom
        plot_layout = QVBoxLayout()
        plot_layout.addWidget(self.ideogram_plot)
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
        self.left_panel.addWidget(self.overlap_table) 

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
        file_menu.addAction("Import BED File", self.bed_viewer.open_file)
        file_menu.addAction("Export Image", self.export_image)
        file_menu.addAction("Quit", self.close)

        # view menu
        view_menu = menu.addMenu("View")
        #view_menu.addAction("Clear Highlights", self.ideogram_plot.clear_highlights)
        view_menu.addAction("Pick Highlight Colour", self.pick_highlight_color)
        view_menu.addAction("Choose Colour By Column", self.open_colour_by_dialog)



    ######################
    # MENU ACTIONS

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
            self.ideogram_plot.save_image(file_path)
            #print(f"Image saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save image:\n{e}")

    ###################
    # ANNOTATION
    def region_clicked(self, chrom, start, end):
        chrom = chrom.lower().replace("chr", "")
        #print("REGION CLICKED", chrom, start, end)  # TEMP
        self.ideogram_plot.add_highlight(chrom, start, end, None)

    def reset_annotation_dropdown(self):
        #print("RESET ANNOTATION DROPDOWN JUST RAN")
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("No annotation")
        self.annotation_dropdown.blockSignals(False)


    def annotation_selection_changed(self, text):
        #print("ANNOTATION SELECTION CHANGED JUST RAN")

        if text == "No annotation":
            # Turn OFF colour-by mode completely
            self.ideogram_plot.colour_by_column = None
            self.ideogram_plot.current_annotation_filter = None

            # Reset palette so next colour-by starts fresh
            self.ideogram_plot.annotation_colors = {}
            self.ideogram_plot.color_index = 0

        elif text == "All annotations":
            # Turn colour-by mode BACK ON using the last selected column
            if self.ideogram_plot.colour_by_column is None:
                # Restore the last chosen column
                self.ideogram_plot.colour_by_column = self.last_colour_by_column

            self.ideogram_plot.current_annotation_filter = "ALL"

        else:
            # Specific annotation selected
            if self.ideogram_plot.colour_by_column is not None:
                self.ideogram_plot.current_annotation_filter = text

        self.ideogram_plot.update_plot()


    def update_annotation_dropdown(self, reset=True):
        #print("UPDATE ANNOTATION DROPDOWN JUST RAN")
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("No annotation")
        self.annotation_dropdown.addItem("All annotations")

        col = self.ideogram_plot.colour_by_column

        if col and col in self.all_column_values:
            for val in sorted(self.all_column_values[col]):
                self.annotation_dropdown.addItem(val)

        if reset:
            self.annotation_dropdown.setCurrentText("No annotation")

        self.annotation_dropdown.blockSignals(False)





    #####################
    # COLOUR

    def open_colour_by_dialog(self):
        #print("OPEN COLOUR BY DIALOG JUST RAN")
        col, ok = QInputDialog.getItem(
            self,
            "Select Column for Highlighting",
            "Choose a column:",
            sorted(self.all_headers),
            editable=False
        )

        if not ok:
            return

        # THIS is the correct variable to set
        self.ideogram_plot.colour_by_column = col
        self.last_colour_by_column = col

        # Reset annotation colours so they rebuild
        self.ideogram_plot.annotation_colors = {}
        self.ideogram_plot.color_index = 0

        # Populate annotation dropdown
        self.annotation_dropdown.blockSignals(True)
        self.annotation_dropdown.clear()
        self.annotation_dropdown.addItem("No annotation")
        self.annotation_dropdown.addItem("All annotations")

        for val in sorted(self.all_column_values[col]):
            self.annotation_dropdown.addItem(val)

        self.annotation_dropdown.setCurrentText("All annotations")
        self.annotation_dropdown.blockSignals(False)

        # Activate filter
        self.ideogram_plot.current_annotation_filter = "ALL"

        # Redraw
        self.ideogram_plot.update_plot()

        # Show overlap table if needed
        values = self.all_column_values[col]
        if {"A-B", "B-C", "A-C"}.issubset(values):
            self.overlap_table.show()
        else:
            self.overlap_table.hide()



    def colour_by_changed(self, text):
        #print("COLOUR BY CHANGED JUST RAN")
        if text == "No colouring":
            self.ideogram_plot.colour_by_column = None
        else:
            self.ideogram_plot.colour_by_column = text

        self.ideogram_plot.annotation_colors = {}
        self.ideogram_plot.color_index = 0
        self.ideogram_plot.update_plot()

    
    def pick_highlight_color(self):
        # highlight colour picker
        color = QColorDialog.getColor()
        if color.isValid():
            self.ideogram_plot.current_highlight_color = color.name()

    
    def update_colour_by_dropdown(self):
        #print("UPDATE COLOUR BY DROPDOWN JUST RAN")
        self.colour_by_dropdown.blockSignals(True)
        self.colour_by_dropdown.clear()
        self.colour_by_dropdown.addItem("No colouring")
        for h in sorted(self.all_headers):
            self.colour_by_dropdown.addItem(h)
        self.colour_by_dropdown.blockSignals(False)


    #################
    # TABLES



    def handle_regions_loaded(self, chrom_dict, headers, full_rows):
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
        self.update_colour_by_dropdown()
        self.update_annotation_dropdown()

        # Send regions to plot
        self.ideogram_plot.set_regions(chrom_dict, headers)

    def update_overlap_table(self, overlaps):
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

        self.overlaps = normalized


    def overlap_clicked(self, row, col):
        #print("AN OVERLAP WAS CLICKED")
        # When user clicks an overlap row, show in popup.
        chrom, start, end, interval_list = self.overlaps[row]
        
        genes = self.get_genes_in_region(chrom, start, end)

        # determine annotation value from highlight column
        annot = interval_list[0][0] if interval_list else None
        self.ideogram_plot.add_highlight(chrom, start, end, annot)

        popup = RegionPopup(chrom, start, end, genes=genes, intervals=interval_list, parent=self)
        popup.exec()


    def open_region_popup(self, chrom, start, end):
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
        results = []
        for gstart, gend, gname in self.gene_index.get(chrom, []):
            if not (gend < start or gstart > end):
                results.append((gname, gstart, gend))
        return results

    
    def load_gene_annotations(self):
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

        for chrom in chroms:
            path = resource_path(f"resources/chromosome_bed/{chrom}.bed")

            if not os.path.exists(path):
                self.gene_index[chrom] = []
                continue

            genes = []
            with open(path) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        _, start, end, name = parts[:4]
                        genes.append((int(start), int(end), name))

            self.gene_index[chrom] = genes








#run application
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
