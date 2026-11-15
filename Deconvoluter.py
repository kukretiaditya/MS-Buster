import tkinter as tk 
from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.ttk import *
import numpy as np 
import os 
import matplotlib.pyplot as plt 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import time 
import threading 
from matplotlib.figure import Figure
import shutil
from matplotlib.widgets import SpanSelector
from multiplierz.mzAPI import mzFile
from copy import deepcopy
from bisect import bisect_left
from scipy.interpolate import interp1d
import mzML
import os
from multiplierz import mzAPI
import traceback
import re
from scipy.signal import find_peaks
from scipy.signal import peak_widths
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, Table, TableStyle, PageBreak
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
import io
import sys

def register():
        print("Trying to Register Interfaces...")
        from multiplierz.mzAPI.management import registerInterfaces
        try:
            print("Registering...")
            registerInterfaces()
        except Exception as e:
            print("Failed Interface Registration:", e)
            print("NOTE: TRY RUNNING AS ADMINISTRATOR")
            pass
def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

class Deconvoluter:
    def __init__(self, root):
        self.root = root 
        self.root.state('zoomed')
        self.root.title("MS-BUSTER")
        try:
            self.root.iconbitmap(default=resource_path("cbt.ico"))
        except:
            pass
        self.menubar = Menu(self.root)

        self.file = Menu(self.menubar, tearoff = 0,relief='ridge')
        self.menubar.add_cascade(label='File', menu=self.file)
        self.file.add_command(label='Load MzML',command=self.load_mzml_file)
        self.file.add_command(label='Load Agilent .d folder', command=self.load_agilent_folder)
        self.file.add_separator()
        self.file.add_command(label='Load AVG Agilent .d folder', command= self.load_ms_from_agilent_folder)
        self.file.add_command(label='Load TXT file',command= self.load_MS_from_txt_file)

        self.tools = Menu(self.menubar, tearoff = 0,relief='ridge')
        self.menubar.add_cascade(label='Tools', menu=self.tools)
        self.tools.add_command(label='Extract and save All MS spectrums between selected time range',command= self.save_ms_spectrum_between_time_range)
        # self.tools.add_command(label='Convert raw .d folder to mzML')

        self.export = Menu(self.menubar, tearoff = 0,relief='ridge')
        self.menubar.add_cascade(label='Export Results', menu=self.export)
        self.export.add_command(label='As PDF',command=self.export_as_pdf)
        self.export.add_separator()
        self.export.add_command(label="Charge table as excel", command=self.export_charge_table_as_excel)

        self.root.config(menu=self.menubar)

        self.gui_start_time = time.time()

        self.file_path = None
        self.importer = None  
        self.mzfile_obj = None 
        self.user_data = None
        self.m_z = None
        self.intensity = None
        self.ms_intensity = None  
        self.tic = None

        self.x_min = None
        self.x_max = None
        self.history = None
        self.ms_spectrum = None
        self.filtered_ms_data = None

        self.time = None
        self.t_min = None
        self.t_max = None

        self.base_name_file = None
        self.output_dir = None
        # self.output_dir_extracted_spectrums = os.path.join(self.output_dir, f"{self.base_name_file}_extracted_spectrum_from {self.t_min} to {self.t_max}")
        # os.makedirs(self.output_dir, exist_ok=True)

        self.spectra_list = []
        self.averaging_in_progress = False
        self.show_peaks_var = tk.BooleanVar(value=True)
    
        self.status_frame = tk.Frame(self.root, bd = 1, relief=tk.RAISED)
        self.status_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.status_var = tk.StringVar()
        self.status_var.set("Status: Ready")

        self.status_bar = tk.Label(self.status_frame, textvariable=self.status_var, bd=1, relief=tk.RAISED, anchor=tk.W)
        self.status_bar.pack(side=tk.LEFT, fill=tk.X)

        # self.clock_var = tk.StringVar()
        # self.clock_label = tk.Label(self.status_frame, textvariable=self.clock_var, bd=1, relief=tk.RAISED, anchor=tk.E)
        # self.clock_label.pack(side=tk.RIGHT, fill=tk.X)

        self.up_time_of_gui_var = tk.StringVar()
        self.up_time_of_gui = tk.Label(self.status_frame, textvariable=self.up_time_of_gui_var)
        self.up_time_of_gui.pack(side=tk.RIGHT, fill = tk.X)

        # self.update_clock()

        self.base_frame = tk.Frame(root)
        self.base_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        v_pane = tk.PanedWindow(self.base_frame, orient=tk.VERTICAL, sashrelief=tk.RAISED)
        v_pane.pack(fill=tk.BOTH, expand=True)

        top_h_pane = tk.PanedWindow(v_pane, orient=tk.HORIZONTAL, sashrelief=tk.RAISED)
        bottom_h_pane = tk.PanedWindow(v_pane, orient=tk.HORIZONTAL, sashrelief=tk.RAISED)

        v_pane.add(top_h_pane)
        v_pane.add(bottom_h_pane)

        self.top_left_frame = tk.Frame(top_h_pane, bd=2, relief=tk.SUNKEN,width=400,height=300)
        self.top_right_frame = tk.Frame(top_h_pane, bd=2, relief=tk.SUNKEN,width=400,height=300)
        self.bottom_left_frame = tk.Frame(bottom_h_pane, bd=2, relief=tk.SUNKEN,width=700,height=300)
        self.bottom_right_frame = tk.Frame(bottom_h_pane, bd=2, relief=tk.SUNKEN,width=400,height=300)

        top_h_pane.add(self.top_left_frame)
        top_h_pane.add(self.top_right_frame)

        bottom_h_pane.add(self.bottom_left_frame)
        bottom_h_pane.add(self.bottom_right_frame)

        # self.step_display = tk.Text(self.top_left_frame,height=9,width=40)
        # self.step_display.pack(pady=3, anchor='e')
        
        # Frame for first text box + scrollbar
        self.text_frame1 = tk.Frame(self.top_left_frame)
        self.text_frame1.pack(fill='both', expand=True, pady=3)

        self.step_display = tk.Text(
            self.text_frame1,
            height=11,
            width=50
        )

        scrollbar1 = tk.Scrollbar(
            self.text_frame1,
            orient="vertical",
            command=self.step_display.yview
        )

        self.step_display.configure(yscrollcommand=scrollbar1.set)

        self.step_display.pack(side=tk.LEFT, fill='both', expand=True)
        scrollbar1.pack(side=tk.RIGHT, fill='y')

        clear_btn = tk.Button(
            self.top_left_frame, 
            text="🗑 Clear", 
            width=6, 
            command=self.clear_step_display
        )
        clear_btn.pack(pady=2, anchor="w")

        # Frame for second text box + scrollbar
        self.text_frame2 = tk.Frame(self.top_left_frame)
        self.text_frame2.pack(fill='both', expand=True, pady=3)

        self.step_display2 = tk.Text(
            self.text_frame2,
            height=7,
            width=40
        )

        scrollbar2 = tk.Scrollbar(
            self.text_frame2,
            orient="vertical",
            command=self.step_display2.yview
        )

        self.step_display2.configure(yscrollcommand=scrollbar2.set)

        self.step_display2.pack(side=tk.LEFT, fill='both', expand=True)
        scrollbar2.pack(side=tk.RIGHT, fill='y')

        self.sub_button_1_frame = Frame(self.top_right_frame)
        self.sub_button_1_frame.pack(anchor="w", pady=5)

        Label(self.sub_button_1_frame, text="Total Ion Chromatogram").pack(side=tk.LEFT, padx=(5, 2))

        self.show_peaks_toggle = Checkbutton(
            self.sub_button_1_frame,
            text="Show Detected Peaks",
            variable=self.show_peaks_var,
            command=self.toggle_peak_display
        )
        self.show_peaks_toggle.pack(side=tk.LEFT, padx=5)
        label = tk.Label(self.sub_button_1_frame, text="Set Time Range:")
        label.pack(side=tk.LEFT, padx=(5, 2)) # 5px padding on left, 2 on right
        self.t_min = tk.Entry(self.sub_button_1_frame, width=10)
        self.t_min.pack(side=tk.LEFT, padx=2)
        self.t_max = tk.Entry(self.sub_button_1_frame, width=10)
        self.t_max.pack(side=tk.LEFT, padx=2)

        self.set_time_button = Button(self.sub_button_1_frame,text="Set time range",command=self.set_time_range)
        self.set_time_button.pack(side = "left", pady=5)

        self.reset_zoom_button = Button(
        self.sub_button_1_frame,
        text="Reset Zoom",
        command=self.reset_time_range
        )

        self.reset_zoom_button.pack(side=tk.LEFT, padx=5)

        self.extract_spectrum = Button(self.sub_button_1_frame, text="Extract AVG MS Spectrum", command=self.data_extractor)
        self.extract_spectrum.pack(side="left", pady=5)

        self.process_data = Combobox(self.sub_button_1_frame, width=27, values=["Centroided", "Default", "Deisotoping"])
        self.process_data.current(1)
        self.process_data.pack(side="left", pady=5)

        self.clear_button = Button(
            self.sub_button_1_frame,
            text="Reset all",
            command=self.clear_all_data
        )

        self.clear_button.pack(side=tk.LEFT, padx=(5, 2))

        self.plot_frame1 = Frame(self.top_right_frame)
        self.plot_frame1.pack(fill='both', expand=True)

        self.fig_right_tab1 = Figure(figsize=(2, 2), dpi=100)
        self.axis_right_tab1 = self.fig_right_tab1.add_subplot(111)
        self.canvas_right_tab1 = FigureCanvasTkAgg(self.fig_right_tab1, master=self.plot_frame1)
        self.canvas_right_tab1.draw()
        self.canvas_right_tab1.get_tk_widget().pack(fill='both', expand=True)
        self.toolbar_right_tab1 = NavigationToolbar2Tk(self.canvas_right_tab1, self.plot_frame1)
        self.toolbar_right_tab1.update()
        self.toolbar_right_tab1.pack()

        # Label(self.bottom_left_frame, text="Mass Spectrum Data").pack(pady=10)

        self.sub_button_02_frame = Frame(self.bottom_left_frame)
        self.sub_button_02_frame.pack(anchor="w", pady=5)

        tk.Label(self.sub_button_02_frame, text="Set m/z Range:").pack(side=tk.LEFT, padx=(5, 2))
        self.x_min = tk.Entry(self.sub_button_02_frame, width=10)
        self.x_min.pack(side=tk.LEFT, padx=(5, 2))
        self.x_max = tk.Entry(self.sub_button_02_frame, width=10)
        self.x_max.pack(side=tk.LEFT, padx=(5, 2))

        self.set_range_button = Button(self.sub_button_02_frame, text="Set m/z range",command=self.set_x_range  )
        self.set_range_button.pack(side=tk.LEFT, padx=(5, 2))

        self.deconvolve_button = Button(self.sub_button_02_frame, text="deconvolve",command=self.Deconvolute)  
        self.deconvolve_button.pack(side=tk.LEFT, padx=(5, 2))

        self.annotatecharge_state = Button(self.sub_button_02_frame,text="Annotate Charge State",command= self.annotate_charge_state)
        self.annotatecharge_state.pack(side=tk.LEFT, padx=(5, 2))
        self.peaks_distance = tk.Entry(self.sub_button_02_frame, width=10)
        self.peaks_distance.pack(side=tk.LEFT, padx=(5, 2))

        self.fig_right_tab2 = Figure(figsize=(6, 4), dpi=100)  # Wider and taller
        self.axis_right_tab2 = self.fig_right_tab2.add_subplot()
        self.canvas_right_tab2 = FigureCanvasTkAgg(self.fig_right_tab2, master=self.bottom_left_frame)
        self.canvas_right_tab2.draw()
        self.canvas_toolbar_frame = tk.Frame(self.bottom_left_frame)
        self.canvas_toolbar_frame.pack(fill='both', expand=True)
        self.toolbar_right_tab2 = NavigationToolbar2Tk(self.canvas_right_tab2, self.canvas_toolbar_frame)
        self.toolbar_right_tab2.update()
        self.toolbar_right_tab2.pack(side=tk.TOP, fill=tk.X)
        self.canvas_right_tab2.get_tk_widget().pack(side=tk.TOP, fill='both', expand=True)

        Label(self.bottom_right_frame, text="Analysis Results").pack(pady=(10, 5))

        table_container = Frame(self.bottom_right_frame)
        table_container.pack(fill='both', expand=True, padx=10, pady=5)
        columns = ('peak1_mz', 'peak2_mz', 'charge1', 'charge2', 'delta', 'mass_kda')
        self.charge_table = Treeview(table_container, columns=columns, show='headings')

        self.charge_table.tag_configure("best_delta", background="#98FB98")   # light green
        self.charge_table.tag_configure("mid_delta",  background="#FFD580")   # light orange/yellow
        self.charge_table.tag_configure("bad_delta",  background="#FFB6C1")   # light pink/red


        self.charge_table.heading('peak1_mz', text='Peak 1 (m/z)')
        self.charge_table.heading('peak2_mz', text='Peak 2 (m/z)')
        self.charge_table.heading('charge1', text='Charge 1')
        self.charge_table.heading('charge2', text='Charge 2')
        self.charge_table.heading('delta', text='Δm')
        self.charge_table.heading('mass_kda', text='Mass (kDa)')

        self.charge_table.column('peak1_mz', width=100, anchor='center')
        self.charge_table.column('peak2_mz', width=100, anchor='center')
        self.charge_table.column('charge1', width=80, anchor='center')
        self.charge_table.column('charge2', width=80, anchor='center')
        self.charge_table.column('delta', width=50,anchor='center')
        self.charge_table.column('mass_kda', width=110, anchor='e') 

        scrollbar = Scrollbar(table_container, orient=VERTICAL, command=self.charge_table.yview)
        self.charge_table.configure(yscrollcommand=scrollbar.set)

        scrollbar.pack(side='right', fill='y')
        self.charge_table.pack(side='left', fill='both', expand=True)
    
    def clear_all_data(self):
        """Completely reset the application state."""

        try:
            self.user_data = None
            self.time = None
            self.intensity = None

            self.ms_spectrum = None
            self.m_z = None
            self.ms_intensity = None
            self.filtered_ms_data = None
            self.charge_analysis_data = []

            self.importer = None
            self.mzfile_obj = None
            self.file_path = None

            self.axis_right_tab1.clear()
            self.axis_right_tab1.set_xlabel("Time")
            self.axis_right_tab1.set_ylabel("Intensity")
            self.axis_right_tab1.set_title("TIC Chromatogram")
            self.canvas_right_tab1.draw()

            self.axis_right_tab2.clear()
            self.axis_right_tab2.set_xlabel("m/z")
            self.axis_right_tab2.set_ylabel("Intensity")
            self.axis_right_tab2.set_title("Mass Spectrum")
            self.canvas_right_tab2.draw()

            self.step_display.delete(1.0, tk.END)
            self.step_display2.delete(1.0, tk.END)

            self.charge_table.delete(*self.charge_table.get_children())


            self.t_min.delete(0, tk.END)
            self.t_max.delete(0, tk.END)

            # Optional: reset dropdown
            # self.dropdown_var.set("Select Option")

            self.status_var.set("Status: Cleared all loaded data")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to clear data: {e}")

    def clear_step_display(self):
        self.step_display.delete(1.0, tk.END)
        self.deconvolution_count = 0

    def export_charge_table_as_excel(self):
        from openpyxl import Workbook
        from openpyxl.styles import PatternFill, Alignment
        import os

        if not self.charge_table.get_children():
            messagebox.showwarning("Warning", "No charge table data to export.")
            return

        if not self.output_dir or not self.base_name_file:
            messagebox.showwarning("Warning", "No output folder available. Please load data first.")
            return

        try:
            # Ensure output directory exists
            os.makedirs(self.output_dir, exist_ok=True)

            # Build descriptive filename
            try:
                t_min = self.t_min.get().strip()
                t_max = self.t_max.get().strip()
            except Exception:
                t_min, t_max = "NA", "NA"

            try:
                peak_dist = self.peaks_distance.get().strip()
            except Exception:
                peak_dist = "NA"

            file_name = (
                f"{self.base_name_file}_"
                f"time_range_{t_min}-to-{t_max}_"
                f"annotated_charge_state_input={peak_dist}.xlsx"
            )
            file_path = os.path.join(self.output_dir, file_name)

            # If file already exists, add a numeric suffix
            counter = 1
            base_file_path = file_path
            while os.path.exists(file_path):
                file_path = base_file_path.replace(".xlsx", f"_{counter}.xlsx")
                counter += 1

            wb = Workbook()
            ws = wb.active
            ws.title = "Charge Table"

            # Headers
            headers = ["Peak 1 (m/z)", "Peak 2 (m/z)", "Charge 1", "Charge 2", "Δm", "Mass (kDa)"]
            ws.append(headers)

            # Define fills for good/bad delta
            good_fill = PatternFill(start_color="98FB98", end_color="98FB98", fill_type="solid")  # pale green
            bad_fill = PatternFill(start_color="FFE4E1", end_color="FFE4E1", fill_type="solid")   # misty rose

            # Append rows with formatting
            for iid in self.charge_table.get_children():
                values = self.charge_table.item(iid)["values"]
                ws.append(values)
                row_idx = ws.max_row
                tags = self.charge_table.item(iid)["tags"]
                best_fill = PatternFill(start_color="98FB98", end_color="98FB98", fill_type="solid")
                mid_fill  = PatternFill(start_color="FFD580", end_color="FFD580", fill_type="solid")
                bad_fill  = PatternFill(start_color="FFB6C1", end_color="FFB6C1", fill_type="solid")

                if "best_delta" in tags:
                    fill = best_fill
                elif "mid_delta" in tags:
                    fill = mid_fill
                else:
                    fill = bad_fill

                for col in range(1, len(values) + 1):
                    ws.cell(row=row_idx, column=col).fill = fill

            # Center align cells
            for col in ws.columns:
                for cell in col:
                    cell.alignment = Alignment(horizontal="center", vertical="center")

            wb.save(file_path)
            self.status_var.set(f"Excel saved: {file_path}")
            messagebox.showinfo("Export Complete", f"Saved Excel file:\n{file_path}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to export Excel:\n{e}")

    def _save_current_figure(self, fig):
        """Return a temporary PNG bytes object from a matplotlib Figure."""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        return buf
    
    def export_as_pdf(self):
        if not self.file_path or not self.output_dir:
            messagebox.showwarning("Warning", "Please load data first.")
            return

        try:
            try:
                t_min = self.t_min.get().strip()
                t_max = self.t_max.get().strip()
            except Exception:
                t_min, t_max = "NA", "NA"
            
            
            avg_method = self.process_data.get().strip()

            pdf_name = (
                f"{self.base_name_file}_"
                f"time_range_{t_min}-to-{t_max}_{avg_method}MS-Buster.pdf"
            )
            pdf_path = os.path.join(self.output_dir, pdf_name)

            doc  = SimpleDocTemplate(pdf_path, pagesize=A4)
            styles = getSampleStyleSheet()
            flow  = []

            logo_path = "header.png"
            if os.path.exists(logo_path):
                flow.append(Image(logo_path, width=450, height=60))
                # flow.append(Spacer(1, 12))


            info = f"<b>Raw file:</b> {os.path.basename(self.file_path)}<br/>" \
                f"<b>Path:</b> {self.file_path}"
            flow.append(Paragraph(info, styles["Normal"]))
            # flow.append(Spacer(1, 12))


            if self.user_data is not None:
                flow.append(Paragraph("<b>Total Ion Chromatogram (TIC)</b>", styles["Heading2"]))
                buf1 = self._save_current_figure(self.fig_right_tab1)
                flow.append(Image(buf1, width=450, height=250))
                # flow.append(Spacer(1, 12))


            if self.ms_spectrum is not None:
                flow.append(Paragraph(f"<b>Raw averaged MS spectrum (selected time range {t_min}--{t_max})</b>", styles["Heading2"]))
                from matplotlib.figure import Figure
                import io
                fig_raw = Figure(figsize=(5, 3), dpi=100)
                ax_raw = fig_raw.add_subplot(111)
                ax_raw.plot(self.ms_spectrum[:, 0], self.ms_spectrum[:, 1], color="blue", linewidth=1)
                ax_raw.set_xlabel("m/z")
                ax_raw.set_ylabel("Intensity")
                ax_raw.set_title("Raw Averaged MS Spectrum")

                buf_raw = io.BytesIO()
                fig_raw.savefig(buf_raw, format="png", dpi=150, bbox_inches="tight")
                buf_raw.seek(0)

                flow.append(Image(buf_raw, width=450, height=250))


            if self.peaks_distance.get():
                flow.append(Paragraph(
                    f"<b>Peak-distance input:</b> {self.peaks_distance.get()}", styles["Normal"]
                ))
            try:
                    x_min = self.x_min.get().strip()
                    x_max = self.x_max.get().strip()
            except:
                    pass
            buf3 = self._save_current_figure(self.fig_right_tab2)
            flow.append(Paragraph(f"<b>Annotated MS plot with range ({x_min}--{x_max})</b>", styles["Heading2"]))
            flow.append(Image(buf3, width=450, height=250))
            # flow.append(Spacer(1, 12))

            if self.charge_table.get_children():
                flow.append(Paragraph("<b>Charge state table</b>", styles["Heading2"]))
                charge_rows = [["Peak1(m/z)","Peak2(m/z)","z1","z2","Δ","Mass(kDa)"]]
                for iid in self.charge_table.get_children():
                    vals = self.charge_table.item(iid)["values"]
                    delta = int(vals[4])
                    if delta == 1:
                        row_color = colors.HexColor("#98FB98")   # green
                    elif delta == 0:
                        row_color = colors.HexColor("#FFD580")   # yellow/orange
                    else:
                        row_color = colors.HexColor("#FFB6C1")   # pink/red

                    charge_rows.append(vals + [row_color])  

                # Keep header intact, strip color only from data rows
                table_data = [charge_rows[0]]  # header row
                for row in charge_rows[1:]:
                    table_data.append(row[:-1])  # drop color column

                tbl = Table(table_data, repeatRows=1)
                ts  = TableStyle([
                    ('BACKGROUND',(0,0),(-1,0),colors.lightblue),
                    ('GRID',(0,0),(-1,-1),0.25,colors.grey)
                ])
                for i, row in enumerate(charge_rows[1:], start=1):
                    ts.add('BACKGROUND', (0,i), (-1,i), row[-1])  # use saved row color
                tbl.setStyle(ts)
                flow.append(tbl)

                # flow.append(Spacer(1, 12))

            txt = self.step_display.get("1.0", tk.END).strip()
            if txt:
                # flow.append(PageBreak())
                flow.append(Paragraph("<b>Manual Deconvolution Results</b>", styles["Heading2"]))
                flow.append(Paragraph("<pre>%s</pre>" % txt.replace("\n","<br/>"), styles["Normal"]))

            doc.build(flow)
            self.status_var.set(f"PDF saved: {pdf_path}")
            messagebox.showinfo("Export Complete", f"Saved PDF:\n{pdf_path}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to export PDF:\n{e}")


    def show_peak_info(self):
        """Display peak information"""
        if not hasattr(self, 'peak_ranges') or not self.peak_ranges:
            return
        
        self.step_display.delete(1.0, tk.END)
        output = f"Detected {len(self.peak_ranges)} Chromatographic Peaks:\n"
        output += "=" * 50 + "\n\n"        
        output += "💡 Tips:\n"
        output += "• Click on any peak to auto-select its range\n"
        output += "• Use span selector for manual selection\n"
        output += "• Selected range will be used for spectrum extraction\n"
        total_time = max(self.time) - min(self.time)
        
        for i, peak_info in enumerate(self.peak_ranges):
            duration = peak_info['end_time'] - peak_info['start_time']
            width_percent = (duration / total_time) * 100
            
            output += f"Peak {i+1}:\n"
            output += f"  • Time Range: {peak_info['start_time']:.3f} - {peak_info['end_time']:.3f} min\n"
            output += f"  • Peak Maximum: {peak_info['peak_time']:.3f} min\n"
            output += f"  • Peak Intensity: {peak_info['peak_intensity']:.2e}\n"
            output += f"  • Peak Width: {duration:.3f} min ({width_percent:.1f}% of total)\n"
            output += f"  • Click to select this peak range\n\n"
        
        
        self.step_display.insert(tk.END, output)

    def detect_peaks_directly(self, time, intensity):
        """
        Direct peak detection without baseline subtraction
        """
        peak_height = np.max(intensity) * 0.03  # 3% of max intensity
        peak_distance = max(len(time) // 50, 5)  # Minimum distance between peaks
        
        peaks, properties = find_peaks(
            intensity,
            height=peak_height,
            distance=peak_distance,
            prominence=peak_height * 0.5
        )
        
        if len(peaks) == 0:
            return [], []
        
        try:
            widths, width_heights, left_ips, right_ips = peak_widths(
                intensity, peaks, rel_height=0.1  # 10% of peak height for boundaries
            )
            
            peak_ranges = []
            for i, peak_idx in enumerate(peaks):
                try:
                    left_idx = int(np.floor(left_ips[i]))
                    right_idx = int(np.ceil(right_ips[i]))
                    
                    left_idx = max(0, left_idx)
                    right_idx = min(len(time) - 1, right_idx)
                    
                    start_time = time[left_idx]
                    end_time = time[right_idx]
                    peak_time = time[peak_idx]
                    peak_intensity = intensity[peak_idx]
                    
                    peak_ranges.append({
                        'peak_idx': peak_idx,
                        'start_time': start_time,
                        'end_time': end_time,
                        'peak_time': peak_time,
                        'peak_intensity': peak_intensity,
                        'start_idx': left_idx,
                        'end_idx': right_idx
                    })
                except (IndexError, ValueError):
                    continue
                    
            return peaks, peak_ranges
            
        except Exception as e:
            print(f"Error calculating peak widths: {e}")
            peak_ranges = []
            for peak_idx in peaks:
                peak_ranges.append({
                    'peak_idx': peak_idx,
                    'start_time': time[max(0, peak_idx-10)],
                    'end_time': time[min(len(time)-1, peak_idx+10)],
                    'peak_time': time[peak_idx],
                    'peak_intensity': intensity[peak_idx],
                    'start_idx': max(0, peak_idx-10),
                    'end_idx': min(len(time)-1, peak_idx+10)
                })
            
            return peaks, peak_ranges
    
    def on_peak_click(self, event):
        """Handle clicks on peak ranges"""
        if event.inaxes != self.axis_right_tab1:
            return
        
        if not hasattr(self, 'peak_ranges') or not self.peak_ranges:
            return
            
        click_time = event.xdata
        if click_time is None:
            return
        
        closest_peak = None
        min_distance = float('inf')
        closest_peak_idx = -1
        
        for idx, peak_info in enumerate(self.peak_ranges):
            peak_time = peak_info['peak_time']
            distance = abs(peak_time - click_time)
            if distance < min_distance:
                min_distance = distance
                closest_peak = peak_info
                closest_peak_idx = idx
        
        if closest_peak and min_distance < (max(self.time) - min(self.time)) * 0.05:
            self.t_min.delete(0, tk.END)
            self.t_min.insert(0, f"{closest_peak['start_time']:.4f}")
            
            self.t_max.delete(0, tk.END)
            self.t_max.insert(0, f"{closest_peak['end_time']:.4f}")
            
            duration = closest_peak['end_time'] - closest_peak['start_time']
            self.status_var.set(f"Selected Peak {closest_peak_idx+1}: {closest_peak['start_time']:.3f}-{closest_peak['end_time']:.3f} min (Width: {duration:.3f} min)")
            
            self.highlight_selected_range(closest_peak)

    def highlight_selected_range(self, peak_info):
        """Highlight the selected peak range"""
        for patch in self.axis_right_tab1.patches:
            if hasattr(patch, 'peak_highlight'):
                patch.remove()
        
        start_time = peak_info['start_time']
        end_time = peak_info['end_time']
        
        highlight = self.axis_right_tab1.axvspan(
            start_time, end_time, 
            alpha=0.3, color='yellow', 
            label='Selected Range'
        )
        highlight.peak_highlight = True
        
        self.canvas_right_tab1.draw()

    def load_mzml_file(self):
        self.file_path = filedialog.askopenfilename(
            title="Select mzML File",
            filetypes=[("mzML Files", "*.mzML"), ("All files", "*.*")]
        )
        
        if not self.file_path:
            return

        self.status_var.set(f"Status: Loading mzML file: {os.path.basename(self.file_path)}...")
        self.root.update_idletasks()

        def process_mzml():
            try:
                # Heavy processing in background thread
                importer = mzML.MZMLImporter(self.file_path)
                tic_data = importer.get_tic()
                
                if tic_data is not None:
                    # Prepare data in background thread
                    user_data = tic_data
                    history = user_data.copy()
                    time_data = user_data[:, 0]
                    intensity_data = user_data[:, 1]
                    
                    base_name_file = os.path.splitext(os.path.basename(self.file_path))[0]
                    output_dir = os.path.join(os.path.dirname(self.file_path), f"{base_name_file}_MS-Buster")
                    
                    # Schedule GUI updates on main thread
                    def update_gui():
                        try:
                            self.reset_analysis_state()
                            self.importer = importer
                            self.mzfile_obj = None
                            self.user_data = user_data
                            self.history = history
                            self.time = time_data
                            self.intensity = intensity_data
                            self.base_name_file = base_name_file
                            self.output_dir = output_dir
                            
                            os.makedirs(self.output_dir, exist_ok=True)
                            output2_file = os.path.join(self.output_dir, f"{self.base_name_file}_raw_tic.txt")
                            np.savetxt(output2_file, self.user_data, delimiter="\t", fmt="%.4f\t%.2f", header="time\tIntensity")
                            
                            self.plot_data_on1()  # Now on main thread
                            self.status_var.set("Status: mzML file loaded successfully.")
                            
                        except Exception as e:
                            self.status_var.set(f"Status: Error updating GUI - {e}")
                    
                    # Schedule GUI update on main thread
                    self.root.after(0, update_gui)
                    
                else:
                    # Handle error on main thread
                    def show_error():
                        messagebox.showerror("Error", "Could not extract TIC data from the mzML file.")
                        self.status_var.set("Status: Error loading mzML.")
                    
                    self.root.after(0, show_error)
                    
            except Exception as e:
                # Handle error on main thread
                def show_error():
                    messagebox.showerror("Error", f"Failed to load mzML file: {e}")
                    self.status_var.set("Status: Error loading mzML.")
                
                self.root.after(0, show_error)

        # Start background thread
        threading.Thread(target=process_mzml, daemon=True).start()

    def load_agilent_folder(self):
        self.file_path = filedialog.askdirectory(title="Select Agilent .d Folder")
        if not self.file_path:
            return

        self.importer = None
        self.status_var.set("Status: Loading Agilent data...")
        self.root.update_idletasks()

        def process_agilent():
            try:
                mzfile_obj = mzAPI.mzFile(self.file_path)
                user_data = np.array(mzfile_obj.xic())
                history = user_data.copy()
                time_data = user_data[:, 0]
                intensity_data = user_data[:, 1]

                base_name_file = os.path.splitext(os.path.basename(self.file_path))[0]
                output_dir = os.path.join(os.path.dirname(self.file_path), f"{base_name_file}_MS-Buster")

                
                def update_gui():
                    try:
                        self.reset_analysis_state()
                        self.mzfile_obj = mzfile_obj
                        self.user_data = user_data
                        self.history = history
                        self.time = time_data
                        self.intensity = intensity_data
                        self.base_name_file = base_name_file
                        self.output_dir = output_dir

                        os.makedirs(self.output_dir, exist_ok=True)
                        output1_file = os.path.join(self.output_dir, f"{self.base_name_file}_raw_tic.txt")
                        np.savetxt(output1_file, self.user_data, delimiter="\t", fmt="%.4f\t%.2f", header="time\tIntensity")

                        self.plot_data_on1()

                        self.step_display2.delete(1.0, tk.END)
                        self.step_display2.insert(tk.END, "time\tIntensity\n")
                        for row in self.user_data:
                            self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")

                        self.status_var.set("Status: Agilent data loaded successfully")

                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to update GUI: {e}")
                        self.status_var.set("Status: Error updating GUI.")

                self.root.after(0, update_gui)

            except Exception as e:
                def show_error():
                    messagebox.showerror("Error", f"Failed to load Agilent folder: {e}")
                    self.status_var.set("Status: Error loading Agilent data.")
                self.root.after(0, show_error)

        threading.Thread(target=process_agilent, daemon=True).start()

    # def load_agilent_folder(self):
    #     self.file_path = filedialog.askdirectory(title="Select Agilent .d Folder")
    #     if not self.file_path:
    #             return
        
    #     self.importer = None
    #     self.status_var.set("Status: Loading Agilent data...")
    #     self.root.update_idletasks()

    #     def process_agilent():
    #         try:
    #             self.mzfile_obj = mzAPI.mzFile(self.file_path)
    #             self.user_data = np.array(self.mzfile_obj.xic())
    #             self.history = self.user_data.copy()
    #             self.time = self.user_data[:, 0]
    #             self.intensity = self.user_data[:, 1]
    #             self.plot_data_on1()
    #             self.step_display2.delete(1.0, tk.END)
    #             self.step_display2.insert(tk.END, "time\tIntensity\n")
    #             for row in self.user_data:
    #                 self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")

    #             self.status_var.set("Status: Data loaded successfully")

    #             self.base_name_file = os.path.splitext(os.path.basename(self.file_path))[0]
    #             self.output_dir = os.path.join(os.path.dirname(self.file_path), f"{self.base_name_file}_MS-Buster")
    #             os.makedirs(self.output_dir, exist_ok=True)
    #             output1_file = os.path.join(self.output_dir, f"{self.base_name_file}_raw_tic.txt")
    #             np.savetxt(output1_file, self.user_data,delimiter="\t", fmt="%.4f\t%.2f", header="time\tIntensity")


    #         except Exception as e:
    #             messagebox.showerror("Error", f"Failed to load data: {e}")
    #             self.status_var.set("Status: Error loading data")
            
    #     threading.Thread(target=process_agilent).start()


    def load_ms_from_agilent_folder(self):
        messagebox.showinfo(title="⚠️ WARNING ⚠️", 
                            message=("This function will average the ENTIRE DATA"
                                    "into a single mass spectrum. If it is a CHROMATOGRAM DATA this can cause noisy "
                                    "peaks to become prominent."))
        
        file_path = filedialog.askdirectory(title="Select Agilent .d Folder")
        if not file_path:
            return
            
        self.importer = None
        self.status_var.set(f"Status: Loading and averaging {os.path.basename(file_path)}...")
        self.root.update_idletasks()

        def loading_and_processing():
            try:
                ms_spectrum = self.get_ms_data_directly_from_agilent(file_path)

                def update_gui():
                    try:
                        self.reset_analysis_state()
                        self.ms_spectrum = ms_spectrum

                        if self.ms_spectrum.size == 0:
                            self.status_var.set("Status: No valid data found in the selected folder.")
                            return

                        self.m_z = self.ms_spectrum[:, 0]
                        self.ms_intensity = self.ms_spectrum[:, 1]
                        self.plot_data_on2()

                        self.step_display2.delete(1.0, tk.END)
                        self.step_display2.insert(tk.END, "m/z\tIntensity\n")
                        for row in self.ms_spectrum:
                            self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")

                        self.status_var.set("Status: Averaged data loaded successfully")

                        self.base_name_file = os.path.splitext(os.path.basename(file_path))[0]
                        self.output_dir = os.path.join(os.path.dirname(file_path), f"{self.base_name_file}_MSBUSTER")
                        os.makedirs(self.output_dir, exist_ok=True)
                        output_file = os.path.join(self.output_dir, f"{self.base_name_file}_averaged_spectrum.txt")
                        np.savetxt(output_file, self.ms_spectrum, delimiter="\t", fmt="%.4f\t%.2f", header="m/z\tIntensity")
                        print(f"Saved averaged spectrum to: {output_file}")

                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to update GUI: {e}")
                        self.status_var.set("Status: Error updating GUI.")

                self.root.after(0, update_gui)

            except Exception as e:
                def show_error():
                    messagebox.showerror("Error", f"Failed to load or process data: {e}")
                    self.status_var.set("Status: Error loading or processing data")
                self.root.after(0, show_error)

        threading.Thread(target=loading_and_processing, daemon=True).start()

    # def load_ms_from_agilent_folder(self):
    #     messagebox.showinfo(title="⚠️ WARNING ⚠️", 
    #                         message=("This function will average the ENTIRE DATA"
    #                                 "into a single mass spectrum. If it is a CHROMATOGRAM DATA this can cause noisy "
    #                                 "peaks to become prominent."))
        
    #     file_path = filedialog.askdirectory(title="Select Agilent .d Folder")
    #     if not file_path:
    #         return
            
    #     self.importer = None
    #     self.status_var.set(f"Status: Loading and averaging {os.path.basename(file_path)}...")
    #     self.root.update_idletasks()

    #     def loading_and_processing():
    #         try:
    #             self.reset_analysis_state()
    #             self.ms_spectrum = self.get_ms_data_directly_from_agilent(file_path)

    #             if self.ms_spectrum.size == 0:
    #                 self.status_var.set("status: No valid data found in the selected folder.")
    #                 return
                
    #             self.m_z = self.ms_spectrum[:,0]
    #             self.ms_intensity = self.ms_spectrum[:,1]
    #             self.plot_data_on2()

    #             self.step_display2.delete(1.0, tk.END)
    #             self.step_display2.insert(tk.END, "m/z\tIntensity\n")
    #             for row in self.ms_spectrum:
    #                 self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")


    #             self.status_var.set("Status: Averaged data loaded successfully")

    #             self.base_name_file = os.path.splitext(os.path.basename(file_path))[0]
    #             self.output_dir = os.path.join(os.path.dirname(self.file_path), f"{self.base_name_file}_MSBUSTER")
    #             os.makedirs(self.output_dir, exist_ok=True)
    #             output_file = os.path.join(self.output_dir, f"{self.base_name_file}_averaged_spectrum.txt")
    #             np.savetxt(output_file, self.ms_spectrum, delimiter="\t", fmt="%.4f\t%.2f", header="m/z\tIntensity")
    #             print(f"Saved averaged spectrum to: {output_file}")

    #         except Exception as e:
    #             messagebox.showerror("Error", f"Failed to load or process data: {e}")
    #             self.status_var.set("Status: Error loading or processing data")
        
    #     threading.Thread(target=loading_and_processing).start()

    def load_MS_from_txt_file(self):
        self.file_path = filedialog.askopenfilename(
            title="Select Text File",
            filetypes=[("Text Files", "*.txt")]
        )

        if not self.file_path:
            return

        self.status_var.set(f"Status: Loading data from {os.path.basename(self.file_path)}")
        self.root.update_idletasks()

        def loading():
            try:
                # Heavy work in background thread
                ms_spectrum = np.loadtxt(self.file_path, skiprows=1)

                def update_gui():
                    try:
                        # GUI work on main thread
                        self.reset_analysis_state()

                        self.ms_spectrum = ms_spectrum

                        if self.ms_spectrum.size == 0:
                            messagebox.showerror("Error", "No valid data found in the selected file.")
                            self.status_var.set("Status: No valid data found in the file.")
                            return

                        self.m_z = self.ms_spectrum[:, 0]
                        self.ms_intensity = self.ms_spectrum[:, 1]

                        self.plot_data_on2()

                        self.step_display2.delete(1.0, tk.END)
                        self.step_display2.insert(tk.END, "m/z\tIntensity\n")

                        for row in self.ms_spectrum:
                            self.step_display2.insert(
                                tk.END,
                                f"{row[0]:.4f}\t{row[1]:.4f}\n"
                            )

                        self.status_var.set("Status: Data loaded successfully")

                        self.base_name_file = os.path.splitext(
                            os.path.basename(self.file_path)
                        )[0]

                        self.output_dir = os.path.join(
                            os.path.dirname(self.file_path),
                            f"{self.base_name_file}_MS-Buster"
                        )

                        os.makedirs(self.output_dir, exist_ok=True)

                        output2_file = os.path.join(
                            self.output_dir,
                            f"{self.base_name_file}.txt"
                        )

                        np.savetxt(
                            output2_file,
                            self.ms_spectrum,
                            delimiter="\t",
                            fmt="%.4f\t%.2f",
                            header="m/z\tIntensity"
                        )

                    except Exception as e:
                        messagebox.showerror("Error", f"Failed to update GUI: {e}")
                        self.status_var.set("Status: Error updating GUI")

                self.root.after(0, update_gui)

            except Exception as e:

                def show_error():
                    messagebox.showerror("Error", f"Failed to load data: {e}")
                    self.status_var.set("Status: Error loading data")

                self.root.after(0, show_error)

        threading.Thread(target=loading, daemon=True).start()

    # def load_MS_from_txt_file(self):
    #     self.file_path = filedialog.askopenfilename(title = "Select Text File", filetypes=[("Text Files","*.txt")])
    #     if not self.file_path:
    #         return
        
    #     self.status_var.set(f"Status: Loading data from {os.path.basename(self.file_path)}")
    #     self.root.update_idletasks()

    #     def loading():
    #         try:
    #             self.reset_analysis_state()
    #             self.ms_spectrum = np.loadtxt(self.file_path, skiprows=1)
    #             self.m_z = self.ms_spectrum[:,0]
    #             self.ms_intensity = self.ms_spectrum[:,1]
    #             if self.ms_spectrum.size == 0:
    #                 messagebox.showerror("Error", "No valid data found in the selected file.")
    #                 self.status_var.set("Status: No valid data found in the file.")
    #                 return
    #             else:
    #                 self.plot_data_on2()
    #             self.step_display2.delete(1.0, tk.END)
    #             self.step_display2.insert(tk.END, "m/z\tIntensity\n")
    #             for row in self.ms_spectrum:
    #                 self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")
    #             self.status_var.set("Status: Data loaded successfully")
    #             self.base_name_file = os.path.splitext(os.path.basename(self.file_path))[0]
    #             self.output_dir = os.path.join(os.path.dirname(self.file_path), f"{self.base_name_file}_MS-Buster")
    #             os.makedirs(self.output_dir, exist_ok=True)
    #             output2_file = os.path.join(self.output_dir, f"{self.base_name_file}.txt")
    #             np.savetxt(output2_file, self.ms_spectrum,delimiter="\t", fmt="%.4f\t%.2f", header="m/z\tIntensity")

    #         except Exception as e:
    #             messagebox.showerror("Error", f"Failed to load data: {e}")
    #             self.status_var.set("Status: Error loading data")

    #     threading.Thread(target=loading).start()

    def apply_deisotoping_to_spectrum(self, spectrum_data, isotope_mass_diff=1.003355, intensity_ratio_threshold=0.8):
            """
            Apply basic deisotoping to remove likely isotope peaks.
            
            Args:
                spectrum_data (numpy.ndarray): 2D array with columns [m/z, intensity]
                isotope_mass_diff (float): Expected mass difference for isotopes (default: 1.003355 Da)
                intensity_ratio_threshold (float): Maximum intensity ratio for isotope peaks
            
            Returns:
                numpy.ndarray: Deisotoped spectrum
            """
            if spectrum_data is None or len(spectrum_data) == 0:
                return spectrum_data

            # First apply centroiding to get clean peaks
            centroided_data = self.apply_centroiding_to_spectrum(spectrum_data)
            
            if len(centroided_data) <= 1:
                return centroided_data

            # Sort by m/z for isotope detection
            sorted_indices = np.argsort(centroided_data[:, 0])
            sorted_data = centroided_data[sorted_indices]
            
            # Track which peaks to keep (start with all True)
            keep_peaks = np.ones(len(sorted_data), dtype=bool)
            
            # Check each peak against previous peaks for isotope patterns
            for i in range(1, len(sorted_data)):
                current_mz = sorted_data[i, 0]
                current_intensity = sorted_data[i, 1]
                
                # Check against previous peaks within isotope range
                for j in range(i):
                    if not keep_peaks[j]:  # Skip already removed peaks
                        continue
                        
                    prev_mz = sorted_data[j, 0]
                    prev_intensity = sorted_data[j, 1]
                    
                    mass_diff = current_mz - prev_mz
                    
                    # Check for isotope patterns (C13, etc.)
                    if 0.9 < mass_diff < 1.1:  # +1 isotope
                        if current_intensity < prev_intensity * intensity_ratio_threshold:
                            keep_peaks[i] = False
                            break
                    elif 1.9 < mass_diff < 2.1:  # +2 isotope
                        if current_intensity < prev_intensity * (intensity_ratio_threshold ** 2):
                            keep_peaks[i] = False
                            break

            # Return only the peaks that are not isotopes
            deisotoped_spectrum = sorted_data[keep_peaks]
            print(f"Deisotoping: {len(centroided_data)} peaks → {len(deisotoped_spectrum)} monoisotopic peaks")
            
            return deisotoped_spectrum

    def apply_centroiding_to_spectrum(self, spectrum_data, intensity_threshold_ratio=0.01, peak_distance=3, prominence_ratio=0.005):
        """
        Apply centroiding to an MS spectrum from mzML file.
        
        Args:
            spectrum_data (numpy.ndarray): 2D array with columns [m/z, intensity]
            intensity_threshold_ratio (float): Minimum intensity as ratio of max intensity (default: 0.01 = 1%)
            peak_distance (int): Minimum distance between peaks in data points (default: 3)
            prominence_ratio (float): Minimum peak prominence as ratio of max intensity (default: 0.005 = 0.5%)
        
        Returns:
            numpy.ndarray: Centroided spectrum containing only peak maxima
        """
        if spectrum_data is None or len(spectrum_data) == 0:
            return spectrum_data

        mz = spectrum_data[:, 0]
        intensity = spectrum_data[:, 1]
        
        if len(intensity) == 0:
            return spectrum_data

        max_intensity = np.max(intensity)
        if max_intensity == 0:
            return spectrum_data

        # Calculate thresholds
        height_threshold = max_intensity * intensity_threshold_ratio
        prominence_threshold = max_intensity * prominence_ratio

        # Find peaks using scipy's find_peaks with multiple criteria
        peaks, properties = find_peaks(
            intensity,
            height=height_threshold,
            distance=peak_distance,
            prominence=prominence_threshold
        )

        if len(peaks) == 0:
            # If no peaks found, return empty spectrum with same structure
            return np.array([]).reshape(0, 2)

        # Return only the peak positions
        centroided_spectrum = spectrum_data[peaks]
        
        print(f"Centroiding: {len(spectrum_data)} points → {len(centroided_spectrum)} peaks")
        return centroided_spectrum

    def _set_averaging_controls(self, running):
        self.averaging_in_progress = running
        state = tk.DISABLED if running else tk.NORMAL
        try:
            self.extract_spectrum.configure(state=state)
            self.process_data.configure(state=state)
        except Exception:
            pass

    def _preview_spectrum_text(self, header, spectrum, max_rows=5000):
        rows = [header, "m/z\tIntensity\n"]
        preview = spectrum[:max_rows]
        rows.extend(f"{row[0]:.4f}\t{row[1]:.4f}\n" for row in preview)
        if len(spectrum) > max_rows:
            rows.append(f"\nShowing first {max_rows:,} of {len(spectrum):,} points. Full spectrum was saved to the output file.\n")
        return ''.join(rows)

    def _downsample_xy_for_display(self, x, y, max_points=20000):
        if x is None or y is None or len(x) <= max_points:
            return x, y
        indices = np.linspace(0, len(x) - 1, max_points, dtype=int)
        return x[indices], y[indices]

    def data_extractor(self, tmin=None, tmax=None):
        if not self.file_path:
            messagebox.showwarning("Warning", "No file loaded. Please load data first.")
            return

        if self.averaging_in_progress:
            messagebox.showinfo("Info", "Averaging is already running. Please wait for it to finish.")
            return

        try:
            t_min = float(tmin) if tmin is not None else float(self.t_min.get())
            t_max = float(tmax) if tmax is not None else float(self.t_max.get())
            if t_min >= t_max:
                raise ValueError("Minimum time must be less than maximum time")
        except ValueError as ve:
            messagebox.showerror("Error", f"Invalid time input: {ve}")
            return

        processing_mode = self.process_data.get()
        importer = self.importer
        mzfile_obj = self.mzfile_obj
        file_path = self.file_path
        output_dir = self.output_dir
        base_name_file = self.base_name_file

        self._set_averaging_controls(True)
        self.status_var.set("Status: Extracting and Averaging data")
        self.root.update_idletasks()

        def work():
            try:
                spectra_to_merge = []
                scan_info = None  # Initialize here
                spectra_count = 0
                ms_spectrum = None

                if importer is not None:  # mzML file processing
                    if processing_mode == "Default":
                        print("Using mzML built-in averaging for Default mode")
                        ms_spectrum = importer.avg_safe(time_range=[t_min, t_max], max_points=50000)
                    else:
                        scan_range = importer.get_scans_from_times([t_min, t_max])
                        if not scan_range:
                            self.root.after(0, lambda: self._finish_averaging_with_message("No scans found in the selected time range."))
                            return

                        importer.reset_reader()
                        
                        for spectrum in importer.msrun:
                            if scan_range[0] <= spectrum.ID <= scan_range[1]:
                                try:
                                    scan_data = mzML.get_data_from_spectrum(spectrum)
                                    if scan_data is not None and len(scan_data) > 0:
                                        if processing_mode == "Centroided":
                                            scan_data = self.apply_centroiding_to_spectrum(scan_data)
                                        elif processing_mode == "Deisotoping":
                                            scan_data = self.apply_deisotoping_to_spectrum(scan_data)
                                        
                                        if scan_data is not None and len(scan_data) > 0:
                                            spectra_to_merge.append(scan_data)
                                            spectra_count += 1
                                except Exception as e:
                                    print(f"Error processing spectrum {spectrum.ID}: {e}")
                                    continue
                        importer.reset_reader()

                        print(f"Processed {spectra_count} spectra in {processing_mode} mode")
                        
                        if spectra_to_merge:
                            ms_spectrum = mzML.merge_spectra(spectra_to_merge, max_points=50000)
                        else:
                            self.root.after(0, lambda: self._finish_averaging_with_message("No valid spectra found in the selected time range."))
                            return

                elif mzfile_obj is not None:  # Agilent .d folder processing
                    local_mzfile_obj = mzFile(file_path)
                    scan_info = local_mzfile_obj.scan_info(t_min, t_max, None, None)
                    
                    if not scan_info:
                        self.root.after(0, lambda: self._finish_averaging_with_message("No scans found in the selected time range."))
                        return

                    for time, _, scan_number, *_ in scan_info:
                        if processing_mode == "Centroided":
                            scan_data = local_mzfile_obj.cscan(scan_number)
                        elif processing_mode == "Deisotoping":
                            scan_data = local_mzfile_obj.deisotope_scan(scan_number)
                        else:
                            scan_data = local_mzfile_obj.scan(scan_number)
                        
                        if scan_data is not None and len(scan_data) > 0:
                            spectra_to_merge.append(np.array(scan_data))

                    if spectra_to_merge:
                        ms_spectrum = mzML.merge_spectra(spectra_to_merge, max_points=50000)
                

                else:
                    self.root.after(0, lambda: self._finish_averaging_with_message("No data loaded.", warning=True))
                    return
                
                if ms_spectrum is not None and len(ms_spectrum) > 0:
                    output_file = None
                    if output_dir and base_name_file:
                        os.makedirs(output_dir, exist_ok=True)
                        
                        if processing_mode == "Default":
                            output_file = os.path.join(output_dir, f"{base_name_file}_{t_min}-{t_max}_mins_normal_averaged_MS_Spectrum.txt")
                        elif processing_mode == "Centroided":
                            output_file = os.path.join(output_dir, f"{base_name_file}_{t_min}-{t_max}_mins_centroided_MS_Spectrum.txt")
                        elif processing_mode == "Deisotoping":
                            output_file = os.path.join(output_dir, f"{base_name_file}__{t_min}-{t_max}_mins_deisotoped_MS_Spectrum.txt")
                        
                        if output_file:
                            np.savetxt(output_file, ms_spectrum, delimiter="\t", fmt="%.4f\t%.2f",
                                    header="m/z\tIntensity", comments='')

                    self.root.after(0, lambda: self._finish_averaging_success(
                        ms_spectrum, processing_mode, t_min, t_max, scan_info, spectra_count, output_file
                    ))
                else:
                    self.root.after(0, lambda: self._finish_averaging_no_data())

            except Exception as e:
                print(traceback.format_exc())
                self.root.after(0, lambda error=e: self._finish_averaging_error(error))

        threading.Thread(target=work, daemon=True).start()

    def _finish_averaging_success(self, ms_spectrum, processing_mode, t_min, t_max, scan_info, spectra_count, output_file):
        try:
            self.ms_spectrum = ms_spectrum
            self.m_z = self.ms_spectrum[:, 0]
            self.ms_intensity = self.ms_spectrum[:, 1].copy()

            max_intensity = np.max(self.ms_intensity)
            if max_intensity > 1000:
                self.ms_intensity = (self.ms_intensity / max_intensity) * 100
                print(f"Normalized intensity from max {max_intensity:.2e} to 100")

            self.plot_data_on2()

            if self.importer:
                result = f"Processed mzML data between {t_min}-{t_max} mins.\n"
                if spectra_count:
                    result += f"Processed {spectra_count} spectra.\n"
                result += f"Averaged {processing_mode} Mass Spectrum\n"
            elif scan_info:
                result = f"Found {len(scan_info)} scans between {t_min}-{t_max} mins.\n"
                result += f"Averaged {processing_mode} Mass Spectrum\n"
            else:
                result = f"Processed data between {t_min}-{t_max} mins.\n"
                result += f"Averaged {processing_mode} Mass Spectrum\n"

            if output_file:
                result += f"Saved full spectrum to: {output_file}\n"

            self.step_display2.delete(1.0, tk.END)
            self.step_display2.insert(tk.END, self._preview_spectrum_text(result, self.ms_spectrum))
            self.status_var.set(f"Status: Averaged {processing_mode} spectrum displayed.")
        finally:
            self._set_averaging_controls(False)

    def _finish_averaging_no_data(self):
        self.status_var.set("Status: No spectral data found in range.")
        self._set_averaging_controls(False)

    def _finish_averaging_with_message(self, message, warning=False):
        if warning:
            messagebox.showwarning("Warning", message)
        else:
            messagebox.showinfo("Info", message)
        self.status_var.set("Status: Ready")
        self._set_averaging_controls(False)

    def _finish_averaging_error(self, error):
        messagebox.showerror("Error", f"Processing failed: {error}")
        self.status_var.set(f"Status: Error - {error}")
        self._set_averaging_controls(False)

    def save_ms_spectrum_between_time_range(self):
        import os 
        if not self.file_path:
            messagebox.showwarning("Warning", "No file loaded. Please load Agilent data first.")
            
        self.status_var.set("Status: Extracting and storing all the MS spectras between time range")
        self.root.update_idletasks()
        
        def saving():
            try:
                t_min = float(self.t_min.get())
                t_max = float(self.t_max.get())


                output_dir_extracted_spectrums = os.path.join(self.output_dir, f"{self.base_name_file}_extracted_spectrum_from_{t_min}_to_{t_max}")
                os.makedirs(output_dir_extracted_spectrums, exist_ok=True)
                
                
                self.mzfile_obj = mzFile(self.file_path)
                scan_info = self.mzfile_obj.scan_info(t_min, t_max, None, None)

                for time, _, scan_number, *_ in scan_info:
                    spectrum = self.mzfile_obj.scan(scan_number)
                    spectrum = np.array(spectrum)

                    if spectrum.size == 0:
                        continue  
                    
                    output_file = os.path.join(output_dir_extracted_spectrums, f"scan_{scan_number}_time_{time:.2f}.txt")
                    np.savetxt(output_file, spectrum, delimiter="\t", fmt="%.4f\t%.2f", header="m/z\tIntensity", comments='')

    
                self.status_var.set(f"Status: {len(scan_info)} scans extracted & saved successfully as txt file between {t_min} and {t_max} minutes.")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to extract and save data: {e}")
                self.status_var.set("Status: Error loading data")
            
        threading.Thread(target=saving).start()

    def annotate_charge_state(self):
        self.peaks, self.properties = find_peaks(
            self.ms_intensity,
            distance=int(self.peaks_distance.get()),
        )
        self.peaks_mz = self.m_z[self.peaks]
        self.peaks_intensity = self.ms_intensity[self.peaks]

        self.sorted_peaks_indices = np.argsort(self.peaks_mz)
        self.sorted_peaks_mz = self.peaks_mz[self.sorted_peaks_indices]
        self.sorted_peaks_intensity = self.peaks_intensity[self.sorted_peaks_indices]
        self.charge_analysis_data = []

        for i in range(len(self.sorted_peaks_mz) - 1):
            m_n = self.sorted_peaks_mz[i]
            m_n1 = self.sorted_peaks_mz[i+1]
            intensity_n = self.sorted_peaks_intensity[i]
            m_p = 1.0073

            delta = m_n1 - m_n
            if delta > 0:
                charge = (m_n1 - m_p) / delta
                peak_data = {
                    'mz': m_n,
                    'intensity': intensity_n,
                    'charge': charge,
                    'rounded_charge': np.round(charge)
                }
                self.charge_analysis_data.append(peak_data)


        xlim = self.axis_right_tab2.get_xlim()
        ylim = self.axis_right_tab2.get_ylim()
        was_zoomed = xlim != (0.0, 1.0) and ylim != (0.0, 1.0)

        self.plot_data_on2()

        if was_zoomed:
            self.axis_right_tab2.set_xlim(xlim)
            self.axis_right_tab2.set_ylim(ylim)

        for peak_info in self.charge_analysis_data:
            mz = peak_info['mz']
            intensity = peak_info['intensity']
            charge_state = peak_info['charge']

            self.axis_right_tab2.annotate(
                f"{charge_state:.1f}+",
                xy=(mz, intensity),
                xytext=(0, 10),
                textcoords='offset points',
                ha='center',
                fontsize=6.5,
                color='red',
                arrowprops=dict(arrowstyle='->', color='blue'),
            )

        self.canvas_right_tab2.draw()
        self.update_charge_table()
  
    # def annotate_charge_state(self):
    #     self.peaks, self.properties = find_peaks(
    #         self.ms_intensity,
    #         distance=int(self.peaks_distance.get()),
    #     )
    #     self.peaks_mz = self.m_z[self.peaks]
    #     self.peaks_intensity = self.ms_intensity[self.peaks]

    #     self.sorted_peaks_indices = np.argsort(self.peaks_mz)

    #     self.sorted_peaks_mz = self.peaks_mz[self.sorted_peaks_indices]
    #     self.sorted_peaks_intensity = self.peaks_intensity[self.sorted_peaks_indices]
    #     self.charge_analysis_data = []


    #     for i in range(len(self.sorted_peaks_mz) - 1):
    #         m_n = self.sorted_peaks_mz[i]
    #         m_n1 = self.sorted_peaks_mz[i+1]
            
    #         intensity_n = self.sorted_peaks_intensity[i]

            
    #         m_p = 1.0073  

    #         delta = m_n1 - m_n
    #         if delta > 0:
    #             charge = (m_n1 - m_p) / delta
    #             peak_data = {
    #                 'mz': m_n,
    #                 'intensity': intensity_n,
    #                 'charge': charge,
    #                 'rounded_charge': np.round(charge)
    #             }
    #             self.charge_analysis_data.append(peak_data)

    #     self.plot_data_on2() 

    #     for peak_info in self.charge_analysis_data:
    #         mz = peak_info['mz']
    #         intensity = peak_info['intensity']
    #         charge_state = peak_info['charge']

    #         self.axis_right_tab2.annotate(
    #             f"{charge_state:.1f}+",  
    #             xy=(mz, intensity),
    #             xytext=(0, 10),
    #             textcoords='offset points',
    #             ha='center',
    #             fontsize=6.5,
    #             color='red',
    #             arrowprops=dict(arrowstyle='->', color='blue'),
    #         )


    #     self.canvas_right_tab2.draw()

    #     self.update_charge_table()   

    def Deconvolute(self):
        start_time = time.time()
        from scipy.signal import find_peaks

        # Initialize counter if not already
        if not hasattr(self, "deconvolution_count"):
            self.deconvolution_count = 0

        self.deconvolution_count += 1

        peaks, properties = find_peaks(
            self.ms_intensity,
            height=0.05,
            distance=3,
            prominence=0.2
        )
        
        mz_peaks_all = self.m_z[peaks]
        intensity_peaks_all = self.ms_intensity[peaks]

        sorted_indices = np.argsort(intensity_peaks_all)[::-1]
        top_n = 2
        top_indices = sorted_indices[:top_n]

        mz_peaks = mz_peaks_all[top_indices]
        intensity_peaks = intensity_peaks_all[top_indices]

        sort_order = np.argsort(mz_peaks)
        mz_peaks = mz_peaks[sort_order]
        intensity_peaks = intensity_peaks[sort_order]

        try:
            t_min = self.t_min.get().strip()
            t_max = self.t_max.get().strip()
        except Exception:
            t_min, t_max = "NA", "NA"

        try:
            process = self.process_data.get().strip()
        except Exception:
            process = "NA"

        
        output = f"Manual Peak Analysis {self.deconvolution_count}\n"
        output += f"Between {t_min} - {t_max} mins\n"
        output += f"Averaged {process} Mass Spectrum\n"
        output += "Selected Peaks (m/z and intensity):\n"
        for i in range(len(mz_peaks)):
            output += f"Peak {i+1}: m/z = {mz_peaks[i]:.4f}, Intensity = {intensity_peaks[i]:.3f}\n"

        if len(mz_peaks) == 2:
            m_n1 = mz_peaks[0]   
            m_n = mz_peaks[1]    
            m_p = 1.0073         

            n = (m_n1 - m_p) / (m_n - m_n1)
            M = n * (m_n - m_p) / 1000

            output += f"\nCalculated charge state (n): {n:.2f}\n"
            output += f"Calculated mass (M): {M:.3f} kDa ({M * 1000:.3f} Da)\n"
        else:
            output += "\nNot enough qualified peaks found for charge/mass calculation.\n"
        
        output += "\n" + "-" * 40 + "\n\n"


        if self.deconvolution_count == 1:
            self.step_display.delete(1.0, tk.END)
            self.step_display.insert(tk.END, output)
        else:
            self.step_display.insert("1.0", output)

        # Finish
        elapsed_time = time.time() - start_time
        self.status_var.set(
            f"Status: Deconvolution completed successfully in {elapsed_time:.4f} seconds"
        )

    def _calculate_mass_from_pair(self, peak1, peak2):
        """
        Calculates the charge and mass from a pair of peak dictionaries.
        Assumes peak1 and peak2 are sorted by m/z.
        """
        try:

            mz1 = min(peak1['mz'], peak2['mz'])
            mz2 = max(peak1['mz'], peak2['mz'])

            m_p = 1.0073  # Mass of a proton
            delta_mz = mz2 - mz1

            if delta_mz == 0:
                return None, None # Cannot divide by zero


            n = (mz1 - m_p) / delta_mz
            M = n * (mz2 - m_p) / 1000  # Mass in kDa

            return M

        except (KeyError, ZeroDivisionError) as e:
            print(f"Could not calculate mass for pair: {e}")
            return None
        
    def update_charge_table(self):
            """
            New pipeline: Filters peaks, groups them into pairs, runs deconvolution
            on each pair, and displays the results in the table and text box.
            """
            self.charge_table.delete(*self.charge_table.get_children())
            # self.step_display.delete(1.0, tk.END) # Clear the text box too

            if not hasattr(self, 'charge_analysis_data') or not self.charge_analysis_data:
                return

            data = self.charge_analysis_data
            valid_indices = set()
            if len(data) > 1:
                for i in range(len(data) - 1):
                    current_charge = data[i]['rounded_charge']
                    next_charge = data[i+1]['rounded_charge']
                    if abs(current_charge - next_charge) <= 1: # Simplified check
                        valid_indices.add(i)
                        valid_indices.add(i + 1)
            peaks_to_process = [data[i] for i in sorted(list(valid_indices))]

            output_string = "Deconvolution Results per Pair:\n"
            output_string += "="*40 + "\n"
            
            for i in range(0, len(peaks_to_process) - 1, 2):
                peak1 = peaks_to_process[i]
                peak2 = peaks_to_process[i+1]
                delta_charge = abs(int(peak2['rounded_charge']) - int(peak1['rounded_charge']))
                charge1_disp = f"{int(round(peak1['rounded_charge']))} ({peak1['charge']:.1f})"
                charge2_disp = f"{int(round(peak2['rounded_charge']))} ({peak2['charge']:.1f})"


                mass_kda = self._calculate_mass_from_pair(peak1, peak2)

                if mass_kda is not None:

                    values = (
                        f"{peak1['mz']:.4f}",
                        f"{peak2['mz']:.4f}",
                        charge1_disp,
                        charge2_disp,
                        f"{delta_charge}",
                        f"{mass_kda:.3f}"
                    )
                    if delta_charge == 1:
                        tag = "best_delta"      # green
                    elif delta_charge == 0:
                        tag = "mid_delta"       # yellow
                    else:
                        tag = "bad_delta"       # red

                    self.charge_table.insert('', tk.END, values=values, tags=(tag,))

                    output_string += f"Pair {i//2 + 1}:\n"
                    output_string += f"  Peak 1: {peak1['mz']:.4f} m/z (Charge: {int(peak1['rounded_charge'])}+)\n"
                    output_string += f"  Peak 2: {peak2['mz']:.4f} m/z (Charge: {int(peak2['rounded_charge'])}+)\n"
                    output_string += f"  Δz (charge diff): {delta_charge}\n"
                    output_string += f"  > Calculated Mass: {mass_kda:.3f} kDa\n\n"

            if len(output_string) < 60: 
                output_string += "No valid pairs found for deconvolution."
    
    def toggle_peak_display(self):
        """
        Called when the 'Show Detected Peaks' checkbox is toggled.
        Triggers a replot of the TIC.
        """
        self.plot_data_on1()

    def plot_data_on1(self):
        if self.time is None or self.intensity is None:
            messagebox.showwarning("Warning", "No data to plot. Please load data first.")
            return

        self.axis_right_tab1.clear()
        
        self.axis_right_tab1.plot(self.time, self.intensity, 'b-', label='TIC', linewidth=1.5)

        if self.show_peaks_var.get():
            # If the toggle is ON, detect and plot peaks
            peaks, peak_ranges = self.detect_peaks_directly(self.time, self.intensity)
            self.peak_ranges = peak_ranges
            
            if len(peaks) > 0:
                for i, peak_info in enumerate(peak_ranges):
                    start_time = peak_info['start_time']
                    end_time = peak_info['end_time']
                    
                    # Highlight peak region with green background
                    self.axis_right_tab1.axvspan(
                        start_time, end_time, 
                        alpha=0.2, color='green',
                        label='Peak Region' if i == 0 else ""
                    )
                    
                    # Mark peak boundaries with red dashed lines
                    self.axis_right_tab1.axvline(start_time, color='red', linestyle='--', alpha=0.7, linewidth=1)
                    self.axis_right_tab1.axvline(end_time, color='red', linestyle='--', alpha=0.7, linewidth=1)
                    
                    # Mark peak maximum with red circle
                    self.axis_right_tab1.plot(
                        peak_info['peak_time'], peak_info['peak_intensity'], 
                        'ro', markersize=6
                    )
                    
                    # Add peak number annotation
                    self.axis_right_tab1.annotate(
                        f'P{i+1}', 
                        xy=(peak_info['peak_time'], peak_info['peak_intensity']),
                        xytext=(0, 15), textcoords='offset points',
                        ha='center', fontsize=10, color='darkblue', weight='bold',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.8)
                    )

                # Get the current y-axis limits after plotting
                ymin, ymax = self.axis_right_tab1.get_ylim()
                
                # Increase the upper limit (ymax) by 20% to create padding for the yellow boxes
                self.axis_right_tab1.set_ylim(ymin, ymax * 1.20) 
                
                self.status_var.set(f"Status: Detected {len(peaks)} peaks. Click a peak to select its range.")
                self.show_peak_info()
            else:
                self.status_var.set("Status: No peaks detected in TIC.")
                self.step_display.delete(1.0, tk.END) # Clear peak info if none are found
                self.peak_ranges = []

        else:
            # If the toggle is OFF, just show the raw TIC
            self.status_var.set("Status: Displaying Raw TIC.")
            self.step_display.delete(1.0, tk.END) # Clear the peak info text box
            self.peak_ranges = [] # Clear peak ranges so clicks do nothing
            
        self.axis_right_tab1.set_xlabel("Time (min)")
        self.axis_right_tab1.set_ylabel('Intensity')
        self.axis_right_tab1.legend().remove()
        
        self.fig_right_tab1.tight_layout()
        
        # Connect events (always needed)
        self.canvas_right_tab1.mpl_connect('button_press_event', self.on_peak_click)
        self.canvas_right_tab1.draw()
        

        self.span_selector = SpanSelector(
            self.axis_right_tab1,
            self.on_select_span,
            direction='horizontal',
            useblit=True,
            props=dict(alpha=0.3, facecolor="grey"),
            interactive=True
        )

    def plot_data_on2(self):
        if self.m_z is None or self.ms_intensity is None:
            messagebox.showwarning("Warning", "No data to plot. Please load data first.")
            return
        
        self.axis_right_tab2.clear()
        plot_mz, plot_intensity = self._downsample_xy_for_display(self.m_z, self.ms_intensity)
        self.axis_right_tab2.plot(plot_mz, plot_intensity)
        self.axis_right_tab2.set_xlabel("m/z")
        self.axis_right_tab2.set_ylabel('Intensity')
        legend = self.axis_right_tab2.get_legend()
        if legend:
            legend.remove()
        self.canvas_right_tab2.draw()
        self.span_selector2 = SpanSelector(
                self.axis_right_tab2,
                self.on_select_span2,
                direction='horizontal',
                useblit=True,
                props=dict(alpha=0.3, facecolor="grey"),
                interactive=True
            )
                    
    def set_x_range(self):
        if self.ms_spectrum is None:
            messagebox.showwarning("Warning", "No data to set range. Please load data first.")
            return
        
        try:
            x_min = float(self.x_min.get())
            x_max = float(self.x_max.get())
            if x_min >= x_max:
                raise ValueError("Minimum value must be less than maximum value.")

            mask = (self.ms_spectrum[:, 0] >= x_min) & (self.ms_spectrum[:, 0] <= x_max)
            self.filtered_ms_data = self.ms_spectrum[mask]

            self.m_z = self.filtered_ms_data[:, 0]
            self.ms_intensity = self.filtered_ms_data[:, 1]
            self.plot_data_on2()

            self.step_display2.delete(1.0, tk.END)
            self.step_display2.insert(tk.END, "Filtered Mass Spectrum\nm/z\t Intensity\n")
            for row in self.filtered_ms_data:
                self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")
            output7_file = os.path.join(self.output_dir, f"{self.base_name_file}_fil_ms_spec_{x_min}-{x_max}.txt")
            np.savetxt(output7_file, self.ms_spectrum,delimiter="\t", fmt="%.4f\t%.2f", header="m/z\tIntensity")
            
            self.status_var.set("Status: m/z range set successfully")
        except ValueError as ve:
            messagebox.showerror("Error", f"Invalid input: {ve}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to set m/z range: {e}")

    def get_scans_from_times(self, time_range):
        if not self.mzfile_obj:
            messagebox.showerror("Error", "File object not found. Please load a file first.")
            return None
            
        start_time, end_time = time_range
        scan_info = self.mzfile_obj.scan_info(start_time, end_time)
        if not scan_info:
            return None
        start_scan = scan_info[0][2]
        end_scan = scan_info[-1][2]
        return (start_scan, end_scan + 1) 

    def get_ms_data_directly_from_agilent(self, file_path):
        self.file_path = file_path
        try:
            self.mzfile_obj = mzAPI.mzFile(self.file_path)
        except Exception as e:
            print(f"Failed to open file {self.file_path}: {e}")
            messagebox.showerror("File Error", f"Could not open or read the file:\n{self.file_path}\n\nError: {e}")
            return np.array([])

        full_scan_range = self.mzfile_obj.scan_range()
        print(f"Processing all scans in: {os.path.basename(self.file_path)}")
        
        try:
            # Get data with improved error handling
            averaged_spectrum = self.get_data_safe(scan_range=full_scan_range)
            
            if averaged_spectrum is None or len(averaged_spectrum) == 0:
                print("No valid data found in the file")
                return np.array([]).reshape(0, 2)
                
            return averaged_spectrum
            
        except Exception as e:
            print(f"Error processing Agilent data: {e}")
            messagebox.showerror("Processing Error", f"Failed to process data: {e}")
            return np.array([]).reshape(0, 2)

    def get_data_safe(self, scan_range=None, time_range=None):
        """Improved version of get_data with better error handling"""
        if time_range:
            scan_range = self.get_scans_from_times(time_range)
            
        data_to_merge = self.grab_data_safe(scan_range)
        
        if len(data_to_merge) == 0:
            print("No valid data found in the specified range")
            return np.array([]).reshape(0, 2)

        if len(data_to_merge) > 1:
            print(f"Merging {len(data_to_merge)} spectra...")
            try:
                merged_data = mzML.merge_spectra(data_to_merge)
                return merged_data
            except Exception as e:
                print(f"Error in merge_spectra: {e}")
                return np.array([]).reshape(0, 2)
        else:
            return data_to_merge[0]

    def grab_data_safe(self, scan_range):
        """Improved version of grab_data with data validation"""
        data_list = []
        
        if scan_range is None:
            return data_list

        print(f"Grabbing data for scan range: {scan_range}")
        
        for s in range(scan_range[0], scan_range[1]):
            try:
                # Get the scan data using the stored file object
                scan_data = self.mzfile_obj.scan(int(s))
                
                if scan_data is None or len(scan_data) == 0:
                    continue
                    
                # Convert to numpy array with proper error handling
                try:
                    impdat = np.array(scan_data, dtype=np.float64)
                except (ValueError, TypeError) as e:
                    print(f"Could not convert scan {s} to numeric data: {e}")
                    continue
                    
                if impdat.size == 0:
                    continue
                    
                # Ensure proper shape
                if impdat.ndim != 2 or impdat.shape[1] != 2:
                    print(f"Invalid shape for scan {s}: {impdat.shape}")
                    continue
                    
                # Filter out m/z values below 10 and invalid values
                valid_mask = ((impdat[:, 0] > 10) & 
                            np.isfinite(impdat[:, 0]) & 
                            np.isfinite(impdat[:, 1]))
                
                if not np.any(valid_mask):
                    continue
                    
                impdat = impdat[valid_mask]
                
                if impdat.size > 0:
                    data_list.append(impdat)
                    
            except Exception as e:
                print(f"Could not import scan {s}. Error: {e}")
                continue

        print(f"Successfully processed {len(data_list)} scans")
        return data_list

    
    def on_select_span(self, tmin, tmax):
        self.t_min.delete(0, tk.END)
        self.t_min.insert(0, f"{tmin:.4f}")
        self.t_max.delete(0, tk.END)
        self.t_max.insert(0, f"{tmax:.4f}")
        self.status_var.set(f"Selected time range: {tmin:.4f} – {tmax:.4f}")  # ✅ Fixed message
        self.data_extractor(tmin, tmax)

    def on_select_span2(self, x_min, x_max):
         if self.ms_spectrum is None:
              return 
         
         self.x_min.delete(0, tk.END)
         self.x_min.insert(0, f"{x_min:.4f}")
         self.x_max.delete(0, tk.END)
         self.x_max.insert(0, f"{x_max:.4f}")

         mask = (self.ms_spectrum[:,0] >= x_min) & (self.ms_spectrum[:,0] <= x_max)
         self.filtered_ms_data = self.ms_spectrum[mask]

         self.m_z = self.filtered_ms_data[:, 0]
         self.ms_intensity = self.filtered_ms_data[:, 1]

         self.status_var.set(f"Analysis: m/z range set to {x_min:.4f} - {x_max:.4f}")

    def reset_analysis_state(self):
        """Clears all analysis results when a new file is loaded."""
        self.ms_spectrum = None
        self.m_z = None
        self.ms_intensity = None
        self.filtered_ms_data = None
        self.charge_analysis_data = []

        self.axis_right_tab2.clear()
        self.axis_right_tab2.set_xlabel("m/z")
        self.axis_right_tab2.set_ylabel("Intensity")
        self.canvas_right_tab2.draw()

        self.charge_table.delete(*self.charge_table.get_children())
        self.step_display2.delete(1.0, tk.END)

    def set_time_range(self):
        try:
            t_min = float(self.t_min.get())
            t_max = float(self.t_max.get())

            if t_min >= t_max:
                raise ValueError("Minimum time must be less than maximum time")

            # Keep original data untouched
            self.plot_data_on1()

            # Zoom x-axis only
            self.axis_right_tab1.set_xlim(t_min, t_max)

            self.canvas_right_tab1.draw()

            # Show filtered values in textbox
            mask = (self.user_data[:, 0] >= t_min) & (self.user_data[:, 0] <= t_max)
            filtered = self.user_data[mask]

            self.step_display2.delete(1.0, tk.END)
            self.step_display2.insert(tk.END, "TIC DATA\ntime\tIntensity\n")

            for row in filtered:
                self.step_display2.insert(
                    tk.END,
                    f"{row[0]:.4f}\t{row[1]:.4f}\n"
                )

            self.status_var.set("Status: Time range set successfully")

        except ValueError as ve:
            messagebox.showerror("Error", f"Invalid time input: {ve}")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to set time range: {e}")


    # def set_time_range(self):
    #     try:
    #         t_min = float(self.t_min.get())
    #         t_max = float(self.t_max.get())
    #         if t_min >= t_max:
    #             raise ValueError("Minimum time must be less than maximum time")

    #         mask = (self.user_data[:, 0] >= t_min) & (self.user_data[:, 0] <= t_max)
    #         filtered = self.user_data[mask]
    #         self.time = filtered[:, 0]
    #         self.intensity = filtered[:, 1]

    #         self.plot_data_on1()

    #         self.step_display2.delete(1.0, tk.END)
    #         self.step_display2.insert(tk.END, "TIC DATA\ntime\tIntensity\n")
    #         for row in filtered:
    #             self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")

    #         self.status_var.set("Status: Time range set successfully")

    #     except ValueError as ve:
    #         messagebox.showerror("Error", f"Invalid time input: {ve}")
    #     except Exception as e:
    #         messagebox.showerror("Error", f"Failed to set time range: {e}")

    # def set_time_range(self):
    #     try:
    #         t_min = float(self.t_min.get())
    #         t_max = float(self.t_max.get())
    #         if t_min >= t_max:
    #                raise ValueError("Minimum time must be less than maximum time")
              
    #         mask = (self.user_data[:, 0] >= t_min) & (self.user_data[:, 0] <= t_max)
    #         self.user_data = self.user_data[mask]
    #         self.time = self.user_data[:, 0]
    #         self.intensity = self.user_data[:, 1]

    #         self.plot_data_on1()
    #         self.step_display2.delete(1.0, tk.END)
    #         self.step_display2.insert(tk.END, "TIC DATA\ntime\tIntensity\n")
    #         for row in self.user_data:
    #             self.step_display2.insert(tk.END, f"{row[0]:.4f}\t{row[1]:.4f}\n")

    #         self.status_var.set("Status: Time range set successfully")

    #     except ValueError as ve:
    #         messagebox.showerror("Error", f"Invalid time input: {ve}")
    #     except Exception as e:
    #         messagebox.showerror("Error", f"Failed to set time range: {e}")

    def reset_time_range(self):

        if self.user_data is None:
            return

        self.axis_right_tab1.set_xlim(
            np.min(self.user_data[:,0]),
            np.max(self.user_data[:,0])
        )

        self.canvas_right_tab1.draw()

        self.status_var.set("Status: Time range reset")

    def update_clock(self):
            # current_time = time.strftime("%H:%M:%S")
            # self.clock_var.set(f"Clock: {current_time}")

            elapsed_time = time.time() - self.gui_start_time
            self.up_time_of_gui_var.set(f"Up Time: {elapsed_time:.2f} seconds")

            self.root.after(1000, self.update_clock)

if __name__ == "__main__":
    
    root = tk.Tk()    
    app = Deconvoluter(root)
    root.mainloop()
