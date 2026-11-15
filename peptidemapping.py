import tkinter as tk 
from tkinter import filedialog, messagebox, ttk
import numpy as np 
import os 
import matplotlib.pyplot as plt 
from PIL import Image, ImageTk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg,NavigationToolbar2Tk
import time 
import threading 
from matplotlib.figure import Figure
import shutil
from matplotlib.widgets import SpanSelector
from copy import deepcopy
from bisect import bisect_left
from scipy.interpolate import interp1d
import pymzml
from pymzml.utils.utils import index_gzip
import pymzml.obo


class PeptideMappingGUI():
    def __init__(self, root):
        self.root = root
        self.root.state('zoomed')  # Maximize window (Windows only)
        self.root.title("Peptide Mapping")
        self.root.iconbitmap("cbt.ico")

        self.gui_start_time = time.time()
        self.file_path = None
        self.user_data = None 
        self.m_z = None
        self.intensity = None
        self.ms_intensity = None  # Intensity specifically for mass spectrum

        self.x_min = None
        self.x_max = None
        self.history = None
        self.ms_spectrum = None
        self.filtered_ms_data = None

        self.time = None
        self.t_min = None
        self.t_max = None

        self.paned = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.paned.pack(fill='both', expand=True)

        self.status_frame = tk.Frame(self.root, bd = 1, relief=tk.RAISED)
        self.status_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.status_var = tk.StringVar()
        self.status_var.set("Status: Ready")

        self.status_bar = tk.Label(self.status_frame, textvariable=self.status_var, bd=1, relief=tk.RAISED, anchor=tk.W)
        self.status_bar.pack(side=tk.LEFT, fill=tk.X)
        


        
        
if __name__ == "__main__":
    root = tk.Tk()
    app = PeptideMappingGUI(root)
    root.mainloop()

