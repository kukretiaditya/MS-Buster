import tkinter as tk
from tkinter import messagebox
import sys
import subprocess
from PIL import Image, ImageTk

def open_deconvolution():
    subprocess.Popen([sys.executable, "Deconvoluter.py"])

def open_peptide_mapping():
    subprocess.Popen([sys.executable, "peptidemapping.py"])

def open_glycan_mapping():
    subprocess.Popen([sys.executable, "glycan.py"])

root = tk.Tk()
root.title("MSbuster")
root.iconbitmap(r"cbt.ico")

screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
window_width = int(screen_width * 0.3)
window_height = int(screen_height * 0.3)

center_x = int((screen_width - window_width) / 2)
center_y = int((screen_height - window_height) / 2)
root.geometry(f"{window_width}x{window_height}+{center_x}+{center_y}")
root.configure(bg="white")

top_y = int(window_height * 0.08)
logo_y = int(window_height * 0.18)
bottom_y = int(window_height * 0.8)
center_x = window_width // 2
left_x = int(window_width * 0.04)
right_x = int(window_width * 0.65)

logo_width = 180
logo_height = 180

logo_label = tk.Label(root, bg="white")

logo_y = top_y + (bottom_y - top_y) // 1.7 - logo_height // 2
logo_label.place(x=center_x - logo_width // 1.7, y=logo_y)

gif_path = r"logo.gif"
gif = Image.open(gif_path)
frames = []

try:
    while True:
        frame = gif.copy().convert("RGBA").resize((logo_width, logo_height), Image.LANCZOS)
        frames.append(ImageTk.PhotoImage(frame))
        gif.seek(len(frames))
except EOFError:
    pass  # End of GIF

# Animate once, then hold final frame
def animate_once(ind=0):
    if ind < len(frames):
        logo_label.configure(image=frames[ind])
        root.after(100, animate_once, ind + 1)
    else:
        logo_label.configure(image=frames[-1])

animate_once()

btn_peptide = tk.Button(root, text="Deconvolution", command=open_deconvolution,
                        bg="#d0f0c0", width=15, height=1)
btn_peptide.place(x=center_x - 75, y=top_y)

btn_deconv = tk.Button(root, text="Peptide mapping", command=open_peptide_mapping,
                       bg="#fdd9b5", width=15, height=1)
btn_deconv.place(x=left_x, y=bottom_y)

btn_glycan = tk.Button(root, text="Glycan analysis", command=open_glycan_mapping,
                       bg="#b0e0e6", width=15, height=1)
btn_glycan.place(x=right_x, y=bottom_y)

root.mainloop()