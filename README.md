<p align="center">
  <img src="header.png" alt="MS-Buster" width="500">
</p>

MS-Buster (Mass Spectrometry Buster) is a Python application for the extraction and deconvolution of intact liquid chromatography-mass spectrometry (LC-MS) data. It provides a graphical interface intended for researchers without programming experience.

<p>
  <img src="https://img.shields.io/badge/Python-3.9-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python 3.9">
  <img src="https://img.shields.io/badge/Platform-Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white" alt="Windows">
  <img src="https://img.shields.io/badge/GUI-Tkinter-FB923C?style=for-the-badge" alt="Tkinter GUI">
  <img src="https://img.shields.io/badge/Input-Agilent%20.d%20%7C%20mzML%20%7C%20text-22D3EE?style=for-the-badge" alt="Agilent .d, mzML, text input">
</p>

## Workflow

1. **Data import.** An Agilent `.d` directory, an mzML file, or a plain-text file is loaded.
2. **Chromatogram display.** Intensity is plotted against retention time.
3. **Peak selection.** A span selector is used to mark the chromatographic peak of interest; its boundaries define the retention-time window.
4. **Spectrum averaging.** All mass spectra within the selected window are extracted and averaged. The averaged spectrum is displayed in a second plot with its own span selector.
5. **Deconvolution.** The intact mass is determined in either manual or automatic mode. All results are written to text files.

## Features

| Feature | Description |
|---|---|
| Agilent `.d` parsing | Reads binary `.d` directories and converts the data into array form |
| Open formats | Reads mzML and plain-text files |
| Interactive chromatogram | Plots intensity against retention time, with peak selection via a matplotlib span selector |
| Spectrum averaging | Averages all spectra within the selected retention-time window into a single m/z spectrum |
| Deconvolution | Provides manual and automatic modes for intact mass determination |
| Result export | Writes all results to text files |
| Graphical interface | Requires no scripting |

## Usage

MS-Buster is distributed as a standalone Windows executable.

1. Download `MS-BUSTER.zip` from the [Releases](https://github.com/kukretiaditya/MS-Buster/releases) page.
2. Extract the entire archive to a local folder. The executable must not be run from within the compressed archive, as it requires the accompanying `_internal/` and `obo/` directories.
3. Open the extracted `MS-BUSTER` folder and launch `MS-BUSTER.exe`.

Proceed as described in the [Workflow](#workflow) section: load the data, select the peak of interest in the chromatogram, select the region of interest in the averaged spectrum, and run manual or automatic deconvolution. Results are saved as text files.
