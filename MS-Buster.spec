# MS-BUSTER.spec - ONEDIR MODE
import os
from PyInstaller.utils.hooks import collect_all, collect_data_files

datas_multiplierz, binaries_multiplierz, hiddenimports_multiplierz = collect_all('multiplierz')
datas_pymzml,      binaries_pymzml,      hiddenimports_pymzml      = collect_all('pymzml')
datas_reportlab,   binaries_reportlab,   hiddenimports_reportlab   = collect_all('reportlab')
datas_scipy,       binaries_scipy,       hiddenimports_scipy       = collect_all('scipy')
datas_matplotlib,  binaries_matplotlib,  hiddenimports_matplotlib  = collect_all('matplotlib')

import pymzml, os
pymzml_dir = os.path.dirname(pymzml.__file__)
obo_dir = os.path.join(pymzml_dir, 'obo')

datas_pymzml_obo = [
    (os.path.join(obo_dir, f), os.path.join('pymzml', 'obo'))
    for f in os.listdir(obo_dir)
    if f.endswith('.obo.gz')
]

a = Analysis(
    ['Deconvoluter.py'],
    pathex=['.'],
    binaries=(
        binaries_multiplierz +
        binaries_pymzml +
        binaries_reportlab +
        binaries_scipy +
        binaries_matplotlib
    ),
    datas=(
        [
            ('mzML.py',               '.'),
            ('launcher.py',           '.'),
            ('peptidemapping.py',     '.'),
            ('cbt.ico',               '.'),
            ('header.png',            '.'),
            ('Untitled design.gif',   '.'),
        ]
        + datas_pymzml_obo
        + datas_multiplierz
        + datas_pymzml
        + datas_reportlab
        + datas_scipy
        + datas_matplotlib
    ),
    hiddenimports=[
        'multiplierz.mzAPI',
        'multiplierz.mzAPI.mzML',
        'multiplierz.mzAPI.mzWiff',
        'multiplierz.mzAPI.mzRaw',
        'multiplierz.mzAPI.mzMGF',
        'multiplierz.mzAPI.management',
        'multiplierz.mzAPI.mzFile',
        'pymzml',
        'pymzml.obo',
        'pymzml.utils',
        'pymzml.utils.utils',
        'scipy.interpolate',
        'scipy.signal',
        'scipy.special._ufuncs_cxx',
        'scipy._lib.messagestream',
        'matplotlib.backends.backend_tkagg',
        'matplotlib.backends.backend_agg',
        'reportlab.platypus',
        'reportlab.lib.pagesizes',
        'reportlab.lib.colors',
        'reportlab.lib.styles',
        'tkinter',
        'tkinter.ttk',
        'tkinter.filedialog',
        'tkinter.messagebox',
        'numpy',
        'numpy.core._multiarray_umath',
        'bisect',
        'copy',
        'threading',
        'io',
        're',
        'shutil',
        'traceback',
    ]
    + hiddenimports_multiplierz
    + hiddenimports_pymzml
    + hiddenimports_reportlab
    + hiddenimports_scipy
    + hiddenimports_matplotlib,
    hookspath=['.'],
    hooksconfig={},
    runtime_hooks=['hook_pymzml.py'],
    excludes=[
        'PyQt5', 'PyQt6', 'PySide2', 'PySide6',
        'IPython', 'jupyter', 'notebook',
        'pytest', 'setuptools',
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='MS-BUSTER',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    icon='cbt.ico',
    uac_admin=True,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='MS-BUSTER',
)