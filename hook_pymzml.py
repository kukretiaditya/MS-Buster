# hook_pymzml.py - ONEFILE MODE
import os
import sys
import shutil

if getattr(sys, 'frozen', False):
    exe_dir = os.path.dirname(sys.executable)
    src_obo = os.path.join(sys._MEIPASS, 'pymzml', 'obo')
    dst_obo = os.path.join(exe_dir, 'obo')

    if os.path.isdir(src_obo) and not os.path.isdir(dst_obo):
        shutil.copytree(src_obo, dst_obo)


    try:
        import gzip
        import pymzml.obo

        def patched_parseOBO(self_obo):
            self_obo._OboTranslator__obo_parsed = True
            obo_root = os.path.join(sys._MEIPASS, 'pymzml', 'obo')
            obo_file = os.path.join(
                obo_root,
                'psi-ms{0}.obo'.format('-' + self_obo.version if self_obo.version else ''),
            )
            if os.path.exists(obo_file):
                pass
            elif os.path.exists(obo_file + '.gz'):
                obo_file = obo_file + '.gz'
            else:
                available = [f for f in os.listdir(obo_root) if f.endswith('.obo.gz')]
                if available:
                    obo_file = os.path.join(obo_root, available[-1])
                else:
                    print(f"OBO file not found: {obo_file}")
                    return

            with gzip.open(obo_file, 'rt', encoding='utf-8') as obo:
                collections = {}
                collect = False
                for line in obo:
                    if line.strip() in ('[Term]', ''):
                        collect = True
                        if not collections:
                            continue
                        self_obo.add(collections)
                        collections = {}
                    else:
                        if line.strip() != '' and collect is True:
                            k = line.find(':')
                            collections[line[:k]] = line[k + 1:].strip()

        pymzml.obo.OboTranslator.parseOBO = patched_parseOBO
        print("OBO patch applied successfully")

    except Exception as e:
        print(f"OBO patch failed: {e}")