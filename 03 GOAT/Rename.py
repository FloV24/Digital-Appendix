import os
import re
import shutil
import sys
from collections import Counter

def gruppiere_liganden(liganden):
    counter = Counter(liganden)
    parts = [
        f"({lig}){counter[lig]}" if counter[lig] > 1 else f"({lig})"
        for lig in sorted(counter)
    ]
    return ''.join(parts)

def rename_and_copy_xyz_files(src_roots, dst_root):
    # Neues Namensmuster
    pattern = re.compile(
        r"^([A-Za-z]+)_([A-Za-z]+)_(\d+)_((?:[A-Za-z0-9]+_?)+?)_Spin_(\d+)\.finalensemble\.xyz$"
    )
    count = 0
    spin_files = {}
    for src_root in src_roots:
        for dirpath, _, filenames in os.walk(src_root):
            for filename in filenames:
                if filename.endswith(".finalensemble.xyz"):
                    match = pattern.match(filename)
                    if not match:
                        print(f"Überspringe ungültigen Dateinamen: {filename}")
                        continue
                    konf_idx = match.group(1)  # GeometrieIndex (Buchstabenkombination)
                    metall = match.group(2)  # Metall
                    ox_stufe = match.group(3)  # Oxidationszahl
                    liganden_raw = match.group(4)  # Ligandenliste
                    mult = match.group(5)  # Multiplizität
                    liganden = liganden_raw.split('_')
                    gruppiert = gruppiere_liganden(liganden)
                    basename = f"{konf_idx}_{metall}_{ox_stufe}_{gruppiert}_Spin_{mult}"
                    if basename not in spin_files:
                        spin_files[basename] = os.path.join(dirpath, filename)
    for basename, src_path in spin_files.items():
        # Remove the _Spin_<mult> part for the subfolder name
        subfolder_name = re.sub(r'_Spin_\d+$', '', basename)
        dst_subfolder = os.path.join(dst_root, subfolder_name)
        os.makedirs(dst_subfolder, exist_ok=True)
        new_name = f"{basename}.xyz"
        dst_path = os.path.join(dst_subfolder, new_name)
        shutil.copy2(src_path, dst_path)
        count += 1
    print(f"{count} Dateien kopiert.")

if __name__ == "__main__":
    quellordner1 = r"D:\GOAT_Out"  # Erster Quellordner
    quellordner2 = r"D:\GOAT_OUT2"  # Zweiter Quellordner
    zielordner = r"D:\GOAT_Sorted"  # Zielordner
    quellordner_liste = [quellordner1, quellordner2]
    rename_and_copy_xyz_files(quellordner_liste, zielordner)