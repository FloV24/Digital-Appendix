import os
import shutil
import re

def extract_energien_from_xyz(filepath):
    # Bestehende Funktion bleibt unverändert
    energien = []
    with open(filepath, 'r', errors='ignore') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        try:
            n_atoms = int(lines[i].strip())
        except ValueError:
            break
        if i + 1 < len(lines):
            kommentar = lines[i + 1].strip()
            try:
                energy = float(kommentar.split()[0])
                energien.append(energy)
            except Exception:
                pass
        i += n_atoms + 2
    return energien

def cleanup_mult_folders(root_dir, target_dir, unsicher_dir, best_spin_dir):
    Grenzfall_Counter = 0
    Verschoben_Counter = 0
    Einzel_counter = 0
    Gesamt_counter = 0
    
    # Statistik-Variablen
    geometrie_kurz = {"L", "OC", "SP", "SPY", "TS", "T", "TBPY", "TP", "TPR", "TPY"}
    geometrie_stats = {}  # {geometrie: {1: x, 2: y, 3: z}}
    geo_mult_counter = {}  # {geo: {mult: count}}
    komplex_mult_count = {1: 0, 2: 0, 3: 0}  # Zählt Komplexe mit n Multiplizitäten
    
    # Erster Durchlauf: Sammle Statistiken zu Geometrien und Multiplizitäten
    print("Sammle Statistiken...")
    for dirpath, _, filenames in os.walk(root_dir):
        # Finde alle Komplexe in diesem Verzeichnis
        komplexe = {}
        for filename in filenames:
            if filename.endswith(".xyz"):
                parts = filename.split("_Spin_")
                if len(parts) == 2:
                    basename = parts[0]
                    mult = parts[1].replace(".xyz", "")
                    
                    # Extrahiere Geometrie-Kürzel (erster Teil des Basename)
                    match = re.match(r"([A-Za-z]+)_", basename)
                    if match and match.group(1) in geometrie_kurz:
                        geo = match.group(1)
                        
                        # Zähle Multiplizitäten pro Komplex
                        if basename not in komplexe:
                            komplexe[basename] = set()
                        komplexe[basename].add(mult)
                        
                        # Zähle Häufigkeit jeder Multiplizität pro Geometrie
                        if geo not in geo_mult_counter:
                            geo_mult_counter[geo] = {}
                        if mult not in geo_mult_counter[geo]:
                            geo_mult_counter[geo][mult] = 0
                        geo_mult_counter[geo][mult] += 1
        
        # Verarbeite die gesammelten Informationen für dieses Verzeichnis
        for basename, mults in komplexe.items():
            # Extrahiere Geometrie
            match = re.match(r"([A-Za-z]+)_", basename)
            if match and match.group(1) in geometrie_kurz:
                geo = match.group(1)
                
                # Zähle Komplexe mit n Multiplizitäten pro Geometrie
                mult_count = len(mults)
                if geo not in geometrie_stats:
                    geometrie_stats[geo] = {1: 0, 2: 0, 3: 0}
                if mult_count in (1, 2, 3):
                    geometrie_stats[geo][mult_count] += 1
                    komplex_mult_count[mult_count] += 1

    # Zweiter Durchlauf: Führe die eigentliche Bereinigung durch
    print("Bereinige Multiplizitäten...")
    for dirpath, _, filenames in os.walk(root_dir):
        mult_files = {}
        for filename in filenames:
            if filename.endswith(".xyz"):
                # Extrahiere Basename und Multiplizität
                parts = filename.split("_Spin_")
                if len(parts) == 2:
                    basename = parts[0]
                    mult = parts[1].replace(".xyz", "")
                    mult_files.setdefault(basename, []).append((int(mult), os.path.join(dirpath, filename)))

        for basename, files in mult_files.items():
            Gesamt_counter += 1
            files.sort(key=lambda x: x[0])  # Sortiere nach Multiplizität
            energies = [(mult, extract_energien_from_xyz(filepath)) for mult, filepath in files]

            # Vergleich aller Multiplizitäten
            best_mult = None
            for i, (mult_i, energies_i) in enumerate(energies):
                is_best = True  # Reset für jede Multiplizität
                for j, (mult_j, energies_j) in enumerate(energies):
                    if i != j and energies_i and energies_j:
                        # Prüfe, ob alle Energien von mult_i kleiner sind als die von mult_j
                        if not all(e_i < e_j for e_i in energies_i for e_j in energies_j):
                            is_best = False
                            break
                if is_best:
                    best_mult = mult_i
                    break

            if best_mult is not None:
                # Verschiebe den besten Spin-Zustand in den Best-Spin-Ordner
                rel_path = os.path.relpath(dirpath, root_dir)
                best_spin_subfolder = os.path.join(best_spin_dir, rel_path)
                os.makedirs(best_spin_subfolder, exist_ok=True)
                for mult, filepath in files:
                    if mult == best_mult:  # Nur den besten Spin-Zustand verschieben
                        new_name = f"{basename}_Spin_{mult}.xyz"
                        shutil.move(filepath, os.path.join(best_spin_subfolder, new_name))
                        Verschoben_Counter += 1
                    else:
                        # Verschiebe die ungünstigeren Multiplizitäten
                        target_subfolder = os.path.join(target_dir, rel_path)
                        os.makedirs(target_subfolder, exist_ok=True)
                        new_name = f"{basename}_Spin_{mult}.xyz"
                        shutil.move(filepath, os.path.join(target_subfolder, new_name))
                        Verschoben_Counter += 1
            else:
                # Verschiebe den gesamten Subordner in den Unsicherheitsordner
                rel_path = os.path.relpath(dirpath, root_dir)
                target_subfolder = os.path.join(unsicher_dir, rel_path)
                os.makedirs(os.path.dirname(target_subfolder), exist_ok=True)
                shutil.move(dirpath, target_subfolder)
                Grenzfall_Counter += len(files)

        if len(mult_files) == 1:
            Einzel_counter += 1

    # Schreibe Statistik in eine Datei im selben Ordner wie das Skript
    statistik_path = os.path.join(os.path.dirname(__file__), "statistik.txt")
    with open(statistik_path, "w", encoding="utf-8") as f:
        # Allgemeine Statistiken
        f.write("=== Allgemeine Statistiken ===\n")
        f.write(f"Vergleichspaare: {Gesamt_counter}\n")
        f.write(f"Verschobene Dateien: {Verschoben_Counter}\n")
        f.write(f"Grenzfälle: {Grenzfall_Counter}\n")
        f.write(f"Einzelne Multiplizitäten: {Einzel_counter}\n\n")
        
        # Statistik zu Anzahl der Multiplizitäten
        f.write("=== Komplexe mit n Multiplizitäten ===\n")
        f.write(f"Komplexe mit 1 Multiplizität: {komplex_mult_count[1]}\n")
        f.write(f"Komplexe mit 2 Multiplizitäten: {komplex_mult_count[2]}\n")
        f.write(f"Komplexe mit 3 Multiplizitäten: {komplex_mult_count[3]}\n\n")
        
        # Geometrie-Statistiken
        f.write("=== Multiplizitäten pro Geometrie ===\n")
        f.write("Geometrie\tKomplexe_mit_1_Mult\tKomplexe_mit_2_Mult\tKomplexe_mit_3_Mult\tGesamt\n")
        for geo in sorted(geometrie_stats):
            gesamt = geometrie_stats[geo][1] + geometrie_stats[geo][2] + geometrie_stats[geo][3]
            row = [geo, geometrie_stats[geo][1], geometrie_stats[geo][2], geometrie_stats[geo][3], gesamt]
            f.write("\t".join(str(x) for x in row) + "\n")
        f.write("\n")
        
        # Häufigkeit jeder einzelnen Multiplizität pro Geometrie
        f.write("=== Häufigkeit jeder Multiplizität pro Geometrie ===\n")
        f.write("Geometrie\tMultiplizität\tAnzahl\n")
        for geo in sorted(geo_mult_counter):
            for mult in sorted(geo_mult_counter[geo], key=lambda x: int(x)):
                f.write(f"{geo}\t{mult}\t{geo_mult_counter[geo][mult]}\n")

    print(f"Bereinigung abgeschlossen. Alle Statistiken wurden in {statistik_path} gespeichert.")

if __name__ == "__main__":
    # Passe die Pfade ggf. an
    sortierter_ordner = r"D:\GOAT_Sorted"  # Basisordner
    ziel_ordner = r"D:\GOAT_Unguenstig"
    unsicher_ordner = r"D:\GOAT_Spin_Unsicher"
    best_spin_ordner = r"D:\GOAT_Best_Spin"
    cleanup_mult_folders(sortierter_ordner, ziel_ordner, unsicher_ordner, best_spin_ordner)