import sqlite3
import numpy as np
from itertools import permutations
import os
import argparse

def fetch_all_ligands(db_path):
    """
    Ruft alle Liganden aus der neuen Datenbank ab.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Abrufen aller Liganden
    cursor.execute("SELECT name FROM liganden")
    ligands = [row[0] for row in cursor.fetchall()]

    conn.close()
    return ligands

def fetch_atoms_for_ligand(db_path, ligand_name):
    """
    Ruft die XYZ-Daten eines Liganden aus der neuen Datenbank ab.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Abrufen der XYZ-Daten
    cursor.execute("SELECT xyz_daten FROM liganden WHERE name = ?", (ligand_name,))
    result = cursor.fetchone()
    if result is None:
        raise ValueError(f"Ligand '{ligand_name}' wurde nicht in der Datenbank gefunden.")
    
    # XYZ-Daten parsen
    xyz_data = result[0].strip().split("\n")
    atoms = []
    for line in xyz_data:
        parts = line.split()
        atom = parts[0]
        x, y, z = map(float, parts[1:])
        atoms.append((atom, x, y, z))

    conn.close()
    return atoms

def transform_ligand(atoms, target_position):
    """
    Transformiert die Ligandenkoordinaten, um sie an die Zielposition anzupassen.
    """
    central_atom_coord = np.array([0.0, 0.0, 0.0])  # Ursprung
    transformed_atoms = []

    for atom, x, y, z in atoms:
        coord = np.array([x, y, z]) - central_atom_coord
        coord += target_position
        transformed_atoms.append((atom, coord[0], coord[1], coord[2]))

    return transformed_atoms

def build_trigonal_pyramidal_complex(db_path, central_atom, ligands):
    """
    Erstellt einen trigonal-pyramidalen Komplex um das Zentralatom mit Koordinationszahl 3.
    """
    # Trigonal-pyramidale Positionen relativ zum Ursprung für Koordinationszahl 3
    positions = [
        np.array([1.0, 0.0, 0.0]),   # +x
        np.array([-0.5, 0.866, 0.0]),  # -x, +y
        np.array([-0.5, -0.866, 0.0])  # -x, -y
    ]

    complex_atoms = [(central_atom, 0.0, 0.0, 0.0)]  # Zentralatom im Ursprung

    # Transformiere und füge jeden Liganden hinzu
    for ligand_name, position in zip(ligands, positions):
        ligand_atoms = fetch_atoms_for_ligand(db_path, ligand_name)
        transformed_ligand = transform_ligand(ligand_atoms, position)
        complex_atoms.extend(transformed_ligand)

    return complex_atoms

def save_trigonal_pyramidal_arrangement(db_path, central_atom, ox, ligands, output_dir, spin=None, chrg=0, mult=1):
    """
    Generiert einen einzelnen trigonal-pyramidalen Komplex mit den angegebenen Liganden und 
    speichert ihn als .inp-Datei.
    """
    # Bestimme den Spin-Ergänzung für den Dateinamen
    spin_suffix = ""
    if spin == "high":
        spin_suffix = "_highS"
    elif spin == "low":
        spin_suffix = "_lowS"

    # Erstelle den Ordner für das Zentralatom innerhalb des spezifischen Unterordners
    central_atom_dir = os.path.join(output_dir, central_atom)
    os.makedirs(central_atom_dir, exist_ok=True)

    # Erstelle den Ordner für den spezifischen Komplex
    ligands_str = "_".join(ligands)
    complex_dir = os.path.join(central_atom_dir, f"{central_atom}_{ox}_{ligands_str}")
    os.makedirs(complex_dir, exist_ok=True)
    
    # Erzeuge einen einfachen Dateinamen ohne Permutationsinformationen
    file_name = f"{central_atom}_{ox}_{ligands_str}{spin_suffix}.inp"
    inp_file_path = os.path.join(complex_dir, file_name)

    # Baue den Komplex mit der gegebenen Ligandenanordnung
    complex_atoms = build_trigonal_pyramidal_complex(db_path, central_atom, ligands)

    # Schreibe die .inp-Datei
    with open(inp_file_path, 'w') as inp_file:
        inp_file.write("!XTB VERYTIGHTSCF LooseOpt\n%geom MaxIter 500 end\n")
        inp_file.write(f"* xyz {chrg} {mult}\n")
        for atom, x, y, z in complex_atoms:
            inp_file.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
        inp_file.write("*\n")

if __name__ == "__main__":
    # Argumentparser einrichten
    parser = argparse.ArgumentParser(description="Generiert trigonal-pyramidale Komplexe mit Liganden.")
    parser.add_argument("--db_path", required=True, help="Pfad zur SQLite-Datenbank.")
    parser.add_argument("--central_atom", required=True, help="Zentralatom des Komplexes (z. B. Fe).")
    parser.add_argument("--ligands", required=True, help="Liste der Liganden, durch Kommas getrennt (z. B. CO,H2O,Cl).")
    parser.add_argument("--output_dir", required=True, help="Verzeichnis für die Ausgabe der .inp-Dateien.")
    parser.add_argument("--chrg", type=int, required=True, help="Ladung des Komplexes (Ganzzahl).")
    parser.add_argument("--mult", type=int, required=True, help="Multiplizität des Komplexes (Ganzzahl).")
    parser.add_argument("--ox", type=int, required=True, help="Oxidationsstufe des Zentralatoms (Ganzzahl).")
    parser.add_argument("--spin", choices=["high", "low"], help="Spin-Zustand des Komplexes (high oder low).")
    
    args = parser.parse_args()

    # Eingaben aus den Argumenten
    db_path = args.db_path
    central_atom = args.central_atom
    ox = args.ox
    selected_names = [name.strip() for name in args.ligands.split(",")]
    output_dir = args.output_dir
    chrg = args.chrg
    mult = args.mult
    spin = args.spin

    # Überprüfen, ob die Anzahl der Liganden korrekt ist (immer 3)
    if len(selected_names) != 3:
        raise ValueError(f"Es müssen genau 3 Liganden angegeben werden.")

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Die Datenbank '{db_path}' wurde nicht gefunden.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_ligands = fetch_all_ligands(db_path)

    for name in selected_names:
        if name not in all_ligands:
            raise ValueError(f"Ligand '{name}' wurde nicht in der Datenbank gefunden.")

    save_trigonal_pyramidal_arrangement(db_path, central_atom, ox, selected_names, output_dir, spin, chrg, mult)