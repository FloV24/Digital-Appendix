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

def rotation_matrix_from_vectors(v1, v2):
    """
    Berechnet die Rotationsmatrix, die den Vektor v1 auf den Vektor v2 abbildet.
    """
    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    if np.isclose(dot, 1.0):  # Kein Unterschied zwischen v1 und v2
        return np.eye(3)
    if np.isclose(dot, -1.0):  # v1 und v2 sind entgegengesetzt
        return -np.eye(3)
    cross_matrix = np.array([
        [0, -cross[2], cross[1]],
        [cross[2], 0, -cross[0]],
        [-cross[1], cross[0], 0]
    ])
    return np.eye(3) + cross_matrix + cross_matrix @ cross_matrix * (1 / (1 + dot))

def transform_ligand(atoms, target_position, distance=2.0):
    """
    Transformiert die Ligandenkoordinaten, um sie an die Zielposition anzupassen.
    Das Atom des Liganden, das ursprünglich im Ursprung liegt, wird auf die Zielposition verschoben.
    Der Ligand wird zusätzlich rotiert, um korrekt ausgerichtet zu sein.
    """
    # Finde das Atom, das ursprünglich im Ursprung liegt (Koordinaten nahe [0, 0, 0])
    central_atom_coord = None
    for atom, x, y, z in atoms:
        if np.isclose([x, y, z], [0.0, 0.0, 0.0]).all():
            central_atom_coord = np.array([x, y, z])
            break

    if central_atom_coord is None:
        raise ValueError("Kein Atom des Liganden liegt im Ursprung.")

    # Zielrichtung (vom Ursprung zur Zielposition)
    target_direction = target_position / np.linalg.norm(target_position)

    # Ursprüngliche Richtung des Liganden (positive z-Richtung)
    original_direction = np.array([0.0, 0.0, 1.0])

    # Rotationsmatrix berechnen
    rotation_matrix = rotation_matrix_from_vectors(original_direction, target_direction)

    transformed_atoms = []

    for atom, x, y, z in atoms:
        # Verschiebe relativ zum zentralen Atom
        coord = np.array([x, y, z]) - central_atom_coord

        # Rotiere die Koordinaten
        coord = rotation_matrix @ coord

        # Verschiebe den Liganden zur Zielposition (nur Translation, keine Skalierung)
        coord += target_position - (rotation_matrix @ central_atom_coord)

        transformed_atoms.append((atom, coord[0], coord[1], coord[2]))

    return transformed_atoms


def build_tetrahedral_complex(db_path, central_atom, ligands):
    """
    Erstellt einen tetraedrischen Komplex um das Zentralatom mit einer Distanz von 1.5.
    """
    # Tetraedrische Positionen relativ zum Ursprung
    positions = [
        np.array([1.0, 1.0, 1.0]),   # Ecke 1
        np.array([-1.0, -1.0, 1.0]), # Ecke 2
        np.array([-1.0, 1.0, -1.0]), # Ecke 3
        np.array([1.0, -1.0, -1.0])  # Ecke 4
    ]
    positions = [pos / np.linalg.norm(pos) * 1.5 for pos in positions]  # Normiere auf Abstand 1.5

    complex_atoms = [(central_atom, 0.0, 0.0, 0.0)]  # Zentralatom im Ursprung

    # Transformiere und füge jeden Liganden hinzu
    for ligand_name, position in zip(ligands, positions):
        ligand_atoms = fetch_atoms_for_ligand(db_path, ligand_name)
        transformed_ligand = transform_ligand(ligand_atoms, position)
        complex_atoms.extend(transformed_ligand)

    return complex_atoms

def save_all_tetrahedral_arrangements(db_path, central_atom, ox, ligands, output_dir, spin=None, chrg=0, mult=1):
    # Nur ein Arrangement mit den übergebenen Liganden
    arrangement = ligands

    # Erstelle einen Unterordner für das Zentralatom
    central_atom_dir = os.path.join(output_dir, central_atom)
    os.makedirs(central_atom_dir, exist_ok=True)

    # Erstelle den Ordner für den spezifischen Komplex
    ligands_str = "_".join(arrangement)
    complex_dir = os.path.join(central_atom_dir, f"{central_atom}_{ox}_{ligands_str}")
    os.makedirs(complex_dir, exist_ok=True)

    # Spin-Ergänzung
    spin_suffix = ""
    if spin == "high":
        spin_suffix = "_highS"
    elif spin == "low":
        spin_suffix = "_lowS"

    # Dateiname ohne 1-1 usw.
    file_name = f"{central_atom}_{ox}_{ligands_str}{spin_suffix}.inp"
    inp_file_path = os.path.join(complex_dir, file_name)
    
    # Baue den Komplex
    complex_atoms = build_tetrahedral_complex(db_path, central_atom, arrangement)

    # Schreibe die .inp-Datei
    with open(inp_file_path, 'w') as inp_file:
        inp_file.write("!XTB VERYTIGHTSCF LooseOpt\n%geom MaxIter 500 end\n")
        inp_file.write(f"* xyz {chrg} {mult}\n")
        for atom, x, y, z in complex_atoms:
            inp_file.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
        inp_file.write("*\n")

if __name__ == "__main__":
    # Verzeichnis des Programms
    program_dir = os.path.dirname(os.path.abspath(__file__))  # Verzeichnis des Programms
    parent_dir = os.path.dirname(program_dir)  # Übergeordnetes Verzeichnis

    # Standardpfade relativ zum übergeordneten Verzeichnis
    default_db_path = os.path.join(parent_dir, "Ligand", "liganden.db")  # Datenbank im Unterordner "Ligand"
    default_output_dir = os.path.join(parent_dir, "Komplexe", "Tetraedrisch")  # Ordner im übergeordneten Verzeichnis

    # Argumentparser einrichten
    parser = argparse.ArgumentParser(description="Generiert tetraedrische Komplexe mit Liganden.")
    parser.add_argument("--db_path", default=default_db_path, help=f"Pfad zur SQLite-Datenbank (Standard: {default_db_path}).")
    parser.add_argument("--central_atom", required=True, help="Zentralatom des Komplexes (z. B. Fe).")
    parser.add_argument("--ox", type=int, required=True, help="Oxidationsstufe des Zentralatoms (Ganzzahl).")
    parser.add_argument("--ligands", required=True, help="Liste der Liganden, durch Kommas getrennt (z. B. Water,Ammonia,Chloride).")
    parser.add_argument("--output_dir", default=default_output_dir, help=f"Verzeichnis für die Ausgabe der XYZ-Dateien (Standard: {default_output_dir}).")
    parser.add_argument("--chrg", type=int, required=True, help="Ladung des Komplexes (Ganzzahl).")
    parser.add_argument("--mult", type=int, required=True, help="Multiplizität des Komplexes (Ganzzahl).")
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

    # Überprüfen, ob die Datenbank existiert
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Die Datenbank '{db_path}' wurde nicht gefunden.")

    # Überprüfen, ob das Ausgabe-Verzeichnis existiert, und erstellen, falls nicht
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Abrufen aller Liganden aus der Datenbank
    all_ligands = fetch_all_ligands(db_path)

    # Überprüfen, ob alle angegebenen Liganden in der Datenbank existieren
    for name in selected_names:
        if name not in all_ligands:
            raise ValueError(f"Ligand '{name}' wurde nicht in der Datenbank gefunden.")

    # Generiere und speichere alle möglichen Anordnungen
    if len(selected_names) == 4:
        save_all_tetrahedral_arrangements(db_path, central_atom, ox, selected_names, output_dir, spin, chrg, mult)
    else:
        raise ValueError("Es müssen genau 4 Liganden angegeben werden.")