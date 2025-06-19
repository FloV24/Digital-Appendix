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
"""
# Only needed if you want to generate isomers for each ligand combinatio


def unique_permutations(ligands):
    # Alle Permutationen der Liganden
    all_permutations = set(permutations(ligands))
    symmetries = [
        # Identität
        (0, 1, 2, 3, 4, 5),  # 1: Identität (keine Veränderung)
        # 90° Rotationen C4
        (0, 1, 4, 5, 3, 2),  # 3: 90° Rotation um die x-Achse
        (5, 4, 2, 3, 0, 1),  # 4: 90° Rotation um die y-Achse
        (2, 3, 1, 0, 4, 5),  # 2: 90° Rotation um die z-Achse
        # 180° Rotationen C4
        (0, 1, 3, 2, 5, 4),  # 5: 180° Rotation um die x-Achse
        (1, 0, 2, 3, 5, 4),  # 6: 180° Rotation um die y-Achse
        (1, 0, 3, 2, 4, 5),  # 7: 180° Rotation um die z-Achse
        # 270° Rotationen C4
        (0, 1, 5, 4, 2, 3),  # 8: 270° Rotation um die x-Achse
        (4, 5, 2, 3, 1, 0),  # 9: 270° Rotation um die y-Achse
        (3, 2, 0, 1, 4, 5),  # 10: 270° Rotation um die z-Achse
        # 120° Rotationen C3 (Angabe von Ausgangsquadrant in -z)
        (4, 5, 1, 0, 3, 2),  # 11: 120° Rotation um die -x/+y-Achse
        (2, 3, 5, 4, 1, 0),  # 12: 120° Rotation um die +x/+y-Achse
        (2, 3, 4, 5, 0, 1),  # 13: 120° Rotation um die -x/-y-Achse
        (5, 4, 1, 0, 2, 3),  # 14: 120° Rotation um die +x/-y-Achse 
        # 240° Rotationen C3 (Angabe von Ausgangsquadrant in -z)
        (3, 2, 5, 4, 0, 1),  # 15: 240° Rotation um die -x/+y-Achse
        (5, 4, 0, 1, 3, 2),  # 16: 240° Rotation um die +x/+y-Achse
        (4, 5, 0, 1, 2, 3),  # 17: 240° Rotation um die -x/-y-Achse
        (3, 2, 4, 5, 1, 0),  # 18: 240° Rotation um die +x/-y-Achse
        # 180° Rotationen C2
        (2, 3, 0, 1, 5, 4),  # 19: 180° Rotation um die xy-Ebene Diagonal nach +x/+y
        (3, 2, 1, 0, 5, 4),  # 20: 180° Rotation um die xy-Ebene Diagonal nach -x/+y
        (5, 4, 3, 2, 1, 0),  # 21: 180° Rotation um die xz-Ebene Diagonal nach +x/+z
        (4, 5, 3, 2, 0, 1),  # 22: 180° Rotation um die xz-Ebene Diagonal nach -x/+z
        (1, 0, 4, 5, 2, 3),  # 23: 180° Rotation um die yz-Ebene Diagonal nach +y/+z
        (1, 0, 5, 4, 3, 2),  # 24: 180° Rotation um die yz-Ebene Diagonal nach -y/+z
        # Inversion
        (1, 0, 3, 2, 5, 4),  # 25: Inversion (alle Atome werden gespiegelt)
        # Spiegelungen Achsen sigma_h
        (0, 1, 2, 3, 5, 4),  # 26: Spiegelung an der xy-Ebene
        (0, 1, 3, 2, 4, 5),  # 27: Spiegelung an der xz-Ebene
        (1, 0, 2, 3, 4, 5),  # 28: Spiegelung an der yz-Ebene
        # Spiegelungen Diagonalen sigma_d
        (3, 2, 1, 0, 4, 5),  # 29: Spiegelung an der -xy / z Ebene
        (2, 3, 0, 1, 4, 5),  # 30: Spiegelung an der +xy / z Ebene
        (5, 4, 2, 3, 1, 0),  # 31: Spiegelung an der -xz / y Ebene
        (4, 5, 2, 3, 0, 1),  # 32: Spiegelung an der +xz / y Ebene
        (0, 1, 5, 4, 3, 2),  # 33: Spiegelung an der -yz / x Ebene
        (0, 1, 4, 5, 2, 3),  # 34: Spiegelung an der +yz / x Ebene
        # 6 S4 Rotationen 90° und 270°
        (1, 0, 4, 5, 3, 2),  # 35: 90° Rotation um die x-Achse + Spiegel
        (5, 4, 3, 2, 0, 1),  # 36: 90° Rotation um die y-Achse + Spiegel
        (2, 3, 1, 0, 5, 4),  # 37: 90° Rotation um die z-Achse + Spiegel
        (1, 0, 5, 4, 2, 3),  # 38: 270° Rotation um die x-Achse + Spiegel
        (4, 5, 3, 2, 1, 0),  # 39: 270° Rotation um die y-Achse + Spiegel
        (3, 2, 0, 1, 5, 4),  # 40: 270° Rotation um die z-Achse + Spiegel
        # 4 S6 Rotationen 60° (Angabe von Ausgangsquadrant in -z)
        (2, 3, 4, 5, 1, 0),  # 41: S6 Rotation in -x/+y Achse
        (4, 5, 1, 0, 2, 5),  # 42: S6 Rotation in +x/+y Achse
        (5, 4, 1, 0, 3, 2),  # 43: S6 Rotation in -x/-y Achse
        (2, 3, 5, 4, 0, 1),  # 44: S6 Rotation in +x/-y Achse
        # 4 S6 Rotationen 180° (Angabe von Ausgangsquadrant in -z)
        (4, 5, 0, 1, 2, 3),  # 45: S6 Rotation in -x/+y Achse
        (3, 4, 2, 5, 0, 1),  # 46: S6 Rotation in +x/+y Achse
        (3, 2, 5, 4, 1, 0),  # 47: S6 Rotation in -x/-y Achse
        (4, 5, 0, 1, 3, 2),  # 48: S6 Rotation in +x/-y Achse
    ]

    unique_structures = set()

    for perm in all_permutations:
        # Schritt 1: Prüfe, ob die aktuelle Permutation durch Symmetrieoperationen abgedeckt ist
        is_unique = True
        for stored_perm in unique_structures:
            for sym in symmetries:
                # Wende die Symmetrieoperation auf die gespeicherte Permutation an
                transformed_perm = tuple(stored_perm[i] for i in sym)
                if transformed_perm == perm:
                    is_unique = False
                    break
            if not is_unique:
                break

        # Schritt 2: Wenn die Permutation einzigartig ist, füge sie hinzu
        if is_unique:
            unique_structures.add(perm)

    return unique_structures
"""
def build_octahedral_complex(db_path, central_atom, ligands):
    """
    Erstellt einen oktaedrischen Komplex um das Zentralatom mit einer Distanz von 1.5.
    """
    # Oktaedrische Positionen relativ zum Ursprung
    positions = [
        np.array([1.5, 0.0, 0.0]),   # +x
        np.array([-1.5, 0.0, 0.0]),  # -x
        np.array([0.0, 1.5, 0.0]),   # +y
        np.array([0.0, -1.5, 0.0]),  # -y
        np.array([0.0, 0.0, 1.5]),   # +z
        np.array([0.0, 0.0, -1.5])   # -z
    ]

    complex_atoms = [(central_atom, 0.0, 0.0, 0.0)]  # Zentralatom im Ursprung

    # Transformiere und füge jeden Liganden hinzu
    for ligand_name, position in zip(ligands, positions):
        ligand_atoms = fetch_atoms_for_ligand(db_path, ligand_name)
        transformed_ligand = transform_ligand(ligand_atoms, position)
        complex_atoms.extend(transformed_ligand)

    return complex_atoms

def save_all_octahedral_arrangements(db_path, central_atom, ox, ligands, output_dir, spin=None, chrg=0, mult=1):
    """
    Speichert nur die direkt übergebenen Liganden als einen Komplex,
    ohne Permutationen zu bilden oder 1-1 Suffixe zu verwenden.
    """
    # Erstelle einen Unterordner für das Zentralatom
    central_atom_dir = os.path.join(output_dir, central_atom)
    os.makedirs(central_atom_dir, exist_ok=True)

    # Erstelle den Ordner für den spezifischen Komplex
    ligands_str = "_".join(ligands)
    complex_dir = os.path.join(central_atom_dir, f"{central_atom}_{ox}_{ligands_str}")
    os.makedirs(complex_dir, exist_ok=True)

    # Spin-Angabe
    spin_suffix = ""
    if spin == "high":
        spin_suffix = "_highS"
    elif spin == "low":
        spin_suffix = "_lowS"

    # Nur ein einziges Arrangement speichern (direkt ligands)
    arrangement_str = "_".join(ligands)
    file_name = f"{central_atom}_{ox}_{arrangement_str}{spin_suffix}.inp"
    inp_file_path = os.path.join(complex_dir, file_name)

    # Komplex erstellen
    complex_atoms = build_octahedral_complex(db_path, central_atom, ligands)

    # .inp file writing
    with open(inp_file_path, 'w') as inp_file:
        inp_file.write("!XTB VERYTIGHTSCF LooseOpt\n%geom MaxIter 500 end\n")
        inp_file.write(f"* xyz {chrg} {mult}\n")
        for atom, x, y, z in complex_atoms:
            inp_file.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")
        inp_file.write("*\n")

if __name__ == "__main__":

    program_dir = os.path.dirname(os.path.abspath(__file__))
    parent_dir = os.path.dirname(program_dir)

    parser = argparse.ArgumentParser(description="Generiert oktaedrische Komplexe mit Liganden.")
    parser.add_argument("--db_path", required=True, help="Pfad zur SQLite-Datenbank.")
    parser.add_argument("--central_atom", required=True, help="Zentralatom des Komplexes (z. B. Fe).")
    parser.add_argument("--ox", type=int, required=True, help="Oxidationsstufe des Zentralatoms (Ganzzahl).")
    parser.add_argument("--ligands", required=True, help="Liste der Liganden, durch Kommas getrennt (z. B. Water,Ammonia,Chloride).")
    parser.add_argument("--output_dir", required=True, help="Verzeichnis für die Ausgabe der XYZ-Dateien.")
    parser.add_argument("--chrg", type=int, required=True, help="Ladung des Komplexes (Ganzzahl).")
    parser.add_argument("--mult", type=int, required=True, help="Multiplizität des Komplexes (Ganzzahl).")
    parser.add_argument("--spin", choices=["high", "low"], help="Spin-Zustand des Komplexes (high oder low).")

    args = parser.parse_args()

    db_path = args.db_path
    central_atom = args.central_atom
    ox = args.ox
    selected_names = [name.strip() for name in args.ligands.split(",")]
    output_dir = args.output_dir
    chrg = args.chrg
    mult = args.mult
    spin = args.spin

    if len(selected_names) != 6:
        raise ValueError("Es müssen genau 6 Liganden angegeben werden.")

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Die Datenbank '{db_path}' wurde nicht gefunden.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_ligands = fetch_all_ligands(db_path)
    for name in selected_names:
        if name not in all_ligands:
            raise ValueError(f"Ligand '{name}' wurde nicht in der Datenbank gefunden.")

    save_all_octahedral_arrangements(db_path, central_atom, ox, selected_names, output_dir, spin, chrg, mult)