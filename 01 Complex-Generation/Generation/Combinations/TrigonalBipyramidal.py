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

def unique_permutations(ligands):
    """
    Generiert eindeutige Permutationen der Liganden unter Berücksichtigung der D3h-Symmetrie
    des trigonal-bipyramidalen Komplexes.
    """
    # Alle Permutationen der Liganden
    all_permutations = set(permutations(ligands))
    
    # Symmetrieoperationen der D3h-Punktgruppe für trigonal-bipyramidale Geometrie
    # Index 0,1: axiale Positionen (z-Achse)
    # Index 2,3,4: äquatoriale Positionen (xy-Ebene)
    symmetries = [
        # Identität
        (0, 1, 2, 3, 4),  # Identität (keine Veränderung)
        
        # C3 Rotationen um die z-Achse (Hauptachse)
        (0, 1, 4, 2, 3),  # 120° Rotation
        (0, 1, 3, 4, 2),  # 240° Rotation
        
        # C2 Rotationen um äquatoriale Achsen (senkrecht zur Hauptachse)
        (1, 0, 2, 4, 3),  # 180° Rotation um die x-Achse
        (1, 0, 4, 3, 2),  # 180° Rotation um eine der anderen äquatorialen Achsen
        (1, 0, 3, 2, 4),  # 180° Rotation um die dritte äquatoriale Achse
        
        # Horizontale Spiegelung (σh) - Spiegelung an der xy-Ebene
        (1, 0, 2, 3, 4),  # Vertauschung der axialen Positionen
        
        # Vertikale Spiegelungen (σv) - Spiegelungen an Ebenen, die die z-Achse enthalten
        (0, 1, 2, 4, 3),  # Spiegelung an der xz-Ebene
        (0, 1, 4, 3, 2),  # Spiegelung an einer weiteren vertikalen Ebene
        (0, 1, 3, 2, 4),  # Spiegelung an der dritten vertikalen Ebene
        
        # S3 improper Rotationen (C3 Rotation gefolgt von Spiegelung an der xy-Ebene)
        (1, 0, 3, 4, 2),  # 240° Rotation + Spiegelung
        (1, 0, 4, 2, 3),  # 120° Rotation + Spiegelung
    ]

    unique_structures = set()

    for perm in all_permutations:
        # Prüfe, ob die aktuelle Permutation durch Symmetrieoperationen abgedeckt ist
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

        # Wenn die Permutation einzigartig ist, füge sie hinzu
        if is_unique:
            unique_structures.add(perm)

    return unique_structures

def build_trigonal_bipyramidal_complex(db_path, central_atom, ligands):
    """
    Erstellt einen trigonal-bipyramidalen Komplex um das Zentralatom.
    """
    # Trigonal-bipyramidale Positionen relativ zum Ursprung
    positions = [
        np.array([0.0, 0.0, 2.0]),  # +z (axial)
        np.array([0.0, 0.0, -2.0]), # -z (axial)
        np.array([2.0, 0.0, 0.0]),  # +x (äquatorial)
        np.array([-1.0, np.sqrt(3), 0.0]),  # -x, +y (äquatorial)
        np.array([-1.0, -np.sqrt(3), 0.0])  # -x, -y (äquatorial)
    ]

    complex_atoms = [(central_atom, 0.0, 0.0, 0.0)]  # Zentralatom im Ursprung

    # Transformiere und füge jeden Liganden hinzu
    for ligand_name, position in zip(ligands, positions):
        ligand_atoms = fetch_atoms_for_ligand(db_path, ligand_name)
        transformed_ligand = transform_ligand(ligand_atoms, position)
        complex_atoms.extend(transformed_ligand)

    return complex_atoms

def save_single_trigonal_bipyramidal_arrangement(db_path, central_atom, ox, ligands, output_dir, spin=None, chrg=0, mult=1):
    """
    Speichert nur die direkt übergebenen 5 Liganden, ohne Permutationen.
    Dateiname ohne 1-1 usw.
    """
    # Nur ein Arrangement mit den übergebenen Liganden
    arrangement = ligands

    # Erstelle einen Unterordner für das Zentralatom
    central_atom_dir = os.path.join(output_dir, central_atom)
    os.makedirs(central_atom_dir, exist_ok=True)

    # Erstelle den Ordner für den spezifischen Komplex
    ligands_str = "_".join(arrangement)
    complex_dir = os.path.join(central_atom_dir, f"{central_atom}_{ox}_{ligands_str}")
    os.makedirs(complex_dir, exist_ok=True)

    # Spin-Suffix
    spin_suffix = ""
    if spin == "high":
        spin_suffix = "_highS"
    elif spin == "low":
        spin_suffix = "_lowS"

    # Dateiname ohne 1-1
    file_name = f"{central_atom}_{ox}_{ligands_str}{spin_suffix}.inp"
    inp_file_path = os.path.join(complex_dir, file_name)

    # Baue den Komplex
    complex_atoms = build_trigonal_bipyramidal_complex(db_path, central_atom, arrangement)

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
    default_output_dir = os.path.join(parent_dir, "Komplexe", "Trigonal Bipyramidal")  # Ordner für trigonal-bipyramidale Geometrie

    # Argumentparser einrichten
    parser = argparse.ArgumentParser(description="Generiert trigonal-bipyramidale Komplexe mit Liganden.")
    parser.add_argument("--db_path", default=default_db_path, help=f"Pfad zur SQLite-Datenbank (Standard: {default_db_path}).")
    parser.add_argument("--central_atom", required=True, help="Zentralatom des Komplexes (z. B. Fe).")
    parser.add_argument("--ox", type=int, required=True, help="Oxidationsstufe des Zentralatoms (Ganzzahl).")
    parser.add_argument("--ligands", required=True, help="Liste der Liganden, durch Kommas getrennt (z. B. Water,Ammonia,Chloride).")
    parser.add_argument("--output_dir", default=default_output_dir, help=f"Verzeichnis für die Ausgabe der .inp-Dateien (Standard: {default_output_dir}).")
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

    if len(selected_names) != 5:
        raise ValueError("Es müssen genau 5 Liganden angegeben werden.")

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Die Datenbank '{db_path}' wurde nicht gefunden.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    all_ligands = fetch_all_ligands(db_path)

    for name in selected_names:
        if name not in all_ligands:
            raise ValueError(f"Ligand '{name}' wurde nicht in der Datenbank gefunden.")

    save_single_trigonal_bipyramidal_arrangement(db_path, central_atom, ox, selected_names, output_dir, spin, chrg, mult)