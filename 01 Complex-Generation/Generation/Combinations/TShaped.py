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
    Generiert eindeutige Permutationen der Liganden unter Berücksichtigung der C2v Symmetrie
    und mehrfach vorkommender Liganden.
    """
    # Alle Permutationen der Liganden
    all_permutations = set(permutations(ligands))
    # Korrekte Symmetrieoperationen für C2v ohne Duplikate
    symmetries = [
        # Identität
        (0, 1, 2), 
        # C2
        (0, 2, 1),
# sigma_v
        (0, 2, 1), #xz 
        (0, 1, 2), #yz
# sigma_v
        (0, 2, 1), #xz 
        (0, 1, 2), #yz
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

def build_t_shaped_complex(db_path, central_atom, ligands):
    """
    Erstellt einen T-förmigen Komplex um das Zentralatom.
    """
    # T-förmige Positionen relativ zum Ursprung
    positions = [
        np.array([0.0, 1.5, 0.0]),   # +y
        np.array([-1.5, 0.0, 0.0]),  # -x
        np.array([1.5, 0.0, 0.0])    # +x
    ]

    complex_atoms = [(central_atom, 0.0, 0.0, 0.0)]  # Zentralatom im Ursprung

    # Transformiere und füge jeden Liganden hinzu
    for ligand_name, position in zip(ligands, positions):
        ligand_atoms = fetch_atoms_for_ligand(db_path, ligand_name)
        transformed_ligand = transform_ligand(ligand_atoms, position)
        complex_atoms.extend(transformed_ligand)

    return complex_atoms

def save_all_t_shaped_arrangements(db_path, central_atom, ox, ligands, output_dir, spin=None, chrg=0, mult=1):
    """
    Speichert nur die direkt übergebenen 3 Liganden als einen T-förmigen Komplex,
    ohne Permutationen zu bilden oder i-total Suffixe zu verwenden.
    """
    # Nur ein Arrangement (direkt die übergebenen Liganden)
    arrangement = ligands

    # Erstelle einen Unterordner für das Zentralatom
    central_atom_dir = os.path.join(output_dir, central_atom)
    os.makedirs(central_atom_dir, exist_ok=True)

    # Erstelle den Ordner für den spezifischen Komplex
    ligands_str = "_".join(arrangement)
    complex_dir = os.path.join(central_atom_dir, f"{central_atom}_{ox}_{ligands_str}")
    os.makedirs(complex_dir, exist_ok=True)

    # Bestimme den Spin-Suffix
    spin_suffix = ""
    if spin == "high":
        spin_suffix = "_highS"
    elif spin == "low":
        spin_suffix = "_lowS"

    # Dateiname vom Typ: Zentralatom_Ox_Ligand1_Ligand2_Ligand3[Spin].inp
    file_name = f"{central_atom}_{ox}_{ligands_str}{spin_suffix}.inp"
    inp_file_path = os.path.join(complex_dir, file_name)

    # Baue den T-förmigen Komplex
    complex_atoms = build_t_shaped_complex(db_path, central_atom, arrangement)

    # Schreibe die .inp-Datei
    with open(inp_file_path, 'w') as inp_file:
        inp_file.write("! B3LYP Def2-SVP Opt D3BJ\n%PAL NPROCS 16 END\n")
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
    default_output_dir_t_shaped = os.path.join(parent_dir, "Komplexe", "TFörmig")  # Ordner für T-förmige Komplexe

    # Argumentparser einrichten
    parser = argparse.ArgumentParser(description="Generiert Komplexe mit Liganden.")
    parser.add_argument("--db_path", default=default_db_path, help=f"Pfad zur SQLite-Datenbank (Standard: {default_db_path}).")
    parser.add_argument("--central_atom", required=True, help="Zentralatom des Komplexes (z. B. Fe).")
    parser.add_argument("--ox", type=int, required=True, help="Oxidationsstufe des Zentralatoms (Ganzzahl).")
    parser.add_argument("--ligands", required=True, help="Liste der Liganden, durch Kommas getrennt (z. B. Water,Ammonia,Chloride).")
    parser.add_argument("--output_dir_t_shaped", default=default_output_dir_t_shaped, help=f"Verzeichnis für die Ausgabe der T-förmigen Komplexe (Standard: {default_output_dir_t_shaped}).")
    parser.add_argument("--chrg", type=int, required=True, help="Ladung des Komplexes (Ganzzahl).")
    parser.add_argument("--mult", type=int, required=True, help="Multiplizität des Komplexes (Ganzzahl).")
    parser.add_argument("--spin", choices=["high", "low"], help="Spin-Zustand des Komplexes (high oder low).")

    args = parser.parse_args()

    # Eingaben aus den Argumenten
    db_path = args.db_path
    central_atom = args.central_atom
    ox = args.ox
    selected_names = [name.strip() for name in args.ligands.split(",")]
    output_dir_t_shaped = args.output_dir_t_shaped
    chrg = args.chrg
    mult = args.mult
    spin = args.spin

    # Überprüfen, ob die Datenbank existiert
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Die Datenbank '{db_path}' wurde nicht gefunden.")

    # Überprüfen, ob die Ausgabe-Verzeichnisse existieren, und erstellen, falls nicht
    if not os.path.exists(output_dir_t_shaped):
        os.makedirs(output_dir_t_shaped)

    # Abrufen aller Liganden aus der Datenbank
    all_ligands = fetch_all_ligands(db_path)

    # Überprüfen, ob alle angegebenen Liganden in der Datenbank existieren
    for name in selected_names:
        if name not in all_ligands:
            raise ValueError(f"Ligand '{name}' wurde nicht in der Datenbank gefunden.")

    # Generiere und speichere alle möglichen Anordnungen für T-förmige Komplexe
    if len(selected_names) == 3:
        save_all_t_shaped_arrangements(db_path, central_atom, ox, selected_names, output_dir_t_shaped, spin, chrg, mult)
    else:
        raise ValueError("Es müssen entweder 3 oder 6 Liganden angegeben werden.")