import sqlite3
import os
import subprocess
import numpy as np
from itertools import combinations_with_replacement, permutations
import random

Skipped_Complexes = 0

def fetch_Octahedrale_metalle(db_path):
    """
    Fetches all metals with the geometry 'Oktaedrisch' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Query metals with geometry 'Oktaedrisch'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Oktaedrisch'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_Tetrahedrale_metalle(db_path):
    """
    Fetches all metals with the geometry 'Tetraedrisch' or 'Trigonal-pyramidal' with coordination number 4 from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Abfrage der Metalle mit Geometrie 'Tetraedrisch' 
    # ODER Metalle mit Geometrie 'Trigonal-pyramidal' und Koordinationszahl 4
    cursor.execute("""
        SELECT name, d_elektronen, oxidation 
        FROM metalle 
        WHERE geometrie = 'Tetraedrisch' 
        OR (geometrie = 'Trigonal-pyramidal' AND koordinationszahl = 4)
    """)

    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_quadratisch_planare_metalle(db_path):
    """
    Fetches all metals with the geometry 'Quadratisch-planar' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Quadratisch-planar'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Quadratisch-planar'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_trigonal_bipyramidale_metalle(db_path):
    """
    Fetches all metals with the geometry 'Trigonal-bipyramidal' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Trigonal-bipyramidal'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Trigonal-bipyramidal'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_trigonal_planare_metalle(db_path):
    """
    Fetches all metals with the geometry 'Trigonal-planar' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Trigonal-planar'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Trigonal-planar'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_gewinkelte_metalle(db_path):
    """
    Fetches all metals with the geometry 'Gewinkelt' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Gewinkelt'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Gewinkelt'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_lineare_metalle(db_path):
    """
    Fetches all metals with the geometry 'Linear' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Linear'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Linear'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_quadratisch_pyramidale_metalle(db_path):
    """
    Fetches all metals with the geometry 'Quadratisch-pyramidal' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Quadratisch-pyramidal'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Quadratisch-pyramidal'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_t_förmige_metalle(db_path):
    """
    Fetches all metals with the geometry 'T-förmig' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'T-förmig'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_trigonal_pyramidale_metalle(db_path):
    """
    Fetches all metals with the geometry 'Trigonal-pyramidal' and coordination number 3 from the database.
    (Those with coordination number 4 are treated as tetrahedral.)
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Nur Trigonal-pyramidale Metalle mit Koordinationszahl 3 abrufen
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Trigonal-pyramidal' AND koordinationszahl = 3")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_trigonal_prismatische_metalle(db_path):
    """
    Fetches all metals with the geometry 'Trigonal-prismatisch' from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Metalle mit Geometrie 'Trigonal-prismatisch'
    cursor.execute("SELECT name, d_elektronen, oxidation FROM metalle WHERE geometrie = 'Trigonal-prismatisch'")
    metalle = [{"name": row[0], "d_elektronen": row[1], "oxidation": row[2]} for row in cursor.fetchall()]
    conn.close()
    return metalle

def fetch_all_ligands(db_path):
    """
    Fetches all ligands from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Query all ligands
    cursor.execute("SELECT name, ladung FROM liganden")
    ligands = [{"name": row[0], "ladung": row[1]} for row in cursor.fetchall()]
    conn.close()
    return ligands

def fetch_atom_count_for_ligand(db_path, ligand_name):
    """
    Fetches the number of atoms of a ligand from the database.
    """
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    # Abfrage der Anzahl der Atome
    cursor.execute("SELECT xyz_daten FROM liganden WHERE name = ?", (ligand_name,))
    result = cursor.fetchone()
    if result is None:
        raise ValueError(f"Ligand '{ligand_name}' wurde nicht in der Datenbank gefunden.")
    # Anzahl der Zeilen in den XYZ-Daten entspricht der Anzahl der Atome
    atom_count = len(result[0].strip().split("\n"))
    conn.close()
    return atom_count

def calculate_multiplicity(d_electrons):
    """
    Returns the possible values for the multiplicity based on oxidation State / number of d-electrons.
    """
    multiplicity_map = {
        0: [1], 1: [2], 2: [3], 3: [4], 4: [3, 5], 5: [2, 6],
        6: [1, 5], 7: [2, 4], 8: [3], 9: [2], 10: [1]
    }
    return multiplicity_map.get(d_electrons, [])

def process_gewinkelte_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates bent complexes with 4 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 2))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            for mult in multiplicities:
                ligand_list = ",".join([l["name"] for l in ligand_set])
                # Erstelle den Namen des Unterordners
                complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
                complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
                # Überprüfen, ob der Ordner bereits existiert
                if os.path.exists(complex_dir_path):
                    Skipped_Complexes += 1
                    continue
                try:
                    # Aufruf von GewinkeltKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "GewinkeltKombi.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ],
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_lineare_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates linear complexes with 2 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 2))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            for mult in multiplicities:
                ligand_list = ",".join([l["name"] for l in ligand_set])
                # Erstelle den Namen des Unterordners
                complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
                complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
                # Überprüfen, ob der Ordner bereits existiert
                if os.path.exists(complex_dir_path):
                    Skipped_Complexes += 1
                    continue
                try:
                    # Aufruf von LinearKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "LinearKombi.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ],
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_linear_gewinkelt_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates linear or bent complexes with 2 ligands.
    """
    # Duplikate anhand der Oxidationsstufe entfernen:
    # Keep track of unique name+oxidation combinations
    unique_metal_ox = set()
    filtered_metalle = []
    for m in metalle:
        # Create a key from both name and oxidation state
        key = (m["name"], m["oxidation"])
        if key not in unique_metal_ox:
            unique_metal_ox.add(key)
            filtered_metalle.append(m)   
    metalle = filtered_metalle

    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 2))
        
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
                
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            
            for mult in multiplicities:
                ligand_list = ",".join([l["name"] for l in ligand_set])
                
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                
                # Erstelle den Namen des Unterordners
                complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
                complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
                
                # Überprüfen, ob der Ordner bereits existiert
                if os.path.exists(complex_dir_path):
                    Skipped_Complexes += 1
                    continue
                    
                try:
                    # Aufruf von LinearKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "LinearKombi.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_tetrahedral_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates tetrahedral complexes with 4 ligands.
    """
    global Skipped_Complexes
    
    # Duplikate anhand der Oxidationsstufe entfernen:
    # Keep track of unique name+oxidation combinations
    unique_metal_ox = set()
    filtered_metalle = []
    for m in metalle:
        # Create a key from both name and oxidation state
        key = (m["name"], m["oxidation"])
        if key not in unique_metal_ox:
            unique_metal_ox.add(key)
            filtered_metalle.append(m)   
    metalle = filtered_metalle

    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 4))
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
                
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            
            # Erstelle Ligandenliste und Ordnernamen außerhalb der Multiplizitätsschleife
            ligand_list = ",".join([l["name"] for l in ligand_set])
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert - VOR der Multiplizitätsschleife
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            # Jetzt durchlaufe die Multiplizitäten
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von TetraederKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TetraederKombi.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_Octahedrale_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir, metall_db_path):
    global Skipped_Complexes
    # Fetch metals with geometry 'Octahedral'
    metalle = fetch_Octahedrale_metalle(metall_db_path)
    if not metalle:
        print("No metals with geometry 'Octahedral' found.")
        return
    # Fetch all ligands
    ligands = fetch_all_ligands(ligand_db_path)
    if len(ligands) < 1:
        print("No ligands found in the database.")
        return
    # Output the names of all found ligands
    ligand_names = [ligand["name"] for ligand in ligands]
    # Counter for created complexes
    # For each metal, create all possible combinations of 6 ligands (with replacement) and select only max. 25 randomly
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Create all possible ligand combinations and shuffle randomly
        all_combos = list(combinations_with_replacement(ligands, 6))
        # Loop over the randomly shuffled combinations
        for ligand_set in all_combos:
            # Calculate the total number of atoms in the combination
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {ligand_names}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Calculate the total charge (sum of ligand charges + oxidation state of the metal)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Calculate the possible multiplicities
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            # Create the name of the subfolder
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Check if the folder already exists
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                ligand_list = ",".join([l["name"] for l in ligand_set])
                # Determine the spin state based on the number of d-electrons
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Lower multiplicity -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Higher multiplicity -> high spin
                        spin_suffix = "_highS"
                try:
                    # Call OCKombi.py as a separate process
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "Octahedral.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_quadratisch_planare_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates square planar complexes with 4 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 4))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von SPKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "SquarePlanar.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir, 
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_trigonal_bipyramidale_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates trigonal-bipyramidal complexes with 5 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 5))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von TrigonalBipyramidalKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TrigonalBipyramidal.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_trigonal_planare_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates trigonal-planar complexes with 3 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 3))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von TrigonalPlanarKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TrigonalPlanar.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_quadratisch_pyramidale_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates square pyramidal complexes with 5 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 5))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von QuadratischPyramidalKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "SquarePyramidal.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir, 
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_t_förmige_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates T-shaped complexes with 3 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 3))
        random.shuffle(all_combos)
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von TFörmigKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TShaped.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_trigonal_pyramidale_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates trigonal-pyramidal complexes with 3 ligands.
    (Those with 4 ligands are treated as tetrahedral)
    """
    global Skipped_Complexes
    for metall in metalle:
        # Die spezifische Unterordner-Erstellung entfernen
        # Direkt den Hauptordner verwenden
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        os.makedirs(central_atom_dir, exist_ok=True)

        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 3))  # Feste Koordinationszahl 3
        random.shuffle(all_combos)
        
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
                
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            
            # Erstelle Ligandenliste und Ordnernamen außerhalb der Multiplizitätsschleife
            ligand_list = ",".join([l["name"] for l in ligand_set])
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            
            # Überprüfen, ob der Ordner bereits existiert - VOR der Multiplizitätsschleife
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
                
            # Erstelle den Ordner jetzt, da wir wissen, dass er noch nicht existiert
            os.makedirs(complex_dir_path, exist_ok=True)
            
            # Jetzt durchlaufe die Multiplizitäten
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                
                try:
                    # Aufruf von TrigonalPyramidalKombi.py als separater Prozess - output_dir direkt verwenden
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TrigonalPyramidal.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,  # Direkt das Hauptverzeichnis statt des spezifischen Unterordners
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def process_trigonal_prismatische_geometry(metalle, ligands, ligand_db_path, output_base_dir, base_dir):
    """
    Creates trigonal-prismatic complexes with 6 ligands.
    """
    global Skipped_Complexes
    for metall in metalle:
        central_atom_dir = os.path.join(output_base_dir, metall['name'])
        if not os.path.exists(central_atom_dir):
            os.makedirs(central_atom_dir)
        # Alle möglichen LigandenCombinationsen erstellen und zufällig mischen
        all_combos = list(combinations_with_replacement(ligands, 6))
        for ligand_set in all_combos:
            # Berechne die Gesamtanzahl der Atome in der Combinations
            total_atoms = sum(fetch_atom_count_for_ligand(ligand_db_path, ligand["name"]) for ligand in ligand_set)
            if total_atoms > 40:
                print(f"Skipping combination {[ligand['name'] for ligand in ligand_set]}, because the total number of atoms is {total_atoms} (more than 40).")
                continue
            # Berechne die Gesamtladung (Summe der Ligandenladungen + Oxidationsstufe des Metalls)
            ligand_charges = [ligand["ladung"] for ligand in ligand_set]
            total_charge = sum(ligand_charges) + metall["oxidation"]
            # Berechne die möglichen Multiplizitäten
            multiplicities = calculate_multiplicity(metall["d_elektronen"])
            ligand_list = ",".join([l["name"] for l in ligand_set])
            # Erstelle den Namen des Unterordners
            complex_dir_name = f"{metall['name']}_{metall['oxidation']}_{'_'.join([ligand['name'] for ligand in ligand_set])}"
            complex_dir_path = os.path.join(central_atom_dir, complex_dir_name)
            # Überprüfen, ob der Ordner bereits existiert
            if os.path.exists(complex_dir_path):
                Skipped_Complexes += 1
                continue
            for mult in multiplicities:
                # Bestimme den Spin-Zustand basierend auf der Anzahl der d-Elektronen
                spin_suffix = ""
                if 4 <= metall["d_elektronen"] <= 7:
                    if mult == min(multiplicities):  # Niedrigere Multiplizität -> low spin
                        spin_suffix = "_lowS"
                    elif mult == max(multiplicities):  # Höhere Multiplizität -> high spin
                        spin_suffix = "_highS"
                try:
                    # Aufruf von TrigonalPrismatischKombi.py als separater Prozess
                    subprocess.run(
                        [
                            "python",
                            os.path.join(base_dir, "Combinations", "TrigonalPrismatic.py"),
                            "--db_path", ligand_db_path,
                            "--central_atom", metall["name"],
                            "--ligands", ligand_list,
                            "--output_dir", output_base_dir,
                            "--chrg", str(total_charge),
                            "--mult", str(mult),
                            "--ox", str(metall["oxidation"]),
                        ] + (["--spin", "low"] if spin_suffix == "_lowS" else ["--spin", "high"] if spin_suffix == "_highS" else []),
                        check=True
                    )
                except subprocess.CalledProcessError as e:
                    print(f"Error creating a complex for metal '{metall['name']}' with ligands '{ligand_list}': {e}")

def count_configurations_and_files(base_dir):
    """
    Counts the number of created configurations (folders) and .inp files in the directory.
    """
    total_folders = 0
    total_inp_files = 0
    
    # Alle Geometrie-Verzeichnisse durchlaufen
    geometry_dirs = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    for geometry in geometry_dirs:
        geometry_path = os.path.join(base_dir, geometry)
        
        # Alle Metall-Verzeichnisse durchlaufen
        metal_dirs = [d for d in os.listdir(geometry_path) if os.path.isdir(os.path.join(geometry_path, d))]
        
        for metal in metal_dirs:
            metal_path = os.path.join(geometry_path, metal)
            
            # Alle Komplex-Verzeichnisse zählen
            complex_dirs = [d for d in os.listdir(metal_path) if os.path.isdir(os.path.join(metal_path, d))]
            total_folders += len(complex_dirs)
            
            # Alle .inp-Dateien in allen Komplex-Verzeichnissen zählen
            for complex_dir in complex_dirs:
                complex_path = os.path.join(metal_path, complex_dir)
                inp_files = [f for f in os.listdir(complex_path) if f.endswith('.inp')]
                total_inp_files += len(inp_files)
    
    return total_folders, total_inp_files

def main():
    global Skipped_Complexes
    # Datenbankpfade
    base_dir = os.path.dirname(os.path.abspath(__file__))
    metall_db_path = os.path.join(base_dir, "Database", "metals.db")
    ligand_db_path = os.path.join(base_dir, "Database", "ligands.db")
    output_base_dir_Octahedral = os.path.join(base_dir, "Complexes ", "Octahedral")
    output_base_dir_Tetrahedral = os.path.join(base_dir, "Complexes", "Tetrahedral")
    output_base_dir_quadratisch_planar = os.path.join(base_dir, "Complexes", "Square Planar")
    output_base_dir_trigonal_bipyramidal = os.path.join(base_dir, "Complexes", "Trigonal Bipyramidal")
    output_base_dir_trigonal_planar = os.path.join(base_dir, "Complexes", "Trigonal Planar")
    output_base_dir_quadratisch_pyramidal = os.path.join(base_dir, "Complexes", "Square Pyramidal")
    output_base_dir_t_förmig = os.path.join(base_dir, "Complexes", "T Shaped")
    output_base_dir_trigonal_pyramidal = os.path.join(base_dir, "Complexes", "Trigonal Pyramidal")
    output_base_dir_trigonal_prismatisch = os.path.join(base_dir, "Complexes", "Trigonal Prismatic")

    # Trigonal Prismatic Geometrie
    metalle_trigonal_prismatisch = fetch_trigonal_prismatische_metalle(metall_db_path)
    if metalle_trigonal_prismatisch:
        ligands = fetch_all_ligands(ligand_db_path)
        process_trigonal_prismatische_geometry(metalle_trigonal_prismatisch, ligands, ligand_db_path, output_base_dir_trigonal_prismatisch, base_dir)
    print("Trigonal Prismatic done")
    # Trigonal Pyramidal Geometrie
    metalle_trigonal_pyramidal = fetch_trigonal_pyramidale_metalle(metall_db_path)
    if metalle_trigonal_pyramidal:
        ligands = fetch_all_ligands(ligand_db_path)
        process_trigonal_pyramidale_geometry(metalle_trigonal_pyramidal, ligands, ligand_db_path, output_base_dir_trigonal_pyramidal, base_dir)
    print("Trigonal Pyramidal done")
    # T-förmige Geometrie
    metalle_t_förmig = fetch_t_förmige_metalle(metall_db_path)
    if metalle_t_förmig:
        ligands = fetch_all_ligands(ligand_db_path)
        process_t_förmige_geometry(metalle_t_förmig, ligands, ligand_db_path, output_base_dir_t_förmig, base_dir)
    print("T-förmig done")
    # Square Pyramidal Geometrie
    metalle_quadratisch_pyramidal = fetch_quadratisch_pyramidale_metalle(metall_db_path)
    if metalle_quadratisch_pyramidal:
        ligands = fetch_all_ligands(ligand_db_path)
        process_quadratisch_pyramidale_geometry(metalle_quadratisch_pyramidal, ligands, ligand_db_path, output_base_dir_quadratisch_pyramidal, base_dir)
    print("Square Pyramidal done")
    # Trigonal planare Geometrie
    metalle_trigonal_planar = fetch_trigonal_planare_metalle(metall_db_path)
    if metalle_trigonal_planar:
        ligands = fetch_all_ligands(ligand_db_path)
        process_trigonal_planare_geometry(metalle_trigonal_planar, ligands, ligand_db_path, output_base_dir_trigonal_planar, base_dir)
    print("Trigonal Planar done")
    # Trigonal Bipyramidale Geometrie
    metalle_trigonal_bipyramidal = fetch_trigonal_bipyramidale_metalle(metall_db_path)
    if metalle_trigonal_bipyramidal:
        ligands = fetch_all_ligands(ligand_db_path)
        process_trigonal_bipyramidale_geometry(metalle_trigonal_bipyramidal, ligands, ligand_db_path, output_base_dir_trigonal_bipyramidal, base_dir)
    print("Trigonal Bipyramidal done")
    # Square Planare Geometrie
    metalle_quadratisch_planar = fetch_quadratisch_planare_metalle(metall_db_path)
    if metalle_quadratisch_planar:
        ligands = fetch_all_ligands(ligand_db_path)
        process_quadratisch_planare_geometry(metalle_quadratisch_planar, ligands, ligand_db_path, output_base_dir_quadratisch_planar, base_dir)
    print("Square Planar done")
    # Octahedrale Geometrie
    metalle_Octahedral = fetch_Octahedrale_metalle(metall_db_path)
    if metalle_Octahedral:
        ligands = fetch_all_ligands(ligand_db_path)
        process_Octahedrale_geometry(metalle_Octahedral, ligands, ligand_db_path, output_base_dir_Octahedral, base_dir, metall_db_path)
    print("Octahedral done")
    # Tetrahedrale Geometrie
    metalle_Tetrahedral = fetch_Tetrahedrale_metalle(metall_db_path)
    if metalle_Tetrahedral:
        ligands = fetch_all_ligands(ligand_db_path)
        process_tetrahedral_geometry(metalle_Tetrahedral, ligands, ligand_db_path, output_base_dir_Tetrahedral, base_dir)
    print("Tetrahedral done")

    # Print the total count of skipped complexes
    print(f"Skipped complexes: {Skipped_Complexes}")
    
    # Zähle erstellte Konfigurationen und Dateien
    Complexes_dir = os.path.join(base_dir, "Complexes")
    folder_count, file_count = count_configurations_and_files(Complexes_dir)
    print(f"Created configurations (folders): {folder_count}")
    print(f"Created .inp files: {file_count}")

if __name__ == "__main__":
    main()







