import os
import re
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('goat_conversion.log')
    ]
)
logger = logging.getLogger('goat_converter')

def read_xyz_file(xyz_path):
    """Read an XYZ file and return atom count, title and coordinates"""
    with open(xyz_path, 'r') as f:
        lines = f.readlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {xyz_path}")
    atom_count = int(lines[0].strip())
    title = lines[1].strip()
    coordinates = [line.strip() for line in lines[2:] if line.strip()]
    return atom_count, title, coordinates

def create_goat_input(xyz_path, original_inp_path=None, charge=None, multiplicity=None):
    """Create GOAT input content from XYZ file"""
    atom_count, title, coordinates = read_xyz_file(xyz_path)
    # Try to determine charge and multiplicity from original input file if available
    if original_inp_path and os.path.exists(original_inp_path):
        try:
            with open(original_inp_path, 'r') as f:
                content = f.read()
                xyz_match = re.search(r'\*\s*xyz\s+(-?\d+)\s+(\d+)', content)
                if xyz_match:
                    charge = int(xyz_match.group(1))
                    multiplicity = int(xyz_match.group(2))
        except Exception as e:
            logger.warning(f"Could not parse original input file: {original_inp_path} ({e})")
    if charge is None or multiplicity is None:
        return None
    input_content = [
        "!GOAT XTB TightSCF",
        "%PAL NPROCS 8 END",
        "%maxcore 3000",
        "%GOAT",
        "NWORKERS 8",
        "AUTOWALL TRUE",
        "END",
        f"* xyz {charge} {multiplicity}",
    ]
    input_content.extend(coordinates)
    input_content.append("*")
    return '\n'.join(input_content)

def find_xyz_file(xyz_dirs, base_filename):
    """Durchsuche alle xyz_dirs rekursiv nach einer .xyz mit base_filename"""
    for xyz_dir in xyz_dirs:
        for root, dirs, files in os.walk(xyz_dir):
            for file in files:
                if file == base_filename:
                    return os.path.join(root, file)
    return None

def process_directory(base_dir, goat_dir, test_mode=False, limit=5):
    goat_output_dir = goat_dir
    os.makedirs(goat_output_dir, exist_ok=True)

    # 1. Alle .inp-Dateien rekursiv finden
    inp_files = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".inp"):
                inp_files.append(os.path.join(root, file))

    # 2. Alle Verzeichnisse merken, um darin nach .xyz zu suchen
    all_dirs = set()
    for root, dirs, files in os.walk(base_dir):
        all_dirs.add(root)

    logger.info(f"Found {len(inp_files)} .inp files to check")
    logger.info(f"GOAT input files will be saved to: {goat_output_dir}")

    successful = 0
    failed = 0

    for inp_path in inp_files:
        try:
            if test_mode and successful >= limit:
                logger.info(f"Test mode: Reached limit of {limit} files. Stopping.")
                break

            base_name = os.path.splitext(os.path.basename(inp_path))[0]
            xyz_filename = base_name + ".xyz"
            xyz_path = find_xyz_file(all_dirs, xyz_filename)
            if not xyz_path:
                logger.warning(f"No matching .xyz file found for {inp_path}")
                failed += 1
                continue

            goat_content = create_goat_input(xyz_path, inp_path)
            if goat_content is None:
                logger.warning(f"Could not create GOAT input for {inp_path} (charge/multiplicity missing)")
                failed += 1
                continue

            if not test_mode:
                new_filename = f"{base_name}.inp"
                new_path = os.path.join(goat_output_dir, new_filename)
                with open(new_path, 'w') as f:
                    f.write(goat_content)
                logger.info(f"Created GOAT input file: {new_path}")

            successful += 1

        except Exception as e:
            logger.error(f"Error processing {inp_path}: {str(e)}")
            failed += 1

    logger.info(f"Anzahl .inp-Dateien gefunden: {len(inp_files)}")
    logger.info(f"Anzahl GOAT-.inp-Dateien {'(theoretisch) ' if test_mode else ''}erstellt: {successful}")

    if test_mode:
        logger.info(f"TEST MODE: Would convert {successful} files successfully (limited to {limit}), {failed} failed")
    else:
        logger.info(f"Conversion complete: {successful} files converted successfully, {failed} failed")

    return successful, failed, len(inp_files)

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Convert XTB input/xyz pairs to GOAT inputs")
    parser.add_argument("--dir", "-d", default="XTB_Data", help="Base directory to search for XTB calculations")
    parser.add_argument("--goat-dir", "-g", default="D:\GOAT_Inp2",
                        help="Output directory for GOAT input files (default: D:\\GOAT_INP)")
    parser.add_argument("--test-mode", "-t", action="store_true", help="Enable test mode (process limited number of files)")
    parser.add_argument("--limit", "-l", type=int, default=5, help="Limit number of files to process in test mode (default: 5)")

    args = parser.parse_args()

    base_dir = os.path.abspath(args.dir)
    if not os.path.isdir(base_dir):
        logger.error(f"Directory not found: {base_dir}")
        return

    logger.info(f"Processing directory: {base_dir}")
    successful, failed, total_inp = process_directory(base_dir, args.goat_dir, args.test_mode, args.limit)
    if successful > 0:
        logger.info(f"Successfully converted {successful} files. GOAT input files are in the directory: {args.goat_dir}")
    else:
        logger.warning("No successful calculations were found or converted.")

if __name__ == "__main__":
    main()