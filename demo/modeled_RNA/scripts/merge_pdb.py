# merge_pdb.py
import subprocess
from pymol import cmd
import sys
import os
import glob
import re

def check_input_directory(input_dir):
    if not os.path.isdir(input_dir):
        print(f"Directory {input_dir} does not exist.")
        sys.exit(1)


def check_output_filename(output_file):
    # Check output_file extension
    if not (output_file.lower().endswith('.pdb') or output_file.lower().endswith('.pse')):
        print(f"Error: output_file '{output_file}' must end with .pdb or .pse for PyMOL.")
        sys.exit(1)
        

def natsorted(filenames : list) -> list:
    """ Sorts filenames by digits in the filenames.
    Supports multiple numbers in a filename.
    Natsorted example: [
        'pocket1_atm.pdb',
        'pocket2_atm.pdb',
        '2pocket11_atm.pdb',
        'pocket10_atm.pdb',
        'pocket11_atm.pdb',
        'pocket11_atm5.pdb',
        'pocket20_atm.pdb'
    ]
    Args:
        filenames (list): List of strings to be sorted by digits.
    Returns:
        list: Sorted list, ordered by digits.
    """
    # return sorted(filenames, key=lambda x: int(re.search(r'\d+', x).group()))
    
    def extract_sort_key(filename):
        # Extract all numbers in the filename as integers, or return the filename as-is for non-numeric parts
        return [int(num) if num.isdigit() else num.lower() 
                for num in re.split(r'(\d+)', filename)]
    
    return sorted(filenames, key=extract_sort_key)


def is_valid_pdb(pdb_path):
    validation_code = (
        "from pymol import cmd\n"
        "cmd.load(r'{}', 'temp')\n"
        "cmd.quit()"
    ).format(pdb_path.replace("'", "\\'"))
    try:
        proc = subprocess.run([
            sys.executable, '-c', validation_code
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return proc.returncode == 0
    except Exception as e:
        print(f"Error running PyMOL subprocess for {pdb_path}: {e}")
        return False


def merge_pdb_files(sample_name, input_dir, output_file, intra_fit_file):
    # Rename files in input_dir without an extension to have .pdb
    for filename in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir, filename)) and not os.path.splitext(filename)[1]:
            os.rename(os.path.join(input_dir, filename), os.path.join(input_dir, filename + '.pdb'))

    pdb_files = natsorted(glob.glob(os.path.join(input_dir, '*.pdb')))
    if not pdb_files:
        print("No PDB files found in the directory.")
        sys.exit(1)

    valid_pdb_files = []
    for pdb_file in pdb_files:
        print(f"Testing PDB file: {pdb_file}")
        if is_valid_pdb(pdb_file):
            valid_pdb_files.append(pdb_file)
        else:
            print(f"PyMOL failed to open {pdb_file}. Skipping this file.")

    if not valid_pdb_files:
        print("No valid PDB files found after validation.")
        sys.exit(1)

    # Find state index by comparing basenames
    state_index = None
    intra_basename = os.path.basename(intra_fit_file)
    for i, pdb_file in enumerate(valid_pdb_files):
        model = os.path.splitext(os.path.basename(pdb_file))[0]
        print(model)
        cmd.load(pdb_file, f'{model}')
        cmd.create(f'{sample_name}_multistate', f'{model}', 1, i+1)
        if os.path.basename(pdb_file) == intra_basename:
            state_index = i + 1

    if state_index is None:
        print(f"Error: intra_fit_file '{intra_fit_file}' not found among valid PDB files.")
        cmd.quit()
        sys.exit(1)

    cmd.alter('all', "chain='A'")
    cmd.alter('all', "segi='A'")
    cmd.intra_fit(f'{sample_name}_multistate', state_index)
    cmd.save(output_file, f'{sample_name}_multistate', state=0)
    cmd.quit()
    print(f'\nFinished merging states for {sample_name}!\n')


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print("Usage: python merge_pdb.py <sample_name> <input_dir> <output_file> <intra_fit_file>\n\
               sample_name(str): Name of the sample to be used in PyMOL.\n\
               input_dir(path): Directory containing PDB files to merge.\n\
               output_file(path): Name of the output file (should end with .pdb or .pse).\n\
               intra_fit_file(path): PDB file to use as the target state for intra_fit alignment.")
        sys.exit(1)

    sample_name = sys.argv[1]
    input_dir = sys.argv[2]
    output_file = sys.argv[3]
    intra_fit_file = sys.argv[4]

    check_input_directory(input_dir)
    check_output_filename(output_file)
    merge_pdb_files(sample_name, input_dir, output_file, intra_fit_file)
    print(f"Merged PDB file saved as {output_file}")