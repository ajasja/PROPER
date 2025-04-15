"""
Utility functions
"""
import Bio
from Bio.PDB import PDBParser, PDBIO, Select
from Bio import SeqIO
import sys
from io import StringIO

resname_3to1_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 'TER':'*',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA':'X', 'MSE': 'M'}

def resname_3to1(resname3: list) -> str:
    """Takes an array of 3 letter names and returns one letter names"""
    result = [resname_3to1_dict[resname] if resname in resname_3to1_dict else 'X' for resname in resname3]
    return result


def print_pymol_selection(df_with_resi_index0):
    """Returns a pymol selection of resids (counted from 1) from the resi_index0 column"""
    resids ='+'.join([str(i+1) for i in df_with_resi_index0.resi_index0])
    print(f'select resi {resids}')

def permute_seq_with_linker(seq: str, position0: int, linker: str) -> str:
    """
    Return a new string composed of the first `position` characters of `seq`,
    followed by `linker`, replacing the rest.
    """
    if position0 < 0 or position0 > len(seq):
        raise ValueError("Position must be within the length of the original sequence.")

    return seq[position0:] + linker + seq[:position0]


# Custom Select class to filter only protein atoms of the first chain
class ProteinFirstChainSelect(Select):
    def __init__(self, target_chain_id):
        self.target_chain_id = target_chain_id

    def accept_chain(self, chain):
        return chain.id == self.target_chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " "  # " " means it's a standard residue (not water/hetero)

    def accept_atom(self, atom):
        return True  # Accept all atoms in the residue

def renumber_residues_structure(structure, start_res_id=1):
    """
    Renumber residues in a Biopython structure object, starting from start_res_id (default 1).
    Each chain is renumbered independently.

    Parameters:
    - structure: Bio.PDB.Structure.Structure
    - start_res_id: int (default 1)

    Returns:
    - structure: renumbered Bio.PDB.Structure.Structure object
    """
    for model in structure:
        for chain in model:
            res_id = start_res_id
            for residue in chain:
                old_id = residue.id
                new_id = (old_id[0], res_id, old_id[2])
                residue.id = new_id
                res_id += 1
    return structure

def renumber_residues(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)
    renumber_residues_structure(structure)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def save_first_chain_protein_atoms(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    model = structure[0]  # First model
    first_chain = next(model.get_chains())  # First chain

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ProteinFirstChainSelect(first_chain.id))

def clean_pdb(input_pdb, output_pdb=None):
    """Takes the protein atoms of the first chain"""
    if output_pdb is None:
        output_pdb = input_pdb
    save_first_chain_protein_atoms(input_pdb, output_pdb)    
    ##TODO: check for breaks!
    renumber_residues(output_pdb, output_pdb)

def get_middle_loop_row(group):
    """
    Return the middle row of a DataFrame group sorted by 'resi_index0'.

    The DataFrame is first sorted in ascending order by the 'resi_index0' column,
    and then the middle row is selected. For groups with an even number of rows,
    the function returns the lower of the two middle rows.

    Parameters:
        group (pd.DataFrame): A pandas DataFrame containing a 'resi_index0' column.

    Returns:
        pd.Series: The middle row of the sorted DataFrame as a pandas Series.
    """
    sorted_group = group.sort_values(by='resi_index0').reset_index(drop=True)
    mid_index = len(sorted_group) // 2
    # For even-sized groups, this returns the lower middle
    return sorted_group.iloc[mid_index - 1 if len(sorted_group) % 2 == 0 else mid_index]

def is_one_line(string):
    """
    Check if a string contains only one line.
    
    Args:
        string (str): The string to check
        
    Returns:
        bool: True if the string contains only one line, False otherwise
    """
    # Remove trailing newlines to avoid counting empty lines at the end
    string = string.rstrip('\n\r')
    
    # Check if there are any newline characters in the string
    return '\n' not in string and '\r' not in string


def load_fasta_dictionary(fasta_str):
    fasta_io = StringIO(fasta_str)
    fasta_dict={}
     # Parse the FASTA data using Biopython
    for record in Bio.SeqIO.parse(fasta_io, "fasta"):
        # Use the sequence ID as the key and the sequence as the value
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict
    
    return fasta_dict

def load_fasta_or_single_line(fasta_str, default_name='linker'):
    if is_one_line(fasta_str):
        return {default_name: fasta_str}
    else:  
        return load_fasta_dictionary(fasta_str)
        

def generate_permuted_sequences(construct_name, seq, linkers, positions1, out_file):
    """
    Generates permuted protein sequences with linkers and writes them to a file.

    Args:
        construct_name (str): The name of the construct.
        seq (str): The original protein sequence.
        linkers (dict): A dictionary of linker names and sequences.
        positions1 (list): A list of positions for permutation.
        out_file (str): The path to the output file.

    Returns:
        None
    """
    with open(out_file, 'w') as f:
        for pos1 in positions1:
            for linker_name, linker in linkers.items():
                name = f"{construct_name}{pos1:03d}-{linker_name}"
                per_seq = permute_seq_with_linker(seq, position0=pos1 - 1, linker=linker)
                f.write(">" + name + "\n")
                f.write(per_seq + "\n")


def file_to_str(file_path):
    """Loads the contents of a file and returns it as a string.

    Args:
        file_path: The path to the file.

    Returns:
        The contents of the file as a string.
    """
    with open(file_path, 'r') as f:
        file_content = f.read()
    return file_content

