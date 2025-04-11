"""
Utility functions
"""

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
    resids ='+'.join([str(i+1) for i in q.resi_index0])
    print(f'select resi {resids}')

def permute_seq_with_linker(seq: str, position0: int, linker: str) -> str:
    """
    Return a new string composed of the first `position` characters of `seq`,
    followed by `linker`, replacing the rest.
    """
    if position0 < 0 or position0 > len(seq):
        raise ValueError("Position must be within the length of the original sequence.")

    return seq[position0:] + linker + seq[:position0]


from Bio.PDB import PDBParser, PDBIO, Select
import sys

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

def save_first_chain_protein_atoms(input_pdb, output_pdb):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", input_pdb)

    model = structure[0]  # First model
    first_chain = next(model.get_chains())  # First chain

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb, select=ProteinFirstChainSelect(first_chain.id))