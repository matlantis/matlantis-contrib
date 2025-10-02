###
#This module was modified from pfcc-extras library.
#Copyright Matlantis Corp. as contributors to Matlantis contrib project
###
import networkx as nx
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.rdchem import BondType
from rdkit.Chem.rdmolfiles import CanonicalRankAtoms
from ase.io import write
from ase.neighborlist import NeighborList, natural_cutoffs


def get_connectivity_matrix(atoms, mult=0.9, self_interaction=False, bothways=True):
    """Get connectivity matrix of an atoms object.

    Args:
        atoms (ase.Atoms): An atoms object
        mult (float, optional): A multiplier for all cutoffs used in ase.neighborlist.natural_cutoffs. Defaults to 0.9.
        self_interaction (bool, optional): Return the atom itself as its own neighbor if set to true. Default to False.
        bothways (bool, optional): Return all neighbors. Default is to return only “half” of the neighbors. Defaults to True.

    Returns:
        matrix (a scipy csr matrix): connectivity matrix
    """
    cutoff = natural_cutoffs(atoms, mult=mult)
    neighborList = NeighborList(
        cutoff, self_interaction=self_interaction, bothways=bothways
    )
    neighborList.update(atoms)
    matrix = neighborList.get_connectivity_matrix()
    return matrix


def atoms_to_nxGraph(atoms, mult=0.9, props=["symbol"], matrix=None, add_bond_info=False):
    """Convert ase.Atoms to a networkx Graph object.
    The graph nodes and edges represent the elements and connectivity of the atoms, respectively.

    Args:
        atoms (ase.Atoms): an atoms object.
        mult (float, optional): A multiplier for all cutoffs used in ase.neighborlist.natural_cutoffs. Defaults to 0.9.
        props (list of strings, optional): properties of the atoms which will be stored in the graph nodes.
            Defaults to None.
        matrix (scipy csr matrix, optional): connectivity matrix of the atoms. Defaults to None.

    Returns:
        graph: a networkx Graph object.
    """
    if props is None:
        props = []
    if matrix is None:
        matrix = get_connectivity_matrix(atoms, mult=mult)

    g = nx.from_scipy_sparse_array(matrix)

    for node in (g.nodes):
        g.nodes[node]["idx"] = node
        for attrib in props:
            if hasattr(atoms[node], attrib):
                g.nodes[node][attrib] = getattr(atoms[node], attrib)
            else:
                # print(f'{atoms} does not have attribute "{attrib}".')
                pass
    
    if add_bond_info:
        write("temp.xyz", atoms)
        mol = Chem.MolFromXYZFile("temp.xyz")
        rdDetermineBonds.DetermineBonds(mol)
        rank_atoms = CanonicalRankAtoms(mol, breakTies=False)
        bond_name_dict ={v:k for k,v in BondType.names.items()}
        for i, j in g.edges:
            if mol.GetBondBetweenAtoms(i, j) is not None:
                g.edges[(i,j)].update({"bond_type": bond_name_dict[mol.GetBondBetweenAtoms(i, j).GetBondType()]})
        
        for node in g.nodes:
            g.nodes[node].update({"is_aromatic": mol.GetAtomWithIdx(node).GetIsAromatic(),
                                  "symmetry_rank": rank_atoms[node]})

    return g