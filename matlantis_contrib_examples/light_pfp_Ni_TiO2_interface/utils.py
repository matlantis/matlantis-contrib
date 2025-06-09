import numpy as np
from ase import Atoms
from scipy.spatial.transform import Rotation as R
from ase.neighborlist import neighbor_list


def remove_duplicate_atoms(atoms: Atoms, threshold: float = 0.1) -> Atoms:
    i_idx, j_idx, shift = neighbor_list('ijS', atoms, threshold)
    
    valid_pairs = i_idx < j_idx
    i_clean = i_idx[valid_pairs]
    j_clean = j_idx[valid_pairs]
    shift_clean = shift[valid_pairs]
    
    keep_mask = np.ones(len(atoms), dtype=bool)
    
    for i, j, s in zip(i_clean, j_clean, shift_clean):
        if np.any(s) and keep_mask[i] and keep_mask[j]:
            keep_mask[i] = False
    
    filtered_atoms = atoms[keep_mask]
    
    return filtered_atoms


def expand_supercell(atoms: Atoms, min_spacing: float = 10.0) -> Atoms:
    cell = atoms.get_cell()
    a1, a2, a3 = cell
    
    volume = np.abs(np.dot(a1, np.cross(a2, a3)))
    
    d1 = volume / np.linalg.norm(np.cross(a2, a3))
    d2 = volume / np.linalg.norm(np.cross(a3, a1))
    d3 = volume / np.linalg.norm(np.cross(a1, a2))
    
    n1 = max(1, int(np.ceil(min_spacing / d1)))
    n2 = max(1, int(np.ceil(min_spacing / d2)))
    n3 = max(1, int(np.ceil(min_spacing / d3)))
    
    supercell = atoms.repeat((n1, n2, n3))
    return supercell


def make_cubic_cut(atoms, length=10.0, remove_duplicate=2.0):
    atoms = expand_supercell(atoms, length*2.0)
    found = False
    max_attempts = 1000
    attempt = 0
    
    while not found and attempt < max_attempts:
        attempt += 1
        
        cell = atoms.get_cell()
        frac_origin = np.random.rand(3)
        origin = np.dot(frac_origin, cell)
        
        rot = R.random()
        rotation_matrix = rot.as_matrix()
        dirs = rotation_matrix.T * length
        
        vertices = []
        for dx in [0, 1]:
            for dy in [0, 1]:
                for dz in [0, 1]:
                    vertex = origin + dx*dirs[:,0] + dy*dirs[:,1] + dz*dirs[:,2]
                    vertices.append(vertex)
        
        all_inside = True
        for v in vertices:
            scaled = np.linalg.solve(cell.T, v)
            if not np.all((scaled >= 0) & (scaled < 1)):
                all_inside = False
                break
        
        if all_inside:
            found = True
    
    if not found:
        raise RuntimeError("Cannot found approperiate position to cut")
    
    atoms.translate(-origin)
    
    atoms_new = Atoms(
        positions = atoms.get_positions().dot(rotation_matrix.T),
        cell = atoms.get_cell()[:].dot(rotation_matrix.T),
        numbers = atoms.get_atomic_numbers(),
        pbc=True,
    )
    
    pos = atoms_new.get_positions()
    mask = np.all((pos >= 0) & (pos <= length), axis=1)
    atoms_cut = atoms_new[mask]
    
    atoms_cut.set_cell(np.eye(3)*length)
    atoms_cut.wrap()
    if remove_duplicate > 0.0:
        atoms_filtered = remove_duplicate_atoms(atoms_cut, remove_duplicate)
    return atoms_filtered

