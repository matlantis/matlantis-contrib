"""
Structure operations for standardization, supercell building, and manipulation.

This module provides functions for manipulating atomic structures including
standardization using spglib, supercell construction, and density calculations.
"""

from typing import Tuple, Optional, Callable
import warnings
import numpy as np
import spglib

from ase import Atoms, units


def standardize_atoms(atoms: Atoms, symprec: float = 1e-4) -> Atoms:
    """
    Standardize the cell of an Atoms object using spglib.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to standardize.
    symprec : float, optional
        Symmetry precision for spglib (default: 1e-4).
        
    Returns
    -------
    ase.Atoms
        Standardized atoms object, or original if standardization fails.
        
    Examples
    --------
    >>> standardized = standardize_atoms(atoms, symprec=1e-5)
    """
    cell, positions, numbers = atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers
    standardized_cell_params = spglib.standardize_cell((cell, positions, numbers), symprec=symprec)
    
    if standardized_cell_params is None:
        msg = (
            'Failed to standardize ase.Atoms using spglib. Please check the structure, '
            'e.g. too many neighbors, too close atom pairs.'
        )
        warnings.warn(msg)
        return atoms
        
    cell, pos, numbers = standardized_cell_params
    atoms_standardized = Atoms(cell=cell, scaled_positions=pos, numbers=numbers, pbc=True)
    return atoms_standardized


def build_supercell(
    atoms: Atoms,
    repeat: Tuple[int, int, int] = (3, 3, 3),
    is_triu: Optional[Callable] = None,
    rotate_cell_to_triu: Optional[Callable] = None
) -> Atoms:
    """
    Build a supercell from a unit cell, sort atoms, and optionally rotate cell to upper-triangular.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to expand.
    repeat : tuple, optional
        Repetition along each axis (default: (3, 3, 3)).
    is_triu : callable, optional
        Function to check if cell is upper-triangular.
    rotate_cell_to_triu : callable, optional
        Function to rotate cell to upper-triangular.
        
    Returns
    -------
    ase.Atoms
        Supercell atoms object.
        
    Examples
    --------
    >>> supercell = build_supercell(atoms, repeat=(2, 2, 2))
    """
    supercell = atoms * repeat
    # Sort atoms by position (z, y, x order)
    indexer = sorted(range(len(supercell)), key=lambda i: supercell[i].position[::-1].tolist())
    supercell = supercell[indexer]
    
    if is_triu is not None and rotate_cell_to_triu is not None:
        cell = supercell.cell
        if not is_triu(cell):
            supercell = rotate_cell_to_triu(supercell)
            
    return supercell


def get_mass_density(atoms: Atoms, volume: Optional[float] = None, units_module=None) -> float:
    """
    Calculate the mass density of an Atoms object in g/cm^3.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object.
    volume : float, optional
        Volume to use instead of atoms.get_volume().
    units_module : module, optional
        ASE units module.
        
    Returns
    -------
    float
        Mass density in g/cm^3.
        
    Examples
    --------
    >>> density = get_mass_density(atoms)  # g/cm^3
    """
    if units_module is None:
        units_module = units
        
    total_mass = np.sum(atoms.get_masses()) * units_module._amu * 1000  # Convert to kg then g
    cell_volume = volume if volume is not None else atoms.get_volume()
    volume_cm3 = cell_volume * 1e-24  # Convert Å^3 to cm^3
    
    return total_mass / volume_cm3


def get_supercell_conversion_matrix(atoms: Atoms, target_atoms: int) -> np.ndarray:
    """
    Calculate an integer conversion matrix for ase.build.make_supercell to generate a supercell
    with a given number of atoms and lattice vectors as close to perpendicular as possible.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object (with N atoms and lattice vectors a1, a2, a3).
    target_atoms : int
        Desired number of atoms in the supercell (must be > N).
        
    Returns
    -------
    np.ndarray
        3x3 integer numpy array for make_supercell.
        
    Raises
    ------
    ValueError
        If target_atoms is not greater than the number of atoms in the input object.
        
    Examples
    --------
    >>> matrix = get_supercell_conversion_matrix(atoms, 216)
    >>> supercell = make_supercell(atoms, matrix)
    """
    N = len(atoms)
    if target_atoms <= N:
        raise ValueError("target_atoms must be greater than the number of atoms in the input Atoms object.")
    
    # Estimate the scaling factor
    scale = (target_atoms / N) ** (1/3)
    
    # Try integer multipliers around the estimated scale
    best_multipliers = None
    best_score = None
    
    for i in range(max(1, int(scale)-2), int(scale)+3):
        for j in range(max(1, int(scale)-2), int(scale)+3):
            for k in range(max(1, int(scale)-2), int(scale)+3):
                n_atoms = N * i * j * k
                if n_atoms < target_atoms:
                    continue
                    
                # Build the supercell lattice
                cell = atoms.cell @ np.diag([i, j, k])
                
                # Calculate angles between lattice vectors
                angles = []
                for m in range(3):
                    for n in range(m+1, 3):
                        dot_product = np.dot(cell[m], cell[n])
                        norm_product = np.linalg.norm(cell[m]) * np.linalg.norm(cell[n])
                        angle = np.arccos(np.clip(dot_product / norm_product, -1, 1))
                        angles.append(angle)
                
                # Score based on how close angles are to 90 degrees
                rect_score = sum(abs(np.degrees(angle) - 90) for angle in angles)
                
                # Total score: prioritize atom count close to target, then rectangularity
                score = (n_atoms - target_atoms) + rect_score
                
                if best_score is None or score < best_score:
                    best_score = score
                    best_multipliers = [i, j, k]
    
    conversion_matrix = np.diag(best_multipliers)
    return conversion_matrix