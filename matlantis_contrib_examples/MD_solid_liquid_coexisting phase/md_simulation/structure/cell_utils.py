"""
Cell manipulation utilities for lattice transformations and scaling.

This module provides functions for manipulating unit cells including
triangularization, expansion based on temperature-volume relationships,
and atom rescaling operations.
"""

from typing import Callable, Optional
import numpy as np
import pandas as pd
from scipy.stats import linregress

from ase import Atoms
from ase.optimize import LBFGS
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary


def is_triu(cell: np.ndarray) -> bool:
    """
    Check if a cell matrix is upper-triangular.
    
    Parameters
    ----------
    cell : np.ndarray
        Cell matrix to check.
        
    Returns
    -------
    bool
        True if upper-triangular, else False.
        
    Examples
    --------
    >>> cell = np.array([[1, 0.5, 0.2], [0, 1, 0.3], [0, 0, 1]])
    >>> is_triu(cell)
    True
    """
    return np.allclose(cell, np.triu(cell))


def rotate_cell_to_triu(atoms: Atoms) -> Atoms:
    """
    Rotate the cell of an Atoms object to be upper-triangular.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to rotate.
        
    Returns
    -------
    ase.Atoms
        Atoms object with rotated cell.
        
    Examples
    --------
    >>> rotated_atoms = rotate_cell_to_triu(atoms)
    """
    atoms = atoms.copy()
    
    # First rotate to align z-axis
    atoms.rotate((0, 0, 1), 'z', rotate_cell=True)
    
    # Then rotate in xy-plane to make cell upper-triangular
    x, y = atoms.cell[1, :2]
    r = np.sqrt(x**2 + y**2)
    
    if r > 1e-10:  # Avoid division by zero
        theta = np.sign(x) * np.arccos(np.clip(y / r, -1, 1)) / np.pi * 180
        atoms.rotate(theta, 'z', rotate_cell=True)
    
    # Round to remove numerical noise
    cell = np.round(atoms.cell, 8)
    atoms.set_cell(cell)
    
    return atoms


def get_expand_cell_func(
    natoms: int, 
    temperatures: np.ndarray, 
    volumes: np.ndarray, 
    window: int = 1, 
    step: int = 10
) -> Callable[[Atoms, float, Optional[int]], Atoms]:
    """
    Create a function to expand unit cells based on temperature-volume relationship.
    
    This function fits a linear relationship between temperature and volume
    from MD simulation data, then returns a function that can expand cells
    to the appropriate volume for a given temperature.
    
    Parameters
    ----------
    natoms : int
        Number of atoms in the reference system.
    temperatures : np.ndarray
        Array of temperatures from MD simulations.
    volumes : np.ndarray
        Array of corresponding volumes from MD simulations.
    window : int, optional
        Rolling window size for smoothing (default: 1).
    step : int, optional
        Step size for rolling window (default: 10).
        
    Returns
    -------
    callable
        Function that expands atoms to target temperature.
        Function signature: expand_cell(atoms, temperature, axis=None)
        
    Examples
    --------
    >>> expand_func = get_expand_cell_func(108, temps, volumes)
    >>> expanded_atoms = expand_func(atoms, 1500.0)  # Expand to 1500K volume
    """
    # Create temperature-volume relationship
    df = pd.DataFrame({'T': temperatures, 'V': volumes})
    df_mean = df.rolling(window, center=True, step=step).mean().dropna()
    
    # Fit linear relationship
    fit_result = linregress(df_mean['T'], df_mean['V'])
    
    def expand_cell(
        atoms: Atoms, 
        temperature: float = 2500, 
        axis: Optional[int] = None
    ) -> Atoms:
        """
        Expand the cell of an atoms object to match the volume at a target temperature.
        
        Parameters
        ----------
        atoms : ase.Atoms
            Atoms object to expand.
        temperature : float, optional
            Target temperature (default: 2500).
        axis : int, optional
            Axis along which to expand. If None, expand isotropically.
            
        Returns
        -------
        ase.Atoms
            Atoms with expanded cell.
        """
        current_volume = atoms.get_volume()
        
        # Calculate target volume based on linear T-V relationship
        volume_per_atom_target = temperature * fit_result.slope + fit_result.intercept
        target_volume = volume_per_atom_target * len(atoms) / natoms
        
        if axis is None:
            # Isotropic expansion
            scale_factor = (target_volume / current_volume) ** (1/3)
            new_cell = atoms.cell * scale_factor
        else:
            # Expand along specific axis
            scale_factor = target_volume / current_volume
            new_cell = atoms.cell.copy()
            new_cell[axis] *= scale_factor
            
        atoms.set_cell(new_cell, scale_atoms=True)
        return atoms
    
    return expand_cell


def rescale_atoms(
    atoms: Atoms, 
    temperature: float, 
    cell_initial: np.ndarray, 
    opt_steps: int = 10
) -> Atoms:
    """
    Rescale atoms based on volume expansion and re-equilibrate.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object to rescale.
    temperature : float
        Target temperature for Maxwell-Boltzmann distribution.
    cell_initial : np.ndarray
        Initial cell matrix for reference.
    opt_steps : int, optional
        Number of optimization steps (default: 10).
        
    Returns
    -------
    ase.Atoms
        Rescaled and equilibrated atoms object.
        
    Examples
    --------
    >>> rescaled = rescale_atoms(atoms, 1000.0, initial_cell)
    """
    atoms = atoms.copy()
    
    # Store scaled positions
    scaled_positions = atoms.get_scaled_positions()
    
    # Calculate volume expansion
    volume_initial = np.linalg.det(cell_initial)
    expansion = atoms.get_volume() / volume_initial
    
    # Apply expansion along z-axis
    new_cell = cell_initial.copy()
    new_cell[2] *= expansion
    
    atoms.set_cell(new_cell)
    atoms.set_scaled_positions(scaled_positions)
    
    # Brief optimization
    if hasattr(atoms, 'calc') and atoms.calc is not None:
        opt = LBFGS(atoms)
        opt.run(steps=opt_steps)
    
    # Set Maxwell-Boltzmann velocities
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature)
    Stationary(atoms)
    
    return atoms