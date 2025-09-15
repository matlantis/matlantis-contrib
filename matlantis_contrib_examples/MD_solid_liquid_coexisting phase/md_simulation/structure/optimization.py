"""
Structure optimization functions for atomic positions and cells.

This module provides functions for optimizing atomic structures using 
various ASE optimizers and filters.
"""

from typing import Callable, Optional
from ase import Atoms
from ase.optimize import LBFGS  
from ase.filters import FrechetCellFilter


def opt_atoms(
    atoms: Atoms, 
    get_calculator: Callable,
    fmax: float = 0.01,
    steps: int = 1000,
    maxstep: float = 0.1,
    filter: bool = True,
    replace: bool = True
) -> Atoms:
    """
    Optimize atomic positions and cell using LBFGS and FrechetCellFilter.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to optimize.
    get_calculator : callable
        Function to get the calculator for atoms.
    fmax : float, optional
        Maximum force for optimization convergence (default: 0.01).
    steps : int, optional
        Maximum number of optimization steps (default: 1000).
    maxstep : float, optional
        Maximum step size (default: 0.1).
    filter : bool, optional
        Whether to use FrechetCellFilter for cell optimization (default: True).
    replace : bool, optional
        Whether to work on a copy of the atoms object (default: True).
        
    Returns
    -------
    ase.Atoms
        Optimized atoms object.
        
    Examples
    --------
    >>> from ase.calculators.emt import EMT
    >>> def get_calc():
    ...     return EMT()
    >>> optimized_atoms = opt_atoms(atoms, get_calc, fmax=0.05)
    """
    if replace:
        atoms = atoms.copy()
    
    atoms.calc = get_calculator()
    
    if filter:
        filtered_atoms = FrechetCellFilter(atoms)
        opt = LBFGS(filtered_atoms, maxstep=maxstep)
    else:
        opt = LBFGS(atoms, maxstep=maxstep)
    
    opt.run(fmax=fmax, steps=steps)
    return atoms