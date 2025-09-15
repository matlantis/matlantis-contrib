"""
Structure analysis functions for order parameter calculation and phase classification.

This module provides functions for analyzing atomic structures including solid/liquid
classification using local order parameters and bulk modulus calculations.
"""

from typing import Callable, Optional, Union, List, Tuple
import numpy as np

from ase import Atoms, units
from ase.optimize import LBFGS
from ase.eos import EquationOfState as EOS
from pyscal3 import System


def classify_solid_atoms(
    atoms: Atoms,
    cutoff: Union[str, float] = 'adaptive',
    padding: float = 1.2,
    q: int = 6,
    threshold: float = 0.5,
    avgthreshold: float = 0.6,
    cluster: bool = True,
    method: str = 'pyscal3',
    **kwargs
) -> np.ndarray:
    """
    Identify solid atoms in the given ASE Atoms object using local order parameters.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object to analyze.
    cutoff : str or float, optional
        Cutoff method for neighbor finding ('adaptive', 'sann', or float value).
        Default: 'adaptive'.
    padding : float, optional
        Additional padding for neighbor search (default: 1.2).
    q : int, optional
        Steinhardt order parameter (default: 6).
    threshold : float, optional
        Local order parameter threshold (default: 0.5).
    avgthreshold : float, optional
        Average order parameter threshold (default: 0.6).
    cluster : bool, optional
        Whether to use clustering for solid detection (default: True).
    method : str, optional
        Classification method ('pyscal3', 'other'). Default: 'pyscal3'.
    **kwargs
        Additional arguments for future methods.
        
    Returns
    -------
    np.ndarray
        Boolean array indicating whether each atom is solid.
        
    Raises
    ------
    ValueError
        If an unknown classification method is specified.
    NotImplementedError
        If an alternative method is requested but not implemented.
        
    Examples
    --------
    >>> solid_mask = classify_solid_atoms(atoms, q=6, threshold=0.5)
    >>> solid_fraction = np.mean(solid_mask)
    """
    if method == 'pyscal3':
        pcsys = System()
        pcsys.read.file(atoms, 'ase')
        
        try:
            if isinstance(cutoff, str):
                pcsys.find.neighbors(method='cutoff', cutoff=cutoff, padding=padding)
            else:
                pcsys.find.neighbors(method='cutoff', cutoff=cutoff, padding=padding)
        except RuntimeError as e:
            print(f"Neighbor finding failed with cutoff={cutoff}, falling back to 'sann'. Exception: {e}")
            pcsys.find.neighbors(method='cutoff', cutoff='sann')
        
        # Find solid atoms using Steinhardt order parameters
        pcsys.find.solids(q=q, threshold=threshold, avgthreshold=avgthreshold, cluster=cluster)
        return pcsys.atoms.solid
        
    elif method == 'other':
        # Placeholder for future alternative classification packages
        raise NotImplementedError("Alternative solid classification method not implemented yet.")
    else:
        raise ValueError(f"Unknown classification method: {method}")


def get_solid_prob(
    atoms: Atoms, 
    q: int = 6, 
    atomic_numbers: Optional[List[int]] = None, 
    **kwargs
) -> np.ndarray:
    """
    Calculate the probability that atoms are in a solid-like environment using local order analysis.
    
    Parameters
    ----------
    atoms : ase.Atoms
        ASE Atoms object to analyze.
    q : int, optional
        Steinhardt order parameter (default: 6).
    atomic_numbers : list of int, optional
        Optional list of atomic numbers to filter atoms.
    **kwargs
        Additional arguments passed to classify_solid_atoms.
        
    Returns
    -------
    np.ndarray
        Array of solid probabilities for the selected atoms.
        
    Examples
    --------
    >>> prob = get_solid_prob(atoms, q=6)
    >>> mean_solid_prob = np.mean(prob)
    """
    if atomic_numbers:
        mask = np.isin(atoms.numbers, atomic_numbers)
        filtered_atoms = atoms[mask]
    else:
        filtered_atoms = atoms
    
    solid_prob = classify_solid_atoms(
        filtered_atoms,
        q=q,
        **kwargs
    )
    return solid_prob


def calculate_bulk_modulus(
    atoms: Atoms,
    volume_ratio_lim: Tuple[float, float],
    get_calculator: Callable,
    num_data_points: int = 10,
    fmax: float = 0.005,
    steps: int = 200,
    maxstep: float = 0.1,
    optimizer=LBFGS,
    return_eos: bool = True
) -> Union[Tuple[float, float, float], Tuple[float, float, float, EOS]]:
    """
    Fit the Birch-Murnaghan equation of state to volume-energy data and return the bulk modulus.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object to analyze.
    volume_ratio_lim : tuple of float
        (min_ratio, max_ratio) for volume scaling.
    get_calculator : callable
        Function that returns an ASE calculator.
    num_data_points : int, optional
        Number of volume points to sample (default: 10).
    fmax : float, optional
        Force convergence criterion for optimization (default: 0.005).
    steps : int, optional
        Maximum optimization steps (default: 200).
    maxstep : float, optional
        Maximum step size for optimization (default: 0.1).
    optimizer : class, optional
        ASE optimizer class (default: LBFGS).
    return_eos : bool, optional
        Whether to return the EOS object (default: True).
        
    Returns
    -------
    tuple
        If return_eos=False: (equilibrium_volume, equilibrium_energy, bulk_modulus_GPa)
        If return_eos=True: (equilibrium_volume, equilibrium_energy, bulk_modulus_GPa, eos)
        
    Examples
    --------
    >>> from ase.calculators.emt import EMT
    >>> def get_calc():
    ...     return EMT()
    >>> V0, E0, B = calculate_bulk_modulus(atoms, (0.9, 1.1), get_calc)
    >>> print(f"Bulk modulus: {B:.2f} GPa")
    """
    atoms = atoms.copy()
    atoms.calc = get_calculator()
    
    original_volume = atoms.get_volume()
    v0, v1 = volume_ratio_lim
    volumes = original_volume * np.linspace(v0, v1, num_data_points)
    energies = np.zeros(num_data_points, dtype=float)
    
    for i, target_volume in enumerate(volumes):
        # Scale cell to target volume
        current_volume = atoms.get_volume()
        scale_factor = (target_volume / current_volume) ** (1/3)
        atoms.set_cell(atoms.get_cell() * scale_factor, scale_atoms=True)
        
        # Optimize structure at this volume
        opt = optimizer(atoms, maxstep=maxstep)
        opt.run(fmax=fmax, steps=steps)
        energies[i] = atoms.get_potential_energy()
    
    # Fit equation of state
    eos = EOS(volumes, energies, eos="birchmurnaghan")
    equilibrium_volume, equilibrium_energy, bulk_modulus_au = eos.fit()
    bulk_modulus_GPa = bulk_modulus_au / units.GPa
    
    result = (equilibrium_volume, equilibrium_energy, bulk_modulus_GPa)
    if return_eos:
        result += (eos,)
    
    return result