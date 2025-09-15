"""
MD Simulation Package for Melting Point Calculations

This package provides modular tools for molecular dynamics simulations
focused on determining melting temperatures through solid-liquid coexistence
simulations and local order parameter analysis.

Modules:
    core: Core simulation engine, parameters, and schedulers
    structure: Structure manipulation, optimization, and analysis tools  
    utils: Input/output and sampling utilities
"""

from .core.parameters import MDSimulationParams, SimulationSettings
from .core.simulation import run_md

# Import commonly used functions from structure module
from .structure.optimization import opt_atoms
from .structure.operations import standardize_atoms, build_supercell
from .structure.analysis import classify_solid_atoms, get_solid_prob
from .structure.cell_utils import is_triu, rotate_cell_to_triu, get_expand_cell_func

__version__ = "1.0.0"
__author__ = "MD Simulation Team"

__all__ = [
    # Core functionality
    "MDSimulationParams",
    "SimulationSettings", 
    "run_md",
    # Structure tools
    "opt_atoms",
    "standardize_atoms",
    "build_supercell", 
    "classify_solid_atoms",
    "get_solid_prob",
    "is_triu",
    "rotate_cell_to_triu",
    "get_expand_cell_func",
]