"""
Core MD simulation engine and high-level simulation functions.

This module provides the main MD simulation runner and higher-level functions
for common simulation workflows like melting, equilibration, and coexistence simulations.
"""

import os
from typing import List, Tuple
from pathlib import Path
from copy import deepcopy

import numpy as np
from joblib import Parallel, delayed

from ase import Atoms, units
from ase.io import read
from ase.io.trajectory import Trajectory
from ase.md.npt import NPT
from ase.md.nptberendsen import NPTBerendsen, Inhomogeneous_NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary, ZeroRotation
from ase.md import MDLogger

from .parameters import MDSimulationParams
from .schedulers import LinearTemperatureScheduler


def run_md(atoms: Atoms, params: MDSimulationParams) -> Tuple[Atoms, MDSimulationParams]:
    """
    Run a single MD simulation for a given Atoms object using specified parameters.
    
    This is the core MD simulation function that handles setup, execution, and
    cleanup for molecular dynamics simulations. It supports both NPT and NVT
    ensembles with optional temperature ramping.
    
    Parameters
    ----------
    atoms : ase.Atoms
        The atoms object to simulate.
    params : MDSimulationParams
        Dataclass containing all MD simulation parameters.
        
    Returns
    -------
    tuple
        (final_atoms, updated_params) where final_atoms is the atoms object
        after simulation and updated_params contains any parameter updates.
        
    Examples
    --------
    >>> from ase.calculators.emt import EMT
    >>> def get_calc():
    ...     return EMT()
    >>> params = MDSimulationParams(
    ...     get_calculator=get_calc,
    ...     num_steps=1000,
    ...     temperature=1000.0
    ... )
    >>> final_atoms, params = run_md(atoms, params)
    """
    # Constants
    amu = getattr(units, "amu", 1.66053906660e-27)
    
    # Prepare atoms
    atoms = atoms.copy()
    atoms.calc = params.get_calculator()
    
    # Convert parameters to ASE units
    timestep_fs = params.timestep * units.fs
    pressure = params.pressure * units.bar
    temperature = params.temperature
    temp_begin = params.temp_begin
    temp_end = params.temp_end
    pfactor_val = params.pfactor * units.GPa * (units.fs**2)
    ttime_fs = params.ttime * units.fs
    taut = params.taut * units.fs
    taup = params.taup * units.fs
    compressibility = params.compressibility / units.GPa
    fixcm = getattr(params, "fixcm", True)
    mask = getattr(params, "mask", (1, 1, 1))
    
    # File handling
    trajfile = params.trajfile or f"md_{int(params.temperature)}K.traj"
    logfile = params.logfile or f"md_{int(params.temperature)}K.log"
    loginterval = params.loginterval
    logger = params.logger
    logger_interval = params.logger_interval
    traj_interval = params.traj_interval
    ensemble = params.ensemble
    overwrite = params.overwrite
    check_density = params.check_density
    attach_functions = params.attach_functions or []
    
    # Check if simulation already exists and skip if not overwriting
    if not overwrite and os.path.exists(trajfile):
        traj = Trajectory(trajfile)
        final_atoms = traj[-1]
        return final_atoms, params
    
    # Initialize velocities
    MaxwellBoltzmannDistribution(atoms, temperature_K=params.temperature)
    Stationary(atoms)
    ZeroRotation(atoms)
    
    # Set up MD dynamics
    if ensemble.lower() == "npt":
        dyn = Inhomogeneous_NPTBerendsen(
            atoms,
            timestep_fs,
            temperature_K=temperature,
            pressure_au=pressure,
            taut=taut,
            taup=taup,
            compressibility_au=compressibility,
            fixcm=fixcm,
            logfile=logfile,
            loginterval=loginterval,
            mask=mask,
        )
    elif ensemble.lower() == "nvt":
        dyn = NVTBerendsen(
            atoms,
            timestep_fs,
            temperature_K=temperature,
            taut=taut,
            logfile=logfile,
            loginterval=loginterval,
        )
    else:
        raise ValueError(f"Unknown ensemble: {ensemble}. Only 'npt' and 'nvt' are supported.")
    
    # Set up trajectory recording
    traj = Trajectory(trajfile, mode="w", atoms=atoms)
    dyn.attach(traj.write, interval=traj_interval)
    
    # Set up logging
    if logger:
        dyn.attach(
            MDLogger(dyn, atoms, logger, header=False, stress=False, peratom=True, mode="a"),
            interval=logger_interval
        )
    
    # Set up temperature scheduling if specified
    if temp_begin is not None and temp_end is not None:
        temperature_scheduler = LinearTemperatureScheduler(
            dyn, 
            temp_start=temp_begin, 
            temp_end=temp_end, 
            timestep_fs=params.timestep,
            nsteps=params.num_steps
        )
        dyn.attach(temperature_scheduler, interval=1)
    
    # Attach additional functions
    for func, interval in attach_functions:
        dyn.attach(func, interval=interval)
    
    # Check initial density if requested
    if check_density:
        rho_initial = np.sum(atoms.get_masses()) * amu * 1000 / (atoms.get_volume() * 1e-24)
    
    # Run simulation
    dyn.run(params.num_steps)
    
    # Final processing
    atoms = dyn.atoms
    atoms.wrap()
    
    # Check final density if requested
    if check_density:
        rho_final = np.sum(atoms.get_masses()) * amu * 1000 / (atoms.get_volume() * 1e-24)
        if rho_final / rho_initial < 0.1:  # System vaporized
            raise ValueError(
                f"The density of final state ({rho_final:.4e} g/cm³) is much lower "
                f"than initial ({rho_initial:.4e} g/cm³). System may have vaporized."
            )
    
    return atoms, params


def run_melt(
    atomses: List[Atoms],
    get_calculator,
    temperatures: List[float],
    num_steps: int,
    timestep: float,
    stress: float,
    pfactor: float,
    ttime: float,
    traj_interval: int = 100,
    loginterval: int = 100,
    check_density: bool = True,
    n_jobs: int = -1,
) -> List[Tuple[Atoms, MDSimulationParams]]:
    """
    Run melting MD simulations for a list of atoms at different temperatures in parallel.
    
    Parameters
    ----------
    atomses : list of ase.Atoms
        List of atoms objects to simulate.
    get_calculator : callable
        Function that returns an ASE calculator.
    temperatures : list of float
        List of temperatures for each atoms object.
    num_steps : int
        Number of MD steps to run.
    timestep : float
        MD timestep in fs.
    stress : float
        External stress in bar.
    pfactor : float
        Pressure factor.
    ttime : float
        Time parameter in fs.
    traj_interval : int, optional
        Trajectory writing interval (default: 100).
    loginterval : int, optional
        Logging interval (default: 100).
    check_density : bool, optional
        Whether to check for vaporization (default: True).
    n_jobs : int, optional
        Number of parallel jobs (default: -1 for all cores).
        
    Returns
    -------
    list
        List of (atoms, params) tuples from each simulation.
    """
    params_list = []
    for atoms, temp in zip(atomses, temperatures):
        trajfile = f"melt_{int(temp)}K.traj"
        logfile = f"melt_{int(temp)}K.log"
        params = MDSimulationParams(
            get_calculator=get_calculator,
            num_steps=num_steps,
            timestep=timestep,
            temperature=temp,
            pressure=stress,  # Note: using stress as pressure
            pfactor=pfactor,
            ttime=ttime,
            trajfile=trajfile,
            logfile=logfile,
            loginterval=loginterval,
            traj_interval=traj_interval,
            ensemble="npt",
            overwrite=False,
            check_density=check_density,
        )
        params_list.append((atoms, params))
    
    # Run in parallel
    parallel = Parallel(n_jobs=n_jobs, backend="threading")
    results = parallel(delayed(run_md)(atoms, params) for atoms, params in params_list)
    return results


def run_equilibrate(
    atoms: Atoms,
    get_calculator,
    num_steps: int,
    timestep: float,
    temperature: float,
    stress: float,
    pfactor: float,
    ttime: float,
    traj_interval: int = 100,
    loginterval: int = 100,
    overwrite: bool = False,
) -> Tuple[Atoms, MDSimulationParams]:
    """
    Run equilibration MD for a single atoms object.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Atoms object to equilibrate.
    get_calculator : callable
        Function that returns an ASE calculator.
    num_steps : int
        Number of equilibration steps.
    timestep : float
        MD timestep in fs.
    temperature : float
        Equilibration temperature in K.
    stress : float
        External stress in bar.
    pfactor : float
        Pressure factor.
    ttime : float
        Time parameter in fs.
    traj_interval : int, optional
        Trajectory writing interval (default: 100).
    loginterval : int, optional
        Logging interval (default: 100).
    overwrite : bool, optional
        Whether to overwrite existing files (default: False).
        
    Returns
    -------
    tuple
        (equilibrated_atoms, params) tuple.
    """
    outfile = "equilib_initial.json"
    logfile = f"equilib_initial_{int(temperature)}K.log"
    trajfile = logfile.replace(".log", ".traj")
    
    if not overwrite and os.path.exists(outfile):
        atoms = read(outfile)
        return atoms, None
    
    params = MDSimulationParams(
        get_calculator=get_calculator,
        num_steps=num_steps,
        timestep=timestep,
        temperature=temperature,
        pressure=stress,  # Note: using stress as pressure
        pfactor=pfactor,
        ttime=ttime,
        trajfile=trajfile,
        logfile=logfile,
        loginterval=loginterval,
        traj_interval=traj_interval,
        ensemble="npt",
        overwrite=overwrite,
        check_density=True,
    )
    
    atoms, params = run_md(atoms, params)
    atoms.write(outfile)
    return atoms, params


def run_md_solid(
    atoms: Atoms,
    get_calculator,
    num_steps: int,
    timestep: float,
    temperature: float,
    temp_begin: float,
    temp_end: float,
    stress: float,
    pfactor: float,
    ttime: float,
    traj_interval: int = 100,
    loginterval: int = 100,
    overwrite: bool = False,
) -> Tuple[Atoms, MDSimulationParams]:
    """
    Run MD for solid phase with temperature ramping.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Solid atoms object to simulate.
    get_calculator : callable
        Function that returns an ASE calculator.
    num_steps : int
        Number of MD steps.
    timestep : float
        MD timestep in fs.
    temperature : float
        Base temperature in K.
    temp_begin : float
        Starting temperature for ramping in K.
    temp_end : float
        Ending temperature for ramping in K.
    stress : float
        External stress in bar.
    pfactor : float
        Pressure factor.
    ttime : float
        Time parameter in fs.
    traj_interval : int, optional
        Trajectory writing interval (default: 100).
    loginterval : int, optional
        Logging interval (default: 100).
    overwrite : bool, optional
        Whether to overwrite existing files (default: False).
        
    Returns
    -------
    tuple
        (final_atoms, params) tuple.
    """
    trajfile = f"solid_{int(temp_begin)}K_{int(temp_end)}K.traj"
    logfile = f"solid_{int(temp_begin)}K_{int(temp_end)}K.log"
    
    params = MDSimulationParams(
        get_calculator=get_calculator,
        num_steps=num_steps,
        timestep=timestep,
        temperature=temperature,
        temp_begin=temp_begin,
        temp_end=temp_end,
        pressure=stress,  # Note: using stress as pressure
        pfactor=pfactor,
        ttime=ttime,
        trajfile=trajfile,
        logfile=logfile,
        loginterval=loginterval,
        traj_interval=traj_interval,
        ensemble="npt",
        overwrite=overwrite,
        check_density=True,
    )
    
    return run_md(atoms, params)