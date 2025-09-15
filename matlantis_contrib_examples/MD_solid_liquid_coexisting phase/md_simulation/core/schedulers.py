"""
Schedulers for controlling MD simulation parameters during runtime.

This module provides classes for dynamically controlling simulation parameters
such as temperature and cell deformation during MD simulations.
"""

from typing import Optional
import numpy as np
import warnings

from ase import Atoms, units
from .parameters import estimate_temp_rate, estimate_steps


class LinearTemperatureScheduler:
    """
    Controls linear temperature scheduling during MD simulation.
    
    This scheduler can be used to increase (melt) or decrease (quench) temperature 
    linearly over a specified number of MD steps. Users can specify either the
    number of steps or the temperature rate.
    
    Parameters
    ----------
    dyn : object
        MD dynamics object (ASE MD class instance).
    temp_start : float
        Starting temperature in K.
    temp_end : float
        Ending temperature in K.
    nsteps : int, optional
        Number of steps for temperature change. Ignored if temp_rate is provided.
    temp_rate : float, optional
        Temperature rate in K/fs. If provided, nsteps is calculated automatically.
    timestep_fs : float, optional
        Timestep in femtoseconds. If None, attempts to get from dyn.timestep.
    start_step : int, optional
        MD step at which to start temperature scheduling (default: 0).
    end_step : int, optional
        MD step at which to end temperature scheduling. If None, calculated from nsteps.
        
    Examples
    --------
    >>> from ase.md.npt import NPT
    >>> dyn = NPT(atoms, timestep=2.0*units.fs, temperature_K=300)
    >>> scheduler = LinearTemperatureScheduler(
    ...     dyn, temp_start=300, temp_end=1500, nsteps=10000
    ... )
    >>> scheduler.attach(interval=1)  # Update temperature every step
    >>> dyn.run(10000)
    
    >>> # Or using temperature rate
    >>> scheduler = LinearTemperatureScheduler(
    ...     dyn, temp_start=300, temp_end=1500, temp_rate=0.12
    ... )
    """
    
    def __init__(
        self,
        dyn,
        temp_start: float,
        temp_end: float,
        nsteps: Optional[int] = None,
        temp_rate: Optional[float] = None,
        timestep_fs: Optional[float] = None,
        start_step: int = 0,
        end_step: Optional[int] = None,
    ):
        self.dyn = dyn
        self.temp_start = temp_start
        self.temp_end = temp_end
        self.start_step = start_step
        
        # Try to get timestep from dynamics object if not provided
        if timestep_fs is None or timestep_fs == 0.0:
            if hasattr(dyn, "timestep"):
                try:
                    timestep_fs = float(dyn.timestep / units.fs)
                except Exception:
                    timestep_fs = 1.0
            else:
                timestep_fs = 1.0
        self.timestep_fs = timestep_fs
        
        # Calculate nsteps and temperature increment
        if temp_rate is not None:
            if nsteps is not None:
                warnings.warn("Both temp_rate and nsteps are provided. nsteps will be ignored.")
            self.nsteps = estimate_steps(self.temp_start, self.temp_end, temp_rate, timestep_fs)
            self.dT = temp_rate * timestep_fs  # Temperature change per MD step
        elif nsteps is not None:
            self.nsteps = nsteps
            temp_rate_fs = estimate_temp_rate(self.temp_start, self.temp_end, nsteps, timestep_fs)
            self.dT = temp_rate_fs * timestep_fs  # Temperature change per MD step
        else:
            raise ValueError("Either nsteps or temp_rate must be specified.")
        
        # Set end step
        if end_step is None:
            self.end_step = self.start_step + self.nsteps
        else:
            self.end_step = end_step
    
    def __call__(self) -> None:
        """
        Update the temperature of the MD simulation at each step.
        
        This method is called automatically when attached to the MD dynamics.
        Only updates temperature within the specified step range.
        """
        current_md_step = self.dyn.nsteps
        
        if self.start_step <= current_md_step < self.end_step:
            # Linear temperature interpolation
            progress = current_md_step - self.start_step
            current_temp = self.temp_start + self.dT * progress
            
            # Ensure we don't overshoot
            if self.temp_start < self.temp_end:
                current_temp = min(current_temp, self.temp_end)
            else:
                current_temp = max(current_temp, self.temp_end)
            
            self.dyn.set_temperature(current_temp)
    
    def update_temperature(self):
        """Alias for __call__ method."""
        self.__call__()
    
    def attach(self, interval: int = 1) -> None:
        """
        Attach the temperature scheduler to the MD dynamics object.
        
        Parameters
        ----------
        interval : int, optional
            How often to update the temperature in MD steps (default: 1).
        """
        self.dyn.attach(self.update_temperature, interval=interval)
    
    def get_current_temperature(self) -> float:
        """
        Get the current target temperature based on current MD step.
        
        Returns
        -------
        float
            Current target temperature in K.
        """
        current_md_step = self.dyn.nsteps
        
        if current_md_step < self.start_step:
            return self.temp_start
        elif current_md_step >= self.end_step:
            return self.temp_end
        else:
            progress = current_md_step - self.start_step
            return self.temp_start + self.dT * progress


class CellDeformationScheduler:
    """
    Scheduler to gradually deform the simulation cell to a target cell during MD.
    
    This scheduler linearly interpolates between the initial cell and target cell
    over a specified number of MD steps, enabling gradual cell deformation.
    
    Parameters
    ----------
    dyn : object
        MD dynamics object.
    atoms : ase.Atoms
        The atoms object whose cell will be deformed.
    target_cell : array-like
        The target cell matrix (3x3).
    nsteps : int
        Number of MD steps over which to perform the deformation.
    start_step : int, optional
        The MD step at which to start the deformation (default: 0).
    end_step : int, optional
        The MD step at which to end the deformation. If None, calculated as
        start_step + nsteps.
        
    Examples
    --------
    >>> target_cell = [[10, 0, 0], [0, 10, 0], [0, 0, 15]]  # Stretch z-axis
    >>> scheduler = CellDeformationScheduler(
    ...     dyn, atoms, target_cell, nsteps=5000
    ... )
    >>> dyn.attach(scheduler, interval=1)
    >>> dyn.run(10000)
    """
    
    def __init__(
        self, 
        dyn, 
        atoms: Atoms, 
        target_cell: np.ndarray, 
        nsteps: int, 
        start_step: int = 0, 
        end_step: Optional[int] = None
    ):
        self.dyn = dyn
        self.atoms = atoms
        self.target_cell = np.array(target_cell)
        self.start_step = start_step
        self.end_step = end_step if end_step is not None else start_step + nsteps
        self.nsteps = min(nsteps, self.end_step - self.start_step)
        self.initial_cell = atoms.get_cell().copy()
    
    def __call__(self):
        """
        Deform the cell at each MD step.
        
        This method linearly interpolates between initial and target cell
        based on the current MD step progress.
        """
        current_md_step = self.dyn.nsteps
        
        if self.start_step <= current_md_step < self.end_step:
            # Calculate interpolation parameter (0 to 1)
            progress = (current_md_step - self.start_step) / self.nsteps
            progress = min(progress, 1.0)  # Ensure we don't overshoot
            
            # Linear interpolation between initial and target cell
            new_cell = (1 - progress) * self.initial_cell + progress * self.target_cell
            self.atoms.set_cell(new_cell, scale_atoms=True)
            
        elif current_md_step >= self.end_step:
            # Ensure final cell is exactly the target
            self.atoms.set_cell(self.target_cell, scale_atoms=True)
    
    def get_current_cell(self) -> np.ndarray:
        """
        Get the current target cell based on current MD step.
        
        Returns
        -------
        np.ndarray
            Current target cell matrix.
        """
        current_md_step = self.dyn.nsteps
        
        if current_md_step < self.start_step:
            return self.initial_cell.copy()
        elif current_md_step >= self.end_step:
            return self.target_cell.copy()
        else:
            progress = (current_md_step - self.start_step) / self.nsteps
            progress = min(progress, 1.0)
            return (1 - progress) * self.initial_cell + progress * self.target_cell