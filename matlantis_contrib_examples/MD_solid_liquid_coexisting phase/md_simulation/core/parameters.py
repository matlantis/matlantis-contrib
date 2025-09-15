"""
Parameter management for MD simulations.

This module provides dataclasses and utilities for managing simulation parameters,
including parameter validation, temperature rate calculations, and settings management.
"""

from dataclasses import dataclass
from typing import Optional, Callable, Tuple, List, Dict, Any


@dataclass
class MDSimulationParams:
    """
    Dataclass for holding all parameters required for a molecular dynamics simulation.
    
    This class centralizes all MD simulation parameters with sensible defaults
    and provides a clean interface for parameter management.
    
    Parameters
    ----------
    get_calculator : Callable
        Function that returns an ASE calculator for the atoms.
    num_steps : int
        Number of MD steps to run (default: 0).
    timestep : float
        Time step for the MD simulation in fs (default: 1.0).
    temperature : float
        Temperature at which to run the simulation in K (default: 298.15).
    temperatures : tuple
        Tuple of temperatures for multi-temperature simulations.
    temp_begin : float, optional
        Starting temperature for temperature ramping.
    temp_end : float, optional
        Ending temperature for temperature ramping.
    pfactor : float
        Pressure factor for the simulation (default: 2e6).
    ttime : float
        Time scale parameter in fs (default: 20.0).
    pressure : float
        External pressure in bar (default: 1.0).
    taut : float
        Thermostat time constant in fs (default: 500.0).
    taup : float
        Barostat time constant in fs (default: 1000.0).
    compressibility : float
        Compressibility in 1/GPa (default: 0.005).
    fixcm : bool
        Whether to fix center of mass motion (default: True).
    mask : tuple
        Mask for cell degrees of freedom (default: (1, 1, 1)).
    trajfile : str, optional
        Path to trajectory file. If None, a default name is used.
    logfile : str, optional
        Path to log file. If None, a default name is used.
    logger : str, optional
        Path to the logging file. If None, a default name is used.
    loginterval : int
        Interval for logging (default: 100).
    traj_interval : int
        Interval for writing trajectory (default: 100).
    logger_interval : int
        Interval for logger output (default: 100).
    ensemble : str
        Type of ensemble to use ('npt' or 'nvt') (default: 'npt').
    overwrite : bool
        Whether to overwrite existing files (default: False).
    check_density : bool
        Whether to check density during the simulation (default: True).
    attach_functions : list, optional
        List of (function, interval) tuples to attach to MD dynamics.
        
    Examples
    --------
    >>> from ase.calculators.emt import EMT
    >>> def get_calc():
    ...     return EMT()
    >>> params = MDSimulationParams(
    ...     get_calculator=get_calc,
    ...     num_steps=10000,
    ...     temperature=1000.0,
    ...     timestep=2.0
    ... )
    """
    
    get_calculator: Callable = None
    num_steps: int = 0
    timestep: float = 1.0
    temperature: float = 298.15
    temperatures: Tuple[float, ...] = (280, 320, 360, 400, 440)
    temp_begin: Optional[float] = None
    temp_end: Optional[float] = None
    pfactor: float = 2e6
    ttime: float = 20.0
    pressure: float = 1.0  # bar
    taut: float = 500.0  # fs
    taup: float = 1000.0  # fs
    compressibility: float = 0.005  # 1/GPa
    fixcm: bool = True
    mask: Tuple[int, int, int] = (1, 1, 1)
    trajfile: Optional[str] = None
    logfile: Optional[str] = None
    logger: Optional[str] = None
    loginterval: int = 100
    traj_interval: int = 100
    logger_interval: int = 100
    ensemble: str = "npt"
    overwrite: bool = False
    check_density: bool = True
    attach_functions: Optional[List[Tuple[Callable, int]]] = None
    
    def __post_init__(self):
        """Validate parameters after initialization."""
        if self.get_calculator is None:
            raise ValueError("get_calculator must be provided")
        if self.num_steps < 0:
            raise ValueError("num_steps must be non-negative")
        if self.timestep <= 0:
            raise ValueError("timestep must be positive")
        if self.temperature <= 0:
            raise ValueError("temperature must be positive")
        if self.ensemble.lower() not in ['npt', 'nvt']:
            raise ValueError("ensemble must be 'npt' or 'nvt'")


class SimulationSettings:
    """
    Helper class to manage simulation settings with flexible parameter handling.
    
    This class allows parameters to be set from an input dictionary (e.g., from file)
    or directly via keyword arguments, with keyword arguments taking precedence.
    
    Parameters
    ----------
    input_dict : dict, optional
        Dictionary of settings loaded from file or other source.
    **kwargs
        Additional settings provided as keyword arguments.
        
    Examples
    --------
    >>> settings_dict = {'temperature': 1000, 'num_steps': 5000}
    >>> settings = SimulationSettings(settings_dict, temperature=1200)
    >>> settings.get('temperature')  # Returns 1200 (kwargs override)
    1200
    >>> settings.get('num_steps')    # Returns 5000 (from input_dict)
    5000
    """
    
    def __init__(self, input_dict: Optional[Dict[str, Any]] = None, **kwargs):
        self._input_dict = input_dict or {}
        self._kwargs = kwargs
    
    def get(self, key: str, default: Any = None) -> Any:
        """
        Get a setting value with fallback priority: kwargs > input_dict > default.
        
        Parameters
        ----------
        key : str
            Setting key to retrieve.
        default : Any, optional
            Default value if key is not found.
            
        Returns
        -------
        Any
            The setting value.
        """
        if key in self._kwargs:
            return self._kwargs[key]
        if key in self._input_dict:
            return self._input_dict[key]
        return default
    
    def as_dict(self) -> Dict[str, Any]:
        """Return merged dictionary with kwargs overriding input_dict."""
        result = dict(self._input_dict)
        result.update(self._kwargs)
        return result
    
    def __getitem__(self, key: str) -> Any:
        """Dictionary-style access."""
        if key in self._kwargs:
            return self._kwargs[key]
        return self._input_dict[key]
    
    def __contains__(self, key: str) -> bool:
        """Check if key exists in either kwargs or input_dict."""
        return key in self._kwargs or key in self._input_dict
    
    def keys(self):
        """Return all available keys."""
        return set(self._input_dict.keys()).union(self._kwargs.keys())
    
    def items(self):
        """Return all key-value pairs."""
        return self.as_dict().items()
    
    def __repr__(self) -> str:
        return f"SimulationSettings({self.as_dict()})"


def estimate_temp_rate(T_start: float, T_end: float, nsteps: int, timestep_fs: float = 1.0) -> float:
    """
    Estimate the temperature rate (K/fs) required to go from T_start to T_end in nsteps.
    
    Parameters
    ----------
    T_start : float
        Initial temperature in K.
    T_end : float
        Final temperature in K.
    nsteps : int
        Number of MD steps.
    timestep_fs : float, optional
        Timestep in femtoseconds (default: 1.0).
        
    Returns
    -------
    float
        Temperature rate in K/fs.
        
    Raises
    ------
    ValueError
        If nsteps is zero.
        
    Examples
    --------
    >>> rate = estimate_temp_rate(300, 1500, 10000, 2.0)  # K/fs
    """
    if nsteps == 0:
        raise ValueError("nsteps must not be zero.")
    return (T_end - T_start) / (nsteps * timestep_fs)


def estimate_steps(T_start: float, T_end: float, temp_rate: float, timestep_fs: float = 1.0) -> int:
    """
    Estimate the number of steps required to go from T_start to T_end with a given temp_rate.
    
    Parameters
    ----------
    T_start : float
        Initial temperature in K.
    T_end : float
        Final temperature in K.
    temp_rate : float
        Temperature rate in K/fs.
    timestep_fs : float, optional
        Timestep in femtoseconds (default: 1.0).
        
    Returns
    -------
    int
        Number of steps required.
        
    Raises
    ------
    ValueError
        If temp_rate is zero.
        
    Examples
    --------
    >>> steps = estimate_steps(300, 1500, 0.12, 2.0)
    """
    if temp_rate == 0:
        raise ValueError("temp_rate must not be zero.")
    return int(abs((T_end - T_start) / (temp_rate * timestep_fs)))