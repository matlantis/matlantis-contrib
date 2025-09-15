"""
Input/output utilities for MD simulations.

This module provides functions for reading configuration files and
handling simulation data I/O operations.
"""

from typing import Dict, Any
from pathlib import Path
import yaml


def read_yaml(yamlfile: str = "md.yaml") -> Dict[str, Any]:
    """
    Read YAML configuration file for MD simulation parameters.
    
    This function loads simulation parameters from a YAML file, with
    graceful handling of missing files.
    
    Parameters
    ----------
    yamlfile : str, optional
        Path to the YAML configuration file (default: "md.yaml").
        
    Returns
    -------
    dict
        Dictionary containing configuration parameters. Returns empty
        dict if file doesn't exist.
        
    Examples
    --------
    >>> config = read_yaml("simulation_config.yaml")
    >>> temperature = config.get("temperature", 298.15)
    
    Example YAML file content:
    ```yaml
    # MD simulation parameters
    num_steps: 10000
    timestep: 2.0
    temperature: 1000
    temperatures: [800, 1000, 1200, 1400, 1600, 1800]
    pressure: 1.0
    ensemble: "npt"
    ```
    """
    yamlfile_path = Path(yamlfile)
    
    if not yamlfile_path.exists():
        return {}
    
    try:
        with open(yamlfile_path, 'r') as f:
            config = yaml.safe_load(f)
        return config if config is not None else {}
    except yaml.YAMLError as e:
        print(f"Error reading YAML file {yamlfile}: {e}")
        return {}
    except Exception as e:
        print(f"Unexpected error reading {yamlfile}: {e}")
        return {}