"""
Adaptive sampling utilities for optimization and function fitting.

This module provides tools for adaptive sampling using Gaussian processes
and Bayesian optimization techniques, particularly useful for equation of
state fitting and parameter optimization.
"""

from typing import Callable, Tuple
import numpy as np
from skopt import Optimizer
from skopt.learning import GaussianProcessRegressor
import inspect


def adaptive_sampling_skopt(
    f: Callable,
    x_bounds: Tuple[float, float],
    n_initial: int = 5,
    n_iter: int = 10,
    num_jobs: int = 1,
    noise: float = 1e-6,
    random_state: int = None,
    min_dist: float = 0.05
) -> Tuple[np.ndarray, np.ndarray, GaussianProcessRegressor]:
    """
    Adaptive sampling using scikit-optimize for fitting analytical functions.
    
    This function uses Gaussian process regression with acquisition functions
    to intelligently sample points for function approximation, with batch 
    sampling and diversity constraints.
    
    Parameters
    ----------
    f : callable
        The target function to sample (e.g., experiment or simulation).
        Should take a single float or 1D array input.
    x_bounds : tuple of float
        (min, max) bounds for the input variable x.
    n_initial : int, optional
        Number of initial random samples (default: 5).
    n_iter : int, optional
        Number of adaptive sampling iterations (default: 10).
    num_jobs : int, optional
        Number of new points to sample per iteration (batch size) (default: 1).
    noise : float, optional
        Noise level for the GP regressor (default: 1e-6).
    random_state : int, optional
        Random seed for reproducibility (default: None).
    min_dist : float, optional
        Minimum allowed distance between new samples and all previous samples
        (default: 0.05).
        
    Returns
    -------
    tuple
        (x_samples, y_samples, gp) where:
        - x_samples: np.ndarray of sampled x values
        - y_samples: np.ndarray of corresponding y values  
        - gp: GaussianProcessRegressor for uncertainty estimation
        
    Examples
    --------
    >>> def expensive_function(x):
    ...     return x**2 + 0.1 * np.sin(10*x)  # Simulate expensive calculation
    >>> x_samples, y_samples, gp = adaptive_sampling_skopt(
    ...     expensive_function, 
    ...     x_bounds=(0, 1),
    ...     n_initial=3,
    ...     n_iter=5
    ... )
    >>> # Use GP for predictions
    >>> x_pred = np.linspace(0, 1, 100)
    >>> y_pred, y_std = gp.predict(x_pred.reshape(-1, 1), return_std=True)
    """
    bounds = [(x_bounds[0], x_bounds[1])]
    rng = np.random.default_rng(random_state)
    
    # Initial random sampling
    x_samples = rng.uniform(x_bounds[0], x_bounds[1], n_initial).reshape(-1, 1)
    y_samples = np.array([f(x[0]) for x in x_samples])
    
    # Set up Gaussian Process and optimizer
    gp = GaussianProcessRegressor(
        alpha=noise**2, 
        normalize_y=True, 
        random_state=random_state
    )
    opt = Optimizer(
        dimensions=bounds, 
        base_estimator=gp, 
        acq_func="LCB",  # Lower Confidence Bound
        random_state=random_state
    )
    
    # Tell optimizer about initial samples
    opt.tell(x_samples.tolist(), y_samples.tolist())
    
    # Adaptive sampling iterations
    for iteration in range(n_iter):
        # Ask for more candidates than needed to allow filtering
        candidate_points = opt.ask(n_points=num_jobs * 2)
        
        # Filter candidates to maintain minimum distance
        selected_points = []
        for candidate in candidate_points:
            # Check distance to all existing samples
            distances = np.abs(x_samples.flatten() - candidate[0])
            if np.all(distances > min_dist):
                selected_points.append(candidate)
            
            # Stop when we have enough points
            if len(selected_points) == num_jobs:
                break
        
        # If no candidates satisfy distance constraint, pick the farthest one
        if not selected_points:
            distances_min = [np.min(np.abs(x_samples.flatten() - candidate[0])) 
                           for candidate in candidate_points]
            best_idx = np.argmax(distances_min)
            selected_points = [candidate_points[best_idx]]
        
        # Evaluate function at selected points
        new_y_values = [f(point[0]) for point in selected_points]
        
        # Update optimizer with new data
        opt.tell(selected_points, new_y_values)
        
        # Add to sample arrays
        x_samples = np.vstack([x_samples, np.array(selected_points)])
        y_samples = np.append(y_samples, new_y_values)
    
    # Fit final GP model
    gp.fit(x_samples, y_samples)
    
    return x_samples.flatten(), y_samples, gp