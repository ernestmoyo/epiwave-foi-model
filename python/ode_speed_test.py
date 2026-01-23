"""
TensorFlow ODE Solver Speed Test for Ross-MacDonald Model
=========================================================
This script compares different ODE solvers in TensorFlow Probability
for the Ross-MacDonald malaria transmission model with:
- Time-varying parameters (bimodal seasonality)
- Multiple sites (N = 10, 50, 100, 500)
- 4 years of simulation

Author: Ernest Moyo
Date: 2025
Issue: https://github.com/ernestmoyo/epiwave-foi-model/issues/4
"""

import tensorflow as tf
import tensorflow_probability as tfp
import numpy as np
import time
import pandas as pd
from typing import Callable, Dict, List, Tuple

# Enable eager execution for debugging
tf.config.run_functions_eagerly(False)

# ============================================================================
# ROSS-MACDONALD MODEL DEFINITION
# ============================================================================

class RossMacDonaldODE:
    """
    Ross-MacDonald malaria transmission model with time-varying mosquito density.
    
    State variables:
        x: Proportion of infected humans
        z: Proportion of infected mosquitoes
    
    Parameters:
        m(t): Mosquito abundance relative to humans (time-varying, seasonal)
        a: Human biting rate
        b: Probability of transmission mosquito -> human
        c: Probability of transmission human -> mosquito
        r: Human recovery rate
        g: Mosquito death rate
    """
    
    def __init__(self, n_sites: int, params: Dict[str, tf.Tensor]):
        self.n_sites = n_sites
        self.params = params
        
    def seasonal_mosquito_density(self, t: tf.Tensor) -> tf.Tensor:
        """
        Bimodal seasonal mosquito density function.
        Two peaks per year (e.g., after rainy seasons).
        
        m(t) = m_base * (1 + amp1 * cos(2*pi*t/365) + amp2 * cos(4*pi*t/365))
        """
        m_base = self.params["m_base"]  # Shape: (n_sites,)
        amp1 = self.params["amp1"]       # Primary seasonal amplitude
        amp2 = self.params["amp2"]       # Secondary seasonal amplitude
        
        # Bimodal seasonality
        seasonal = (1.0 + 
                   amp1 * tf.cos(2.0 * np.pi * t / 365.0) + 
                   amp2 * tf.cos(4.0 * np.pi * t / 365.0))
        
        return m_base * seasonal
    
    def __call__(self, t: tf.Tensor, state: tf.Tensor) -> tf.Tensor:
        """
        ODE right-hand side function.
        
        Args:
            t: Current time (scalar)
            state: Current state [x, z] with shape (2, n_sites)
            
        Returns:
            Derivatives [dx/dt, dz/dt] with shape (2, n_sites)
        """
        x = state[0]  # Human infection proportion (n_sites,)
        z = state[1]  # Mosquito infection proportion (n_sites,)
        
        # Get time-varying mosquito density
        m = self.seasonal_mosquito_density(t)
        
        # Fixed parameters
        a = self.params["a"]  # Biting rate
        b = self.params["b"]  # Transmission prob mosquito -> human
        c = self.params["c"]  # Transmission prob human -> mosquito
        r = self.params["r"]  # Human recovery rate
        g = self.params["g"]  # Mosquito death rate
        
        # Ross-MacDonald equations
        dx_dt = m * a * b * z * (1.0 - x) - r * x
        dz_dt = a * c * x * (1.0 - z) - g * z
        
        return tf.stack([dx_dt, dz_dt], axis=0)


# ============================================================================
# SOLVER IMPLEMENTATIONS
# ============================================================================

def solve_with_dormand_prince(ode_fn: Callable, 
                               initial_state: tf.Tensor,
                               t_span: Tuple[float, float],
                               t_eval: tf.Tensor) -> tf.Tensor:
    """Solve ODE using Dormand-Prince (adaptive RK45) method."""
    solver = tfp.math.ode.DormandPrince(atol=1e-6, rtol=1e-6)
    results = solver.solve(ode_fn, t_span[0], initial_state, 
                          solution_times=t_eval)
    return results.states


def solve_with_bdf(ode_fn: Callable,
                   initial_state: tf.Tensor,
                   t_span: Tuple[float, float],
                   t_eval: tf.Tensor) -> tf.Tensor:
    """Solve ODE using BDF (Backward Differentiation Formula) method."""
    solver = tfp.math.ode.BDF(atol=1e-6, rtol=1e-6)
    results = solver.solve(ode_fn, t_span[0], initial_state,
                          solution_times=t_eval)
    return results.states


def solve_with_fixed_euler(ode_fn: Callable,
                           initial_state: tf.Tensor,
                           t_span: Tuple[float, float],
                           dt: float = 0.1) -> tf.Tensor:
    """Solve ODE using fixed-step Euler method."""
    t_start, t_end = t_span
    n_steps = int((t_end - t_start) / dt)
    
    states = [initial_state]
    state = initial_state
    t = t_start
    
    for _ in range(n_steps):
        derivative = ode_fn(t, state)
        state = state + dt * derivative
        t = t + dt
        states.append(state)
    
    return tf.stack(states, axis=0)


def solve_with_rk4(ode_fn: Callable,
                   initial_state: tf.Tensor,
                   t_span: Tuple[float, float],
                   dt: float = 0.5) -> tf.Tensor:
    """Solve ODE using fixed-step RK4 method."""
    t_start, t_end = t_span
    n_steps = int((t_end - t_start) / dt)
    
    states = [initial_state]
    state = initial_state
    t = t_start
    
    for _ in range(n_steps):
        k1 = ode_fn(t, state)
        k2 = ode_fn(t + dt/2, state + dt/2 * k1)
        k3 = ode_fn(t + dt/2, state + dt/2 * k2)
        k4 = ode_fn(t + dt, state + dt * k3)
        
        state = state + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)
        t = t + dt
        states.append(state)
    
    return tf.stack(states, axis=0)


# ============================================================================
# SPEED TEST FUNCTIONS
# ============================================================================

def create_test_parameters(n_sites: int) -> Dict[str, tf.Tensor]:
    """Create realistic test parameters for n_sites locations."""
    
    # Base mosquito density varies by site (log-normal variation)
    np.random.seed(42)
    m_base = tf.constant(np.random.lognormal(mean=1.0, sigma=0.5, size=n_sites), 
                        dtype=tf.float32)
    
    return {
        "m_base": m_base,
        "amp1": tf.constant(0.6, dtype=tf.float32),   # Primary seasonal amplitude
        "amp2": tf.constant(0.3, dtype=tf.float32),   # Secondary amplitude
        "a": tf.constant(0.3, dtype=tf.float32),      # Biting rate
        "b": tf.constant(0.5, dtype=tf.float32),      # Transmission mosquito->human
        "c": tf.constant(0.8, dtype=tf.float32),      # Transmission human->mosquito
        "r": tf.constant(1/14, dtype=tf.float32),     # Recovery rate (14 day duration)
        "g": tf.constant(1/10, dtype=tf.float32),     # Mosquito death rate (10 day lifespan)
    }


def run_speed_test(n_sites_list: List[int] = [10, 50, 100, 500],
                   n_years: int = 4,
                   n_repeats: int = 5) -> pd.DataFrame:
    """
    Run speed tests for different solvers and number of sites.
    
    Returns DataFrame with timing results.
    """
    
    results = []
    t_span = (0.0, 365.0 * n_years)
    t_eval = tf.linspace(0.0, 365.0 * n_years, 365 * n_years)
    
    solvers = {
        "DormandPrince": solve_with_dormand_prince,
        "BDF": solve_with_bdf,
        "RK4_fixed": solve_with_rk4,
        "Euler_fixed": solve_with_fixed_euler,
    }
    
    for n_sites in n_sites_list:
        print(f"\nTesting with N = {n_sites} sites...")
        
        # Create model
        params = create_test_parameters(n_sites)
        model = RossMacDonaldODE(n_sites, params)
        
        # Initial state: low infection levels
        initial_state = tf.constant(
            [[0.01] * n_sites,  # x: human infection
             [0.0] * n_sites],  # z: mosquito infection
            dtype=tf.float32
        )
        
        for solver_name, solver_fn in solvers.items():
            print(f"  Testing {solver_name}...")
            
            times = []
            for rep in range(n_repeats):
                try:
                    start = time.perf_counter()
                    
                    if "fixed" in solver_name.lower():
                        _ = solver_fn(model, initial_state, t_span)
                    else:
                        _ = solver_fn(model, initial_state, t_span, t_eval)
                    
                    end = time.perf_counter()
                    times.append(end - start)
                    
                except Exception as e:
                    print(f"    Error: {e}")
                    times.append(np.nan)
            
            results.append({
                "n_sites": n_sites,
                "solver": solver_name,
                "mean_time": np.nanmean(times),
                "std_time": np.nanstd(times),
                "min_time": np.nanmin(times),
                "max_time": np.nanmax(times),
            })
            
            print(f"    Mean time: {np.nanmean(times):.4f}s")
    
    return pd.DataFrame(results)


def run_gradient_speed_test(n_sites_list: List[int] = [10, 50, 100],
                            n_years: int = 1,
                            n_repeats: int = 3) -> pd.DataFrame:
    """
    Speed test differentiation through ODE solvers.
    This is critical for Bayesian inference with greta/TensorFlow.
    """
    
    results = []
    t_span = (0.0, 365.0 * n_years)
    
    for n_sites in n_sites_list:
        print(f"\nGradient test with N = {n_sites} sites...")
        
        # Create model with trainable parameters
        m_base = tf.Variable(np.random.lognormal(1.0, 0.5, n_sites).astype(np.float32))
        
        params = {
            "m_base": m_base,
            "amp1": tf.constant(0.6, dtype=tf.float32),
            "amp2": tf.constant(0.3, dtype=tf.float32),
            "a": tf.constant(0.3, dtype=tf.float32),
            "b": tf.constant(0.5, dtype=tf.float32),
            "c": tf.constant(0.8, dtype=tf.float32),
            "r": tf.constant(1/14, dtype=tf.float32),
            "g": tf.constant(1/10, dtype=tf.float32),
        }
        
        initial_state = tf.constant(
            [[0.01] * n_sites, [0.0] * n_sites],
            dtype=tf.float32
        )
        
        t_eval = tf.linspace(0.0, 365.0 * n_years, 100)
        
        # Test gradient computation with DormandPrince
        print("  Testing gradient with DormandPrince...")
        times = []
        
        for rep in range(n_repeats):
            try:
                start = time.perf_counter()
                
                with tf.GradientTape() as tape:
                    model = RossMacDonaldODE(n_sites, params)
                    solver = tfp.math.ode.DormandPrince(atol=1e-5, rtol=1e-5)
                    solution = solver.solve(model, t_span[0], initial_state,
                                           solution_times=t_eval)
                    # Loss: sum of final human infection proportions
                    loss = tf.reduce_sum(solution.states[-1, 0, :])
                
                grads = tape.gradient(loss, m_base)
                
                end = time.perf_counter()
                times.append(end - start)
                
            except Exception as e:
                print(f"    Error: {e}")
                times.append(np.nan)
        
        results.append({
            "n_sites": n_sites,
            "test_type": "gradient",
            "mean_time": np.nanmean(times),
            "std_time": np.nanstd(times),
        })
        
        print(f"    Mean gradient time: {np.nanmean(times):.4f}s")
    
    return pd.DataFrame(results)


# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("TensorFlow ODE Solver Speed Test")
    print("Ross-MacDonald Model with Time-Varying Parameters")
    print("=" * 60)
    
    # Check TensorFlow version
    print(f"\nTensorFlow version: {tf.__version__}")
    print(f"TensorFlow Probability version: {tfp.__version__}")
    
    # Run forward solve speed tests
    print("\n" + "=" * 60)
    print("FORWARD SOLVE SPEED TEST")
    print("=" * 60)
    forward_results = run_speed_test(
        n_sites_list=[10, 50, 100, 500],
        n_years=4,
        n_repeats=5
    )
    
    print("\n" + "-" * 40)
    print("Forward Solve Results Summary:")
    print("-" * 40)
    print(forward_results.to_string(index=False))
    
    # Run gradient speed tests
    print("\n" + "=" * 60)
    print("GRADIENT COMPUTATION SPEED TEST")
    print("=" * 60)
    gradient_results = run_gradient_speed_test(
        n_sites_list=[10, 50, 100],
        n_years=1,
        n_repeats=3
    )
    
    print("\n" + "-" * 40)
    print("Gradient Results Summary:")
    print("-" * 40)
    print(gradient_results.to_string(index=False))
    
    # Save results
    forward_results.to_csv("forward_solve_results.csv", index=False)
    gradient_results.to_csv("gradient_results.csv", index=False)
    print("\nResults saved to CSV files.")

