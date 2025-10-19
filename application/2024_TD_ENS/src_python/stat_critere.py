import numpy as np
import pandas as pd
from collections import namedtuple
from dataclasses import dataclass
from random import uniform, gauss
from numpy import inf, nansum, log, size, var, mean, isclose, sqrt, zeros, all
from typing import Callable, Tuple, Sequence
from scipy import stats
from Direct_model import setup_ginette, generate_zone_parameters, initial_conditions, boundary_conditions, run_direct_model, remove_first_two_days

# =============================================================================
# STATISTICAL CRITERIA FUNCTIONS
# =============================================================================

def rmse(predictions, ref, sigma=1.0):
    """
    Calculate Root Mean Squared Error with optional normalization by sigma.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values  
    - sigma: normalization factor (default=1.0)
    
    Returns:
    - RMSE value
    """
    differences = ((predictions - ref) / sigma)
    differences_squared = differences ** 2
    mean_of_differences_squared = differences_squared.mean()
    rmse_val = np.sqrt(mean_of_differences_squared)
    return rmse_val

def mae(predictions, ref):
    """
    Calculate Mean Absolute Error.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - MAE value
    """
    return np.mean(np.abs(predictions - ref))

def bias(predictions, ref):
    """
    Calculate Bias (mean error).
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - Bias value
    """
    return np.mean(predictions - ref)

def pbias(predictions, ref):
    """
    Calculate Percent Bias.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - PBIAS in percentage
    """
    if np.sum(ref) == 0:
        return np.inf
    return 100 * np.sum(predictions - ref) / np.sum(ref)

def nse(predictions, ref):
    """
    Calculate Nash-Sutcliffe Efficiency.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - NSE value
    """
    ss_res = np.sum((ref - predictions) ** 2)
    ss_tot = np.sum((ref - np.mean(ref)) ** 2)
    if ss_tot == 0:
        return -np.inf
    return 1 - (ss_res / ss_tot)

def kge(predictions, ref):
    """
    Calculate Kling-Gupta Efficiency.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - KGE value
    """
    r = np.corrcoef(predictions, ref)[0, 1]
    alpha = np.std(predictions) / np.std(ref) if np.std(ref) != 0 else np.inf
    beta = np.mean(predictions) / np.mean(ref) if np.mean(ref) != 0 else np.inf
    
    kge_val = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (beta - 1)**2)
    return kge_val

def d_index(predictions, ref):
    """
    Calculate Index of Agreement (Willmott, 1981).
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - Index of Agreement value
    """
    numerator = np.sum((ref - predictions)**2)
    ref_mean = np.mean(ref)
    denominator = np.sum((np.abs(predictions - ref_mean) + np.abs(ref - ref_mean))**2)
    
    if denominator == 0:
        return 0
    return 1 - (numerator / denominator)

def modified_nse(predictions, ref, j=1):
    """
    Calculate Modified Nash-Sutcliffe Efficiency.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    - j: exponent for modification (default=1)
    
    Returns:
    - Modified NSE value
    """
    numerator = np.sum(np.abs(ref - predictions)**j)
    denominator = np.sum(np.abs(ref - np.mean(ref))**j)
    
    if denominator == 0:
        return -np.inf
    return 1 - (numerator / denominator)

def log_nse(predictions, ref):
    """
    Calculate Log Nash-Sutcliffe Efficiency.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - Log NSE value
    """
    # Add small constant to avoid log(0)
    epsilon = 1e-10
    log_pred = np.log(predictions + epsilon)
    log_ref = np.log(ref + epsilon)
    
    return nse(log_pred, log_ref)

def correlation_coefficient(predictions, ref):
    """
    Calculate Pearson correlation coefficient.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - Correlation coefficient and p-value
    """
    return stats.pearsonr(predictions, ref)

def volumetric_efficiency(predictions, ref):
    """
    Calculate Volumetric Efficiency.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - VE value
    """
    numerator = np.sum(np.abs(predictions - ref))
    denominator = np.sum(ref)
    
    if denominator == 0:
        return -np.inf
    return 1 - (numerator / denominator)

def calculate_all_metrics(predictions, ref):
    """
    Calculate all available statistical metrics.
    
    Parameters:
    - predictions: simulated/predicted values
    - ref: reference/observed values
    
    Returns:
    - Dictionary containing all metrics
    """
    # Remove NaN values
    mask = ~(np.isnan(predictions) | np.isnan(ref))
    pred_clean = predictions[mask]
    ref_clean = ref[mask]
    
    if len(pred_clean) < 2:
        return None
    
    # Calculate correlation
    r, p_value = correlation_coefficient(pred_clean, ref_clean)
    
    metrics = {
        'RMSE': rmse(pred_clean, ref_clean),
        'MAE': mae(pred_clean, ref_clean),
        'Bias': bias(pred_clean, ref_clean),
        'PBIAS': pbias(pred_clean, ref_clean),
        'NSE': nse(pred_clean, ref_clean),
        'KGE': kge(pred_clean, ref_clean),
        'D_index': d_index(pred_clean, ref_clean),
        'Modified_NSE': modified_nse(pred_clean, ref_clean),
        'Log_NSE': log_nse(pred_clean, ref_clean),
        'VE': volumetric_efficiency(pred_clean, ref_clean),
        'R': r,
        'R_squared': r**2,
        'p_value': p_value,
        'n_points': len(pred_clean)
    }
    
    return metrics

def calculate_metrics(data):
    """
    Calculate metrics for temperature sensors in merged dataset.
    
    Parameters:
    - data: DataFrame with columns like 'Temp1_sim', 'Temp1_obs', etc.
    
    Returns:
    - Dictionary of metrics for each sensor
    """
    metrics = {}

    for i in range(1, 5):  # Check Temp1 through Temp4
        sensor = f'Temp{i}'
        sim_col = f'{sensor}_sim'
        obs_col = f'{sensor}_obs'
        
        if sim_col in data.columns and obs_col in data.columns:
            sim_data = data[sim_col].dropna()
            obs_data = data[obs_col].dropna()
            
            # Ensure same length
            min_len = min(len(sim_data), len(obs_data))
            if min_len < 10:  # Need minimum data points
                continue
                
            sim_data = sim_data.iloc[:min_len]
            obs_data = obs_data.iloc[:min_len]
            
            # Calculate all metrics
            sensor_metrics = calculate_all_metrics(sim_data.values, obs_data.values)
            if sensor_metrics is not None:
                metrics[sensor] = sensor_metrics

    return metrics

def categorize_performance(nse_val, r2_val, kge_val=None):
    """
    Categorize model performance based on statistical criteria.
    
    Parameters:
    - nse_val: Nash-Sutcliffe Efficiency value
    - r2_val: R-squared value
    - kge_val: Kling-Gupta Efficiency value (optional)
    
    Returns:
    - Performance category string
    """
    # Use KGE if provided, otherwise use NSE and R2
    if kge_val is not None:
        if nse_val > 0.75 and r2_val > 0.75 and kge_val > 0.75:
            return "Excellent"
        elif nse_val > 0.65 and r2_val > 0.65 and kge_val > 0.65:
            return "Very Good"
        elif nse_val > 0.50 and r2_val > 0.50 and kge_val > 0.50:
            return "Good"
        elif nse_val > 0.40 and r2_val > 0.40 and kge_val > 0.40:
            return "Satisfactory"
        else:
            return "Unsatisfactory"
    else:
        if nse_val > 0.75 and r2_val > 0.75:
            return "Excellent"
        elif nse_val > 0.65 and r2_val > 0.65:
            return "Very Good"
        elif nse_val > 0.50 and r2_val > 0.50:
            return "Good"
        elif nse_val > 0.40 and r2_val > 0.40:
            return "Satisfactory"
        else:
            return "Unsatisfactory"

def get_performance_thresholds():
    """
    Return standard performance thresholds for hydrological model evaluation.
    
    Returns:
    - Dictionary with performance criteria and their thresholds
    """
    return {
        'NSE': {
            'Excellent': '>0.75',
            'Very Good': '0.65-0.75',
            'Good': '0.50-0.65',
            'Satisfactory': '0.40-0.50',
            'Unsatisfactory': '<0.40'
        },
        'R²': {
            'Excellent': '>0.75',
            'Very Good': '0.65-0.75',
            'Good': '0.50-0.65',
            'Satisfactory': '0.40-0.50',
            'Unsatisfactory': '<0.40'
        },
        'KGE': {
            'Excellent': '>0.75',
            'Very Good': '0.65-0.75',
            'Good': '0.50-0.65',
            'Satisfactory': '0.40-0.50',
            'Unsatisfactory': '<0.40'
        },
        'PBIAS': {
            'Excellent': '<±10%',
            'Very Good': '±10-15%',
            'Good': '±15-25%',
            'Satisfactory': '±25-40%',
            'Unsatisfactory': '>±40%'
        }
    }

def evaluate_model_performance(metrics_dict):
    """
    Evaluate overall model performance from a dictionary of metrics.
    
    Parameters:
    - metrics_dict: Dictionary containing 'NSE', 'R_squared', 'KGE', etc.
    
    Returns:
    - Dictionary with performance category and individual assessments
    """
    nse_val = metrics_dict.get('NSE', -np.inf)
    r2_val = metrics_dict.get('R_squared', 0)
    kge_val = metrics_dict.get('KGE', None)
    pbias_val = abs(metrics_dict.get('PBIAS', 100))
    
    # Overall category
    overall_category = categorize_performance(nse_val, r2_val, kge_val)
    
    # Individual assessments
    individual_assessments = {
        'NSE': categorize_performance(nse_val, 1.0),  # Use dummy R2=1.0 for NSE-only assessment
        'R²': categorize_performance(1.0, r2_val),   # Use dummy NSE=1.0 for R2-only assessment
    }
    
    if kge_val is not None:
        individual_assessments['KGE'] = categorize_performance(1.0, 1.0, kge_val)
    
    # PBIAS assessment
    if pbias_val < 10:
        individual_assessments['PBIAS'] = "Excellent"
    elif pbias_val < 15:
        individual_assessments['PBIAS'] = "Very Good"
    elif pbias_val < 25:
        individual_assessments['PBIAS'] = "Good"
    elif pbias_val < 40:
        individual_assessments['PBIAS'] = "Satisfactory"
    else:
        individual_assessments['PBIAS'] = "Unsatisfactory"
    
    return {
        'overall_category': overall_category,
        'individual_assessments': individual_assessments,
        'metrics_values': {
            'NSE': nse_val,
            'R²': r2_val,
            'KGE': kge_val,
            'PBIAS': pbias_val
        }
    }

# --- Parameter Classes ---
# --- Model Parameters, Priors, and Perturbation ---
class Param:
    """Stores parameters for a layer in a named tuple-like format."""
    def __init__(self, REF_k, REF_n, REF_l, REF_r):
        self.REF_k = REF_k
        self.REF_n = REF_n
        self.REF_l = REF_l
        self.REF_r = REF_r

    def __repr__(self):
        return f"Param(REF_k={self.REF_k}, REF_n={self.REF_n}, REF_l={self.REF_l}, REF_r={self.REF_r})"


class Prior:
    """Defines a prior distribution for a model parameter."""
    def __init__(self, range: Tuple[float, float], sigma: float, density: Callable[[float], float] = lambda x: 1.0):
        self.range = range
        self.sigma = sigma
        self.density = density

    def perturb(self, val: float) -> float:
        """Perturbs a parameter value, adding random noise while staying within bounds."""
        new_val = val + gauss(0, self.sigma)
        while new_val > self.range[1]:
            new_val = val + gauss(0, self.sigma)
        while new_val < self.range[0]:
            new_val = val + gauss(0, self.sigma)
        return new_val

# --- Model Execution ---
@dataclass
class SimulationResult:
    """Stores the result of a model simulation."""
    time: np.ndarray
    output: np.ndarray
    params: Param
    success: bool
    message: str

def run_simulation_with_params(param_values: Tuple[float, float, float, float], 
                                time: np.ndarray, 
                                forcing: np.ndarray, 
                                obs: np.ndarray, 
                                burn_in: int = 0) -> SimulationResult:
    """
    Run the model simulation with given parameter values.
    
    Parameters:
    - param_values: Tuple of parameter values (k, n, l, r)
    - time: Time array for the simulation
    - forcing: Forcing data (e.g., precipitation, temperature)
    - obs: Observed data for comparison
    - burn_in: Number of initial timesteps to discard (default=0)
    
    Returns:
    - SimulationResult object containing time, output, parameters, success flag, and message
    """
    # Unpack parameter values
    k, n, l, r = param_values
    
    # Create Param object
    params = Param(k, n, l, r)
    
    # Run the model (replace with actual model function)
    try:
        output = run_direct_model(params, time, forcing)
        
        # Apply burn-in removal if needed
        if burn_in > 0:
            output = output[burn_in:]
            time = time[burn_in:]
        
        return SimulationResult(time=time, output=output, params=params, success=True, message="Simulation run successful.")
    except Exception as e:
        return SimulationResult(time=np.array([]), output=np.array([]), params=params, success=False, message=str(e))

def run_multiple_simulations(param_grid: pd.DataFrame, 
                             time: np.ndarray, 
                             forcing: np.ndarray, 
                             obs: np.ndarray, 
                             burn_in: int = 0) -> pd.DataFrame:
    """
    Run multiple simulations over a grid of parameters.
    
    Parameters:
    - param_grid: DataFrame where each row is a set of parameters (k, n, l, r)
    - time: Time array for the simulation
    - forcing: Forcing data (e.g., precipitation, temperature)
    - obs: Observed data for comparison
    - burn_in: Number of initial timesteps to discard (default=0)
    
    Returns:
    - DataFrame with simulation results, one row per simulation
    """
    results = []

    for i, row in param_grid.iterrows():
        # Extract parameter values
        param_values = (row['k'], row['n'], row['l'], row['r'])
        
        # Run simulation
        result = run_simulation_with_params(param_values, time, forcing, obs, burn_in)
        
        # Store result
        results.append({
            'params': result.params,
            'time': result.time,
            'output': result.output,
            'success': result.success,
            'message': result.message,
            'nse': nse(result.output, obs),
            'r_squared': correlation_coefficient(result.output, obs)[0]**2,
            'kge': kge(result.output, obs),
            'pbias': pbias(result.output, obs)
        })

    return pd.DataFrame(results)

def summarize_simulation_results(results: pd.DataFrame, 
                                 obs: np.ndarray, 
                                 metrics: list = None) -> pd.DataFrame:
    """
    Summarize the results of simulations, calculating statistics and metrics.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metrics: List of metrics to calculate (default=None calculates all)
    
    Returns:
    - DataFrame with summary statistics and metrics for each simulation
    """
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    summary = {}

    for metric in metrics:
        summary[metric] = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    
    # Convert to DataFrame
    summary_df = pd.DataFrame(summary)
    
    # Add parameter columns
    for param in ['k', 'n', 'l', 'r']:
        summary_df[param] = results['params'].apply(lambda p: getattr(p, f'REF_{param}'))
    
    return summary_df

def plot_simulation_results(results: pd.DataFrame, obs: np.ndarray, metric: str = 'nse'):
    """
    Plot the results of simulations, showing the distribution of a metric.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metric: Metric to plot (default='nse')
    """
    import matplotlib.pyplot as plt
    
    # Extract metric values
    metric_values = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    
    plt.figure(figsize=(10, 6))
    plt.hist(metric_values, bins=20, color='skyblue', edgecolor='black')
    plt.axvline(np.nanmean(metric_values), color='red', linestyle='dashed', linewidth=1)
    plt.title(f'Distribution of {metric.upper()}')
    plt.xlabel(metric.upper())
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    plt.show()

def plot_comparison_with_observed(results: pd.DataFrame, obs: np.ndarray, index: int = 0):
    """
    Plot the simulated vs observed values for a specific simulation result.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - index: Index of the simulation result to plot (default=0)
    """
    import matplotlib.pyplot as plt
    
    # Extract time and output for the selected result
    time = results.loc[index, 'time']
    output = results.loc[index, 'output']
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, output, label='Simulated', color='skyblue')
    plt.plot(time, obs, label='Observed', color='red', linestyle='dashed')
    plt.title('Simulated vs Observed Values')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_distributions(results: pd.DataFrame):
    """
    Plot the distributions of parameters from the simulation results.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    params = ['k', 'n', 'l', 'r']
    
    plt.figure(figsize=(12, 8))
    
    for i, param in enumerate(params, 1):
        plt.subplot(2, 2, i)
        plt.hist(results[param], bins=20, color='skyblue', edgecolor='black')
        plt.title(f'Distribution of Parameter {param.upper()}')
        plt.xlabel(param.upper())
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_correlation_matrix(results: pd.DataFrame):
    """
    Plot the correlation matrix of parameters and metrics.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    # Select numeric columns for correlation
    numeric_cols = results.select_dtypes(include=[np.number]).columns.tolist()
    
    # Compute the correlation matrix
    corr = results[numeric_cols].corr()
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap='coolwarm')
    plt.title('Correlation Matrix')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

def plot_metric_trends_over_time(results: pd.DataFrame, obs: np.ndarray, metric: str = 'nse'):
    """
    Plot the trends of a metric over time for the best and worst simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metric: Metric to plot (default='nse')
    """
    import matplotlib.pyplot as plt
    
    # Extract time, metric values, and indices of best/worst simulations
    time = results['time'].iloc[0]
    metric_values = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    best_index = metric_values.idxmax()
    worst_index = metric_values.idxmin()
    
    # Extract best and worst simulation data
    best_simulation = results.loc[best_index]
    worst_simulation = results.loc[worst_index]
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, best_simulation['output'], label='Best Simulation', color='green')
    plt.plot(time, worst_simulation['output'], label='Worst Simulation', color='red')
    plt.plot(time, obs, label='Observed', color='blue', linestyle='dashed')
    plt.title(f'Trend of {metric.upper()} Over Time')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_sensitivity(results: pd.DataFrame, obs: np.ndarray, parameter: str = 'k'):
    """
    Plot the sensitivity of model performance to changes in a specific parameter.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - parameter: Parameter to analyze (default='k')
    """
    import matplotlib.pyplot as plt
    
    # Extract the parameter values and corresponding NSE
    param_values = results[parameter]
    nse_values = results['nse']
    
    plt.figure(figsize=(10, 6))
    plt.scatter(param_values, nse_values, alpha=0.7, color='skyblue')
    plt.title(f'Sensitivity of NSE to {parameter.upper()}')
    plt.xlabel(parameter.upper())
    plt.ylabel('NSE')
    plt.grid()
    plt.show()

def plot_best_worst_simulations(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot the best and worst simulations based on NSE, including observed data.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Extract the best and worst results
    best_result = results.loc[results['nse'].idxmax()]
    worst_result = results.loc[results['nse'].idxmin()]
    
    # Time array (assuming same for all simulations)
    time = best_result['time']
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, best_result['output'], label='Best Simulation (NSE)', color='green')
    plt.plot(time, worst_result['output'], label='Worst Simulation (NSE)', color='red')
    plt.plot(time, obs, label='Observed', color='blue', linestyle='dashed')
    plt.title('Best and Worst Simulations Based on NSE')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_metric_boxplots(results: pd.DataFrame, metrics: list = None):
    """
    Plot boxplots of metrics for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to plot (default=None plots all)
    """
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        plt.boxplot(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                            patch_artist=True, 
                                            boxprops=dict(facecolor='skyblue', color='black'),
                                            medianprops=dict(color='red'))
        plt.title(f'Boxplot of {metric.upper()}')
        plt.ylabel(metric.upper())
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_simulation_success_rate(results: pd.DataFrame):
    """
    Plot the success rate of simulations over time.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    # Calculate success rate (percentage of successful simulations)
    success_rate = results['success'].value_counts(normalize=True).sort_index()
    
    plt.figure(figsize=(8, 6))
    plt.bar(success_rate.index.astype(str), success_rate.values, color='skyblue', edgecolor='black')
    plt.title('Simulation Success Rate')
    plt.xlabel('Success (1=True, 0=False)')
    plt.ylabel('Proportion')
    plt.grid(axis='y', alpha=0.75)
    plt.xticks(rotation=0)
    plt.show()

def plot_metric_correlations(results: pd.DataFrame, metrics: list = None):
    """
    Plot scatter plots and correlation coefficients between pairs of metrics.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to analyze (default=None analyzes all)
    """
    import matplotlib.pyplot as plt
    from itertools import combinations
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    # Generate all pairs of metrics
    metric_pairs = list(combinations(metrics, 2))
    
    n_cols = 2
    n_rows = (len(metric_pairs) + n_cols - 1) // n_cols  # Ceiling division
    
    plt.figure(figsize=(12, 6 * n_rows))
    
    for i, (metric1, metric2) in enumerate(metric_pairs, 1):
        plt.subplot(n_rows, n_cols, i)
        plt.scatter(results[metric1], results[metric2], alpha=0.7, color='skyblue')
        plt.title(f'{metric1.upper()} vs {metric2.upper()}')
        plt.xlabel(metric1.upper())
        plt.ylabel(metric2.upper())
        
        # Calculate and display correlation coefficient
        corr_coef = results[[metric1, metric2]].corr().iloc[0, 1]
        plt.text(0.05, 0.95, f'ρ = {corr_coef:.2f}', transform=plt.gca().transAxes,
                 fontsize=12, verticalalignment='top')
        
        plt.grid()
    
    plt.tight_layout()
    plt.show()

def plot_best_worst_parameters(results: pd.DataFrame):
    """
    Plot the parameter values for the best and worst simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    # Extract the best and worst results
    best_result = results.loc[results['nse'].idxmax()]
    worst_result = results.loc[results['nse'].idxmin()]
    
    params = ['k', 'n', 'l', 'r']
    x = np.arange(len(params))
    
    plt.figure(figsize=(10, 6))
    plt.bar(x - 0.2, best_result[params], 0.4, label='Best Simulation', color='green')
    plt.bar(x + 0.2, worst_result[params], 0.4, label='Worst Simulation', color='red')
    plt.xticks(x, params)
    plt.title('Parameter Values for Best and Worst Simulations')
    plt.ylabel('Parameter Value')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    plt.show()

def plot_metric_violinplots(results: pd.DataFrame, metrics: list = None):
    """
    Plot violin plots of metrics for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to plot (default=None plots all)
    """
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        plt.violinplot(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                              showmeans=True, 
                                              showmedians=True)
        plt.title(f'Violin Plot of {metric.upper()}')
        plt.ylabel(metric.upper())
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_simulation_time_series(results: pd.DataFrame, obs: np.ndarray, index: int = 0):
    """
    Plot the time series of simulated values for a specific simulation result.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - index: Index of the simulation result to plot (default=0)
    """
    import matplotlib.pyplot as plt
    
    # Extract time and output for the selected result
    time = results.loc[index, 'time']
    output = results.loc[index, 'output']
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, output, label='Simulated', color='skyblue')
    plt.plot(time, obs, label='Observed', color='red', linestyle='dashed')
    plt.title('Simulated Time Series')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_correlation_scatter(results: pd.DataFrame, parameter: str = 'k'):
    """
    Scatter plot of parameter values against NSE, colored by simulation success.
    
    Parameters:
    - results: DataFrame with simulation results
    - parameter: Parameter to plot (default='k')
    """
    import matplotlib.pyplot as plt
    
    # Extract the parameter values and corresponding NSE
    param_values = results[parameter]
    nse_values = results['nse']
    success_values = results['success']
    
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(param_values, nse_values, c=success_values, cmap='coolwarm', alpha=0.7)
    plt.title(f'Scatter Plot of NSE vs {parameter.upper()}')
    plt.xlabel(parameter.upper())
    plt.ylabel('NSE')
    plt.colorbar(scatter, label='Success (1=True, 0=False)')
    plt.grid()
    plt.show()

def plot_metric_histograms(results: pd.DataFrame, metrics: list = None):
    """
    Plot histograms of metrics for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to plot (default=None plots all)
    """
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        plt.hist(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                        bins=20, 
                                        color='skyblue', 
                                        edgecolor='black')
        plt.title(f'Histogram of {metric.upper()}')
        plt.xlabel(metric.upper())
        plt.ylabel('Frequency')
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_best_worst_metric_trends(results: pd.DataFrame, obs: np.ndarray, metric: str = 'nse'):
    """
    Plot the trends of a metric over time for the best and worst simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metric: Metric to plot (default='nse')
    """
    import matplotlib.pyplot as plt
    
    # Extract time, metric values, and indices of best/worst simulations
    time = results['time'].iloc[0]
    metric_values = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    best_index = metric_values.idxmax()
    worst_index = metric_values.idxmin()
    
    # Extract best and worst simulation data
    best_simulation = results.loc[best_index]
    worst_simulation = results.loc[worst_index]
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, best_simulation['output'], label='Best Simulation', color='green')
    plt.plot(time, worst_simulation['output'], label='Worst Simulation', color='red')
    plt.plot(time, obs, label='Observed', color='blue', linestyle='dashed')
    plt.title(f'Trend of {metric.upper()} Over Time')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_importance(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot the importance of parameters based on their sensitivity.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Calculate absolute sensitivity (partial derivatives) for each parameter
    param_sensitivities = results[['k', 'n', 'l', 'r']].apply(lambda row: np.abs(np.gradient(row)), axis=0)
    
    # Calculate mean sensitivity for each parameter
    mean_sensitivities = param_sensitivities.mean().sort_values(ascending=False)
    
    plt.figure(figsize=(8, 6))
    plt.bar(mean_sensitivities.index, mean_sensitivities.values, color='skyblue', edgecolor='black')
    plt.title('Parameter Importance Based on Sensitivity')
    plt.xlabel('Parameter')
    plt.ylabel('Mean Absolute Sensitivity')
    plt.grid(axis='y', alpha=0.75)
    plt.xticks(rotation=0)
    plt.show()

def plot_simulation_residuals(results: pd.DataFrame, obs: np.ndarray, index: int = 0):
    """
    Plot the residuals (differences between simulated and observed values) for a specific simulation result.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - index: Index of the simulation result to plot (default=0)
    """
    import matplotlib.pyplot as plt
    
    # Extract time, output, and observed data for the selected result
    time = results.loc[index, 'time']
    output = results.loc[index, 'output']
    
    # Calculate residuals
    residuals = output - obs
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, residuals, label='Residuals', color='skyblue')
    plt.axhline(0, color='red', linestyle='dashed', linewidth=1)
    plt.title('Residuals of Simulation')
    plt.xlabel('Time')
    plt.ylabel('Residual')
    plt.legend()
    plt.grid()
    plt.show()

def plot_metric_scatter_matrix(results: pd.DataFrame, metrics: list = None):
    """
    Plot a scatter matrix of metrics for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to include (default=None includes all)
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    # Select only the specified metrics
    subset = results[metrics]
    
    # Plot scatter matrix
    sns.pairplot(subset, diag_kind='kde', markers='o', palette='coolwarm')
    plt.suptitle('Scatter Matrix of Metrics', y=1.02)
    plt.show()

def plot_best_worst_simulation_outputs(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot the output of the best and worst simulations based on NSE, including observed data.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Extract the best and worst results
    best_result = results.loc[results['nse'].idxmax()]
    worst_result = results.loc[results['nse'].idxmin()]
    
    # Time array (assuming same for all simulations)
    time = best_result['time']
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, best_result['output'], label='Best Simulation (NSE)', color='green')
    plt.plot(time, worst_result['output'], label='Worst Simulation (NSE)', color='red')
    plt.plot(time, obs, label='Observed', color='blue', linestyle='dashed')
    plt.title('Best and Worst Simulation Outputs')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_pairwise_relationships(results: pd.DataFrame, parameters: list = None):
    """
    Plot pairwise relationships between parameters using scatter plots.
    
    Parameters:
    - results: DataFrame with simulation results
    - parameters: List of parameters to include (default=None includes all)
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if parameters is None:
        parameters = ['k', 'n', 'l', 'r']
    
    # Select only the specified parameters
    subset = results[parameters]
    
    # Plot pairwise relationships
    sns.pairplot(subset, diag_kind='kde', markers='o', palette='coolwarm')
    plt.suptitle('Pairwise Relationships Between Parameters', y=1.02)
    plt.show()

def plot_metric_trends_across_simulations(results: pd.DataFrame, obs: np.ndarray, metric: str = 'nse'):
    """
    Plot the trends of a metric across all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metric: Metric to plot (default='nse')
    """
    import matplotlib.pyplot as plt
    
    # Extract time and metric values
    time = results['time'].iloc[0]
    metric_values = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, metric_values, label='Metric Trend', color='skyblue')
    plt.title(f'Trend of {metric.upper()} Across Simulations')
    plt.xlabel('Time')
    plt.ylabel(metric.upper())
    plt.legend()
    plt.grid()
    plt.show()

def plot_simulation_success_overview(results: pd.DataFrame):
    """
    Plot an overview of simulation success rates and performance categories.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    # Calculate success rate (percentage of successful simulations)
    success_rate = results['success'].value_counts(normalize=True).sort_index()
    
    # Performance categorization
    performance_categories = ['Excellent', 'Very Good', 'Good', 'Satisfactory', 'Unsatisfactory']
    category_counts = results['nse'].apply(lambda x: categorize_performance(x, 1.0)).value_counts().reindex(performance_categories, fill_value=0)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    # Success rate plot
    axs[0].bar(success_rate.index.astype(str), success_rate.values, color='skyblue', edgecolor='black')
    axs[0].set_title('Simulation Success Rate')
    axs[0].set_xlabel('Success (1=True, 0=False)')
    axs[0].set_ylabel('Proportion')
    axs[0].grid(axis='y', alpha=0.75)
    axs[0].set_xticklabels(success_rate.index.astype(str), rotation=0)
    
    # Performance category plot
    axs[1].bar(category_counts.index, category_counts.values, color='lightgreen', edgecolor='black')
    axs[1].set_title('Performance Category Distribution')
    axs[1].set_xlabel('Category')
    axs[1].set_ylabel('Count')
    axs[1].grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_metric_density_estimates(results: pd.DataFrame, metrics: list = None):
    """
    Plot density estimates of metrics for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to plot (default=None plots all)
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        sns.kdeplot(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                            fill=True, 
                                            color='skyblue', 
                                            alpha=0.7)
        plt.title(f'Density Estimate of {metric.upper()}')
        plt.xlabel(metric.upper())
        plt.ylabel('Density')
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_best_worst_simulation_characteristics(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot the characteristics of the best and worst simulations based on NSE.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Extract the best and worst results
    best_result = results.loc[results['nse'].idxmax()]
    worst_result = results.loc[results['nse'].idxmin()]
    
    params = ['k', 'n', 'l', 'r']
    x = np.arange(len(params))
    
    plt.figure(figsize=(12, 8))
    
    # Plot parameter values
    plt.subplot(2, 2, 1)
    plt.bar(x - 0.2, best_result[params], 0.4, label='Best Simulation', color='green')
    plt.bar(x + 0.2, worst_result[params], 0.4, label='Worst Simulation', color='red')
    plt.xticks(x, params)
    plt.title('Parameter Values for Best and Worst Simulations')
    plt.ylabel('Parameter Value')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    
    # Plot metric values
    plt.subplot(2, 2, 2)
    plt.bar(x - 0.2, best_result[['nse', 'r_squared', 'kge', 'pbias']], 0.4, label='Best Simulation', color='green')
    plt.bar(x + 0.2, worst_result[['nse', 'r_squared', 'kge', 'pbias']], 0.4, label='Worst Simulation', color='red')
    plt.xticks(x, ['NSE', 'R²', 'KGE', 'PBIAS'])
    plt.title('Metric Values for Best and Worst Simulations')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(axis='y', alpha=0.75)
    
    # Plot observed vs simulated for best and worst
    plt.subplot(2, 2, 3)
    plt.plot(obs, label='Observed', color='blue', linestyle='dashed')
    plt.plot(best_result['output'], label='Best Simulation', color='green')
    plt.plot(worst_result['output'], label='Worst Simulation', color='red')
    plt.title('Observed vs Simulated Values')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    
    plt.tight_layout()
    plt.show()

def plot_metric_trends_for_parameter(results: pd.DataFrame, obs: np.ndarray, parameter: str = 'k'):
    """
    Plot the trends of a metric over time for simulations grouped by a specific parameter.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - parameter: Parameter to group by (default='k')
    """
    import matplotlib.pyplot as plt
    
    # Extract unique parameter values
    param_values = results[parameter].unique()
    
    plt.figure(figsize=(10, 6))
    
    # Plot for each parameter value
    for value in param_values:
        subset = results[results[parameter] == value]
        plt.plot(subset['time'].iloc[0], subset['nse'], label=f'{parameter.upper()} = {value}')
    
    plt.title(f'Trend of NSE Over Time by {parameter.upper()}')
    plt.xlabel('Time')
    plt.ylabel('NSE')
    plt.legend()
    plt.grid()
    plt.show()

def plot_simulation_histories(results: pd.DataFrame, obs: np.ndarray, indices: list = None):
    """
    Plot the histories of selected simulations, including observed data.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - indices: List of indices of the simulations to plot (default=None plots best and worst)
    """
    import matplotlib.pyplot as plt
    
    if indices is None:
        # Default to best and worst simulations
        indices = [results['nse'].idxmax(), results['nse'].idxmin()]
    
    plt.figure(figsize=(10, 6))
    
    for index in indices:
        result = results.loc[index]
        plt.plot(result['time'], result['output'], label=f'Simulation {index}')
    
    plt.plot(obs, label='Observed', color='red', linestyle='dashed')
    plt.title('Simulation Histories')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid()
    plt.show()

def plot_parameter_effects_on_metrics(results: pd.DataFrame, metrics: list = None):
    """
    Plot the effects of parameters on model metrics using boxplots.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to analyze (default=None analyzes all)
    """
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    params = ['k', 'n', 'l', 'r']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        sns.boxplot(x=results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                             y=results['k'], 
                                             palette='coolwarm')
        plt.title(f'Effect of Parameter K on {metric.upper()}')
        plt.xlabel(metric.upper())
        plt.ylabel('Parameter K')
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_simulation_error_analysis(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot an analysis of simulation errors (residuals) for all simulations.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Calculate residuals for all simulations
    results['residuals'] = results.apply(lambda row: row['output'] - obs, axis=1)
    
    plt.figure(figsize=(10, 6))
    plt.hist(results['residuals'].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                         bins=30, 
                                         color='skyblue', 
                                         edgecolor='black')
    plt.title('Distribution of Simulation Residuals')
    plt.xlabel('Residual')
    plt.ylabel('Frequency')
    plt.grid(axis='y', alpha=0.75)
    plt.show()

def plot_metric_boxplots_by_parameter(results: pd.DataFrame, metrics: list = None):
    """
    Plot boxplots of metrics for simulations grouped by parameter values.
    
    Parameters:
    - results: DataFrame with simulation results
    - metrics: List of metrics to plot (default=None plots all)
    """
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    params = ['k', 'n', 'l', 'r']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        sns.boxplot(x=results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                             y=results['k'], 
                                             palette='coolwarm')
        plt.title(f'Boxplot of {metric.upper()} by Parameter K')
        plt.xlabel(metric.upper())
        plt.ylabel('Parameter K')
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_simulation_performance_overview(results: pd.DataFrame):
    """
    Plot an overview of simulation performance metrics.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    metrics = ['nse', 'r_squared', 'kge', 'pbias']
    performance_categories = ['Excellent', 'Very Good', 'Good', 'Satisfactory', 'Unsatisfactory']
    
    # Calculate performance metrics
    performance_metrics = results[metrics].apply(lambda row: categorize_performance(row['nse'], row['r_squared'], row['kge']), axis=1)
    results['performance_category'] = performance_metrics
    
    # Count of simulations in each category
    category_counts = results['performance_category'].value_counts().reindex(performance_categories, fill_value=0)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    # Performance category plot
    axs[0].bar(category_counts.index, category_counts.values, color='lightgreen', edgecolor='black')
    axs[0].set_title('Simulation Performance Category Distribution')
    axs[0].set_xlabel('Performance Category')
    axs[0].set_ylabel('Count')
    axs[0].grid(axis='y', alpha=0.75)
    axs[0].set_xticklabels(category_counts.index, rotation=45)
    
    # Metric distribution plots
    for i, metric in enumerate(metrics):
        axs[1].hist(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                           bins=20, 
                                           alpha=0.7, 
                                           label=metric.upper())
    
    axs[1].set_title('Distribution of Performance Metrics')
    axs[1].set_xlabel('Value')
    axs[1].set_ylabel('Frequency')
    axs[1].legend()
    axs[1].grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_metric_trends_with_observed(results: pd.DataFrame, obs: np.ndarray, metric: str = 'nse'):
    """
    Plot the trends of a metric over time for all simulations, with observed data.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metric: Metric to plot (default='nse')
    """
    import matplotlib.pyplot as plt
    
    # Extract time and metric values
    time = results['time'].iloc[0]
    metric_values = results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan)
    
    plt.figure(figsize=(10, 6))
    plt.plot(time, metric_values, label='Metric Trend', color='skyblue')
    plt.plot(time, obs, label='Observed', color='red', linestyle='dashed')
    plt.title(f'Trend of {metric.upper()} Over Time')
    plt.xlabel('Time')
    plt.ylabel(metric.upper())
    plt.legend()
    plt.grid()
    plt.show()

def plot_simulation_success_rate_by_category(results: pd.DataFrame):
    """
    Plot the success rate of simulations by performance category.
    
    Parameters:
    - results: DataFrame with simulation results
    """
    import matplotlib.pyplot as plt
    
    # Calculate success rate (percentage of successful simulations)
    success_rate = results['success'].value_counts(normalize=True).sort_index()
    
    # Performance categorization
    performance_categories = ['Excellent', 'Very Good', 'Good', 'Satisfactory', 'Unsatisfactory']
    category_counts = results['nse'].apply(lambda x: categorize_performance(x, 1.0)).value_counts().reindex(performance_categories, fill_value=0)
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    # Success rate plot
    axs[0].bar(success_rate.index.astype(str), success_rate.values, color='skyblue', edgecolor='black')
    axs[0].set_title('Simulation Success Rate')
    axs[0].set_xlabel('Success (1=True, 0=False)')
    axs[0].set_ylabel('Proportion')
    axs[0].grid(axis='y', alpha=0.75)
    axs[0].set_xticklabels(success_rate.index.astype(str), rotation=0)
    
    # Performance category plot
    axs[1].bar(category_counts.index, category_counts.values, color='lightgreen', edgecolor='black')
    axs[1].set_title('Performance Category Distribution')
    axs[1].set_xlabel('Category')
    axs[1].set_ylabel('Count')
    axs[1].grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_metric_density_with_observed(results: pd.DataFrame, obs: np.ndarray, metrics: list = None):
    """
    Plot density estimates of metrics for all simulations, with observed data.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    - metrics: List of metrics to plot (default=None plots all)
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if metrics is None:
        metrics = ['nse', 'r_squared', 'kge', 'pbias']
    
    plt.figure(figsize=(12, 8))
    
    for i, metric in enumerate(metrics, 1):
        plt.subplot(2, 2, i)
        sns.kdeplot(results[metric].apply(lambda x: x if isinstance(x, (int, float)) else np.nan), 
                                            fill=True, 
                                            color='skyblue', 
                                            alpha=0.7)
        sns.kdeplot(obs, fill=True, color='red', alpha=0.3, label='Observed')
        plt.title(f'Density Estimate of {metric.upper()}')
        plt.xlabel(metric.upper())
        plt.ylabel('Density')
        plt.legend()
        plt.grid(axis='y', alpha=0.75)
    
    plt.tight_layout()
    plt.show()

def plot_best_worst_simulation_outputs_comparison(results: pd.DataFrame, obs: np.ndarray):
    """
    Plot the output of the best and worst simulations based on NSE, in a comparative manner.
    
    Parameters:
    - results: DataFrame with simulation results
    - obs: Observed data for comparison
    """
    import matplotlib.pyplot as plt
    
    # Extract the best and worst results
    best_result = results.loc[results['nse'].idxmax()]
    worst_result = results.loc[results['nse'].idxmin()]
    
    # Time array (assuming same for all simulations)
    time = best_result['time']
    
    fig, axs = plt.subplots(2, 1, figsize=(10, 10))
    
    # Plot observed vs simulated for best and worst
    axs[0].plot(obs, label='Observed', color='blue', linestyle='dashed')
    axs[0].plot(best_result['output'], label='Best Simulation', color='green')
    axs[0].plot(worst_result['output'], label='Worst Simulation', color='red')
    axs[0].set_title('Observed vs Simulated Values')
    axs[0].set_xlabel('Time')
    axs[0].set_ylabel('Value')
    axs[0].legend()
    axs[0].grid()
    
    # Residuals plot
    axs[1].plot(best_result['time'], best_result['output'] - obs, label='Best Simulation Residuals', color='green')
    axs[1].plot(worst_result['time'], worst_result['output'] - obs, label='Worst Simulation Residuals', color='red')
    axs[1].axhline(0, color='blue', linestyle='dashed', linewidth=1, label='Zero Line')
    axs[1].set_title('Residuals of Best and Worst Simulations')
    axs[1].set_xlabel('Time')
    axs[1].set_ylabel('Residual')
    axs[1].legend()
    axs[1].grid()
    
    plt.tight_layout()
    plt.show()

def plot_parameter_correlation_scatter_matrix(results: pd.DataFrame, parameters: list = None):
    """
    Plot scatter plots and correlation coefficients between pairs of parameters.
    
    Parameters:
    - results: DataFrame with simulation results
    - parameters: List of parameters to include (default=None includes all)
    
    Returns:
    None: Displays the plots
    """
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    if parameters is None:
        parameters = ['k', 'n', 'l', 'r']
    
    # Select only the specified parameters
    subset = results[parameters]
    
    # Plot pairwise relationships
    sns.pairplot(subset, diag_kind='kde', markers='o', palette='coolwarm')
    plt.suptitle('Parameter Correlation Scatter Matrix', y=1.02)
    plt.show()

# =============================================================================
# MCMC CLASSES AND FUNCTIONS
# =============================================================================

class ParamsPriors:
    """Manages priors for all parameters in a layer."""
    def __init__(self, priors: Sequence[Prior]):
        self.prior_list = priors

    def sample(self):
        """Sample values for all parameters from their respective priors."""
        return Param(*(prior.sample() for prior in self.prior_list))

    def perturb(self, param: Param) -> Param:
        """Perturb each parameter using its respective prior."""
        return Param(*(prior.perturb(getattr(param, name)) for prior, name in zip(self.prior_list, ["REF_k", "REF_n", "REF_l", "REF_r"])))


# --- Layer Class ---
class Layer:
    """
    Represents a single layer with a name, depth (z_obs), and parameters.
    This `z_obs` corresponds to the observation depths, such as 0.1m, 0.2m, and 0.3m.
    """
    def __init__(self, name: str, z_obs: float, REF_k: float, REF_n: float, REF_l: float, REF_r: float):
        self.name = name  # Name of the layer (e.g., "Layer 1")
        self.z_obs = z_obs  # Depth of the observation, corresponds to zob (e.g., 0.1, 0.2, 0.3 meters)
        self.params = Param(REF_k, REF_n, REF_l, REF_r)  # Parameters for this layer

    def __repr__(self) -> str:
        return f"{self.name}: observed at {self.z_obs} m, with parameters {self.params}"


# --- MCMC Core Functions ---
def compute_energy(sim_temp: pd.DataFrame, obs_temp: pd.DataFrame, burn_in: int, sigma2: float) -> float:
    """
    Calculates the energy of the system, measuring the discrepancy between model and observed data.
    This function compares the simulated temperature at each depth to the observed data.
    
    :param sim_temp: Simulated temperature data
    :param obs_temp: Observed temperature data
    :param burn_in: The number of burn-in iterations
    :param sigma2: The variance of the noise in the temperature measurements
    :return: The energy of the system (i.e., the discrepancy between simulation and observation)
    """
    discrepancies = sim_temp.iloc[burn_in:].values - obs_temp.iloc[burn_in:].values
    norm2 = np.nansum(discrepancies ** 2)  # Sum of squared discrepancies
    return 0.5 * norm2 / sigma2  # Return energy, scaled by variance


def compute_log_acceptance(current_energy: float, previous_energy: float) -> float:
    """
    Computes the log acceptance ratio for the Metropolis-Hastings algorithm.
    This is used to decide whether to accept or reject the new set of parameters.
    
    :param current_energy: The energy computed with the new parameters
    :param previous_energy: The energy computed with the current parameters
    :return: The log acceptance ratio
    """
    return previous_energy - current_energy


def mcmc_step(layers: Sequence[Layer], observed_data: pd.DataFrame, sigma2: float, burn_in: int, priors: ParamsPriors, date_simul_bg, z_bottom, dz, nb_zone, alt_thk):
    """
    Performs a single MCMC iteration, including perturbation and acceptance decision.
    
    :param layers: List of Layer objects representing the model layers
    :param observed_data: DataFrame containing observed temperature data
    :param sigma2: The variance of the observed temperature data
    :param burn_in: The number of burn-in iterations
    :param priors: ParamPriors object containing the priors for model parameters
    :param date_simul_bg: Background date or time for the simulation
    :param z_bottom: Bottom depth for the simulation
    :param dz: Depth step for discretizing the model
    :param nb_zone: Number of zones or layers in the simulation
    :param alt_thk: Altered thickness for the model zones
    :return: The updated layers after perturbation
    """
    # 1. Perturb each layer's parameters
    perturbed_layers = [priors.perturb(layer.params) for layer in layers]

    # 2. Simulate temperature for the perturbed parameters
    perturbed_params = np.array([[layer.REF_k, layer.REF_n, layer.REF_l, layer.REF_r] for layer in perturbed_layers])
    
    # Assuming only 1 set of parameters (REF_k, REF_n, REF_l, REF_r) are used for now
    perturbed_model = run_direct_model(date_simul_bg, z_bottom, dz, nb_zone, alt_thk, *perturbed_params[0])  # Simulate temperature with perturbed parameters

    # 3. Compute energy for the current simulation and the perturbed one
    current_model = run_direct_model(date_simul_bg, z_bottom, dz, nb_zone, alt_thk, *np.array([[layer.params.REF_k, layer.params.REF_n, layer.params.REF_l, layer.params.REF_r] for layer in layers])[0])
    
    # Calculate the energy (discrepancy between observed and simulated temperatures)
    current_energy = compute_energy(current_model, observed_data, burn_in, sigma2)
    perturbed_energy = compute_energy(perturbed_model, observed_data, burn_in, sigma2)

    # 4. Metropolis-Hastings acceptance criterion
    log_acceptance = compute_log_acceptance(perturbed_energy, current_energy)
    acceptance_probability = min(1, np.exp(log_acceptance))

    # 5. Accept or reject the perturbed parameters based on the acceptance probability
    if uniform(0, 1) < acceptance_probability:
        # Accept perturbed parameters
        for i, layer in enumerate(layers):
            layer.params = perturbed_layers[i]
        accepted = True
    else:
        # Reject perturbed parameters
        accepted = False

    return layers, accepted


# --- Running MCMC ---
def run_mcmc(num_iterations: int, layers: Sequence[Layer], observed_data: pd.DataFrame, sigma2: float, burn_in: int, priors: ParamsPriors):
    """
    Runs the MCMC process for a number of iterations to sample the parameter space.
    
    :param num_iterations: Number of MCMC iterations
    :param layers: List of Layer objects representing the model layers
    :param observed_data: DataFrame containing observed temperature data
    :param sigma2: The variance of the observed temperature data
    :param burn_in: The number of burn-in iterations
    :param priors: ParamPriors object containing the priors for model parameters
    :return: The sampled parameter values from the MCMC process
    """
    accepted_parameters = []

    for iteration in range(num_iterations):
        print(f"Iteration {iteration + 1}/{num_iterations}")

        # Perform MCMC step
        layers, accepted = mcmc_step(layers, observed_data, sigma2, burn_in, priors, date_simul_bg, z_bottom, dz, nb_zone, alt_thk)
        
        # Store accepted parameters
        if accepted:
            accepted_parameters.append({layer.name: layer.params for layer in layers})
    
    return accepted_parameters

def evaluate_detailed_performance(merged_data, available_sensors):
    """
    Calculate detailed performance metrics for temperature sensors using manual calculations.
    
    This is a comprehensive fallback function when other statistical functions are not available.
    
    Parameters:
    - merged_data: DataFrame with columns like 'Temp1_sim', 'Temp1_obs', etc.
    - available_sensors: List of available sensor names (e.g., ['Temp1', 'Temp2', 'Temp3'])
    
    Returns:
    - Dictionary of metrics for each sensor
    """
    print("DETAILED PERFORMANCE METRICS (using stat_critere fallback implementation):")
    print("-" * 80)
    
    all_metrics = {}
    depth_map = {'Temp1': '-10 cm', 'Temp2': '-20 cm', 'Temp3': '-30 cm', 'Temp4': '-40 cm'}
    
    for sensor in available_sensors:
        obs_data = merged_data[f'{sensor}_obs'].dropna()
        sim_data = merged_data[f'{sensor}_sim'].dropna()
        
        # Ensure same length
        min_len = min(len(obs_data), len(sim_data))
        obs_data = obs_data.iloc[:min_len]
        sim_data = sim_data.iloc[:min_len]
        
        if len(obs_data) < 10:
            print(f"{sensor}: Insufficient data points ({len(obs_data)})")
            continue
            
        # Calculate basic statistics
        obs_mean = obs_data.mean()
        sim_mean = sim_data.mean()
        obs_std = obs_data.std()
        sim_std = sim_data.std()
        
        # Error metrics using stat_critere functions when available
        try:
            # Use the functions from this module
            rmse_val = rmse(sim_data.values, obs_data.values)
            mae_val = mae(sim_data.values, obs_data.values)
            bias_val = bias(sim_data.values, obs_data.values)
            pbias_val = pbias(sim_data.values, obs_data.values)
            nse_val = nse(sim_data.values, obs_data.values)
            kge_val = kge(sim_data.values, obs_data.values)
            d_index_val = d_index(sim_data.values, obs_data.values)
            correlation, p_value = correlation_coefficient(sim_data.values, obs_data.values)
        except:
            # Manual fallback calculations
            diff = sim_data - obs_data
            rmse_val = np.sqrt(np.mean(diff**2))
            mae_val = np.mean(np.abs(diff))
            bias_val = np.mean(diff)
            
            # Correlation metrics
            correlation, p_value = stats.pearsonr(obs_data, sim_data)
            
            # Nash-Sutcliffe Efficiency
            ss_res = np.sum((obs_data - sim_data)**2)
            ss_tot = np.sum((obs_data - obs_mean)**2)
            nse_val = 1 - (ss_res / ss_tot) if ss_tot != 0 else -np.inf
            
            # Percent Bias (PBIAS)
            pbias_val = 100 * np.sum(sim_data - obs_data) / np.sum(obs_data) if np.sum(obs_data) != 0 else np.inf
            
            # Index of Agreement (Willmott)
            numerator = np.sum((obs_data - sim_data)**2)
            denominator = np.sum((np.abs(sim_data - obs_mean) + np.abs(obs_data - obs_mean))**2)
            d_index_val = 1 - (numerator / denominator) if denominator != 0 else 0
            
            # Kling-Gupta Efficiency (KGE)
            r_component = correlation
            alpha_component = sim_std / obs_std if obs_std != 0 else np.inf
            beta_component = sim_mean / obs_mean if obs_mean != 0 else np.inf
            kge_val = 1 - np.sqrt((r_component - 1)**2 + (alpha_component - 1)**2 + (beta_component - 1)**2)
        
        # Calculate relative metrics
        r_squared = correlation**2
        rmse_rel = rmse_val / obs_mean * 100 if obs_mean != 0 else np.inf
        mae_rel = mae_val / obs_mean * 100 if obs_mean != 0 else np.inf
        
        # Store all metrics
        all_metrics[sensor] = {
            'RMSE': rmse_val,
            'MAE': mae_val,
            'Bias': bias_val,
            'RMSE_rel': rmse_rel,
            'MAE_rel': mae_rel,
            'R': correlation,
            'R²': r_squared,
            'NSE': nse_val,
            'PBIAS': pbias_val,
            'D_index': d_index_val,
            'KGE': kge_val,
            'p_value': p_value,
            'obs_mean': obs_mean,
            'sim_mean': sim_mean,
            'obs_std': obs_std,
            'sim_std': sim_std,
            'n_points': len(obs_data)
        }
        
        # Print detailed results for each sensor
        depth = depth_map.get(sensor, 'Unknown')
        
        print(f"{sensor} (at {depth}):")
        print(f"  Data points:           {len(obs_data)}")
        print(f"  RMSE:                  {rmse_val:.3f} °C ({rmse_rel:.1f}%)")
        print(f"  MAE:                   {mae_val:.3f} °C ({mae_rel:.1f}%)")
        print(f"  Bias:                  {bias_val:+.3f} °C")
        print(f"  Correlation (R):       {correlation:.3f} (p={p_value:.3e})")
        print(f"  R²:                    {r_squared:.3f}")
        print(f"  Nash-Sutcliffe (NSE):  {nse_val:.3f}")
        print(f"  Percent Bias (PBIAS):  {pbias_val:+.1f}%")
        print(f"  Index Agreement (d):   {d_index_val:.3f}")
        print(f"  Kling-Gupta (KGE):     {kge_val:.3f}")
        print()
    
    return all_metrics

def print_performance_summary(all_metrics, use_stat_critere=False):
    """
    Print a comprehensive performance summary table and overall assessment.
    
    Parameters:
    - all_metrics: Dictionary of metrics for each sensor
    - use_stat_critere: Boolean indicating which column names to use
    
    Returns:
    - Dictionary with overall statistics
    """
    if not all_metrics:
        print("No metrics available for summary.")
        return {}
    
    depth_map = {'Temp1': '-10 cm', 'Temp2': '-20 cm', 'Temp3': '-30 cm', 'Temp4': '-40 cm'}
    
    print("="*80)
    print("PERFORMANCE SUMMARY TABLE")
    print("="*80)
    print(f"{'Sensor':<8} {'Depth':<8} {'RMSE':<6} {'R²':<6} {'NSE':<6} {'PBIAS':<7} {'KGE':<6}")
    print("-" * 80)
    
    for sensor, metrics in all_metrics.items():
        depth = depth_map.get(sensor, 'Unknown')
        rmse_val = metrics.get('RMSE', 0)
        r2_val = metrics.get('R_squared', metrics.get('R²', 0))
        nse_val = metrics.get('NSE', 0)
        pbias_val = metrics.get('PBIAS', 0)
        kge_val = metrics.get('KGE', 0)
        
        print(f"{sensor:<8} {depth:<8} {rmse_val:.3f}  "
              f"{r2_val:.3f}  {nse_val:.3f}  "
              f"{pbias_val:+6.1f}% {kge_val:.3f