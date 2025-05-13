# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 14:57:23 2023
"""
import numpy as np
import pandas as pd
import os
from src.src_python.Direct_model import setup_ginette, generate_zone_parameters, initial_conditions, boundary_conditions, run_direct_model, remove_first_two_days



def calcul_marginale(df, parametre, erreur_colonne='Erreur'):
    # Regrouper par le paramètre choisi et sommer les erreurs
    marginale = (
        df.groupby(parametre)[erreur_colonne]
        .sum()
        .reset_index()
        .sort_values(by=parametre)  # Trier par ordre croissant
    )
    return marginale



def generer_string(pert,sigma_post, Lchaine, variation_k):
  
    # Obtenir la date et l'heure actuelles
    date_heure = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Générer la chaîne de paramètres
    variation_params = ", ".join([k for k in variation_k])
    
    # Construire la chaîne finale
    resultat = (
        f"{date_heure}|"
        f"pertu_{pert} | "
        f"sigma_post{sigma_post}"
        f"Lchaine_{Lchaine} | "
        f"var_{variation_params} | "
        
    )
    
    return resultat

def norme_2(vect1,vect2,sigma):
    """
    Calcule la somme des carrés des différences entre deux vecteurs.
    
    Parameters:
        vect1 (list ou numpy.array) : Premier vecteur
        vect2 (list ou numpy.array) : Deuxième vecteur
        
    Returns:
        float : La somme des carrés des différences
    """
    if len(vect1) != len(vect2):
        raise ValueError("Les deux vecteurs doivent avoir la même taille.")
    
    return (sum((a - b) ** 2 for a, b in zip(vect1, vect2)))/sigma**2



def read_input_metro():
    """
    Reads the file 'read_input_metro.txt'.

    The file should contain:
    - First line: Number of iterations (integer)
    - Second line: Uncertainties of the temperature measurements (space-separated floats)
    - Third line: Number of zones (integer)
    - Fourth line: List of parameter names (space-separated strings, e.g., ['k', 'n', 'lambda', 'c'])

    Returns:
    - iterations (int): Number of iterations
    - uncertainties (list of float): List of uncertainties
    - nb_zone (int): Number of zones
    - parameters (list of str): List of parameter names
    """
    try:
        with open("read_input_metro.txt", 'r') as f:
            lines = f.readlines()
            if len(lines) < 4:
                raise ValueError("Error: 'read_input_metro.txt' must contain at least 4 lines.")
            iterations = int(lines[0].strip())
            uncertainties = list(map(float, lines[1].strip().split()))
            nb_zone = int(lines[2].strip())
            parameters = lines[3].strip().split()
    except FileNotFoundError:
        raise FileNotFoundError("Error: 'read_input_metro.txt' file not found.")
    except ValueError as e:
        raise ValueError(f"Error processing file: {e}")

    return iterations, uncertainties, nb_zone, parameters




def read_bounds_params(nb_zone, parameters):
    """
    Reads the file 'bound_params.txt' containing zone, parameter, min_value, max_value, and std.

    Parameters:
    - nb_zone: Number of zones
    - parameters: List of parameter names (e.g., ['k', 'n', 'lambda', 'c'])

    Returns:
    - A list of dictionaries, one for each zone, where each dictionary maps parameters to their bounds and std.
    """
    params_dict = {zone: {param: None for param in parameters} for zone in range(1, nb_zone + 1)}

    try:
        with open("bound_params.txt", 'r') as f:
            for line in f:
                if not line.strip():
                    continue

                columns = line.strip().split()
                if len(columns) != 5:
                    print(f"Invalid line format: {line}")
                    continue

                zone, param, min_value, max_value, std = int(columns[0]), columns[1], float(columns[2]), float(columns[3]), float(columns[4])

                if 1 <= zone <= nb_zone and param in parameters:
                    params_dict[zone][param] = (min_value, max_value, std)
                else:
                    print(f"Invalid entry: {line}")

    except FileNotFoundError:
        print("Error: 'bound_params.txt' file not found.")
    except ValueError as e:
        print(f"Error processing line: {line}. {e}")

    return [params_dict[zone] for zone in range(1, nb_zone + 1)]




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
            new_val = self.range[0] + (new_val - self.range[1]) % (self.range[1] - self.range[0])
        while new_val < self.range[0]:
            new_val = self.range[1] - (self.range[0] - new_val) % (self.range[1] - self.range[0])
        return new_val

    def sample(self) -> float:
        """Samples a random value within the specified range."""
        return uniform(*self.range)

    def __repr__(self) -> str:
        return f"Prior(range={self.range}, sigma={self.sigma})"


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

        # Perform a single MCMC step
        layers, accepted = mcmc_step(layers, observed_data, sigma2, burn_in, priors)
        
        if accepted:
            accepted_parameters.append([layer.params for layer in layers])

    return accepted_parameters

