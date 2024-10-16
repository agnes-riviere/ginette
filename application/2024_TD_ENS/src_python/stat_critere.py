import numpy as np
import pandas as pd



def rmse(predictions, ref,sigma):
    differences = ((predictions - ref)/sigma)                      #the DIFFERENCEs.
    differences_squared = differences ** 2                    #the SQUAREs of ^
    mean_of_differences_squared = differences_squared.mean()  #the MEAN of ^
    rmse_val = np.sqrt(mean_of_differences_squared)           #ROOT of ^

    return rmse_val


import numpy as np

def calculate_metrics(data):
    metrics = {}

    for i in range(1, 4):
        sensor = f'Temp{i}'
        diff = data[f'{sensor}_sim'] - data[f'{sensor}_obs']
        
        rmse = np.sqrt(np.mean(diff**2))  # RMSE
        mae = np.mean(np.abs(diff))  # MAE
        bias = 100 * np.sum(diff) / np.sum(np.abs(data[f'{sensor}_obs']))  # PBias

        std_obs = data[f'{sensor}_obs'].std()
        std_sim = data[f'{sensor}_sim'].std()
        mean_obs = data[f'{sensor}_obs'].mean()
        mean_sim = data[f'{sensor}_sim'].mean()

        kge = 1 - np.sqrt((np.square(std_sim/std_obs - 1)) +
                            (np.square(mean_sim/mean_obs - 1)) +
                            (np.square(rmse/std_obs - 1)))

        metrics[sensor] = {'RMSE': rmse, 'MAE': mae, 'PBias': bias, 'KGE': kge}

    return metrics
def compute_energy(temp_obs, temp_sim, variance: float, remanence):
    """
    Calcule l'énergie entre deux séries temporelles de température.

    Parameters:
    - temp_obs:  série temporelle de température observée.
    - temp2:  série temporelle de température simulée.
    - variance variance.
    - remanence: facteur de rémanence.

    Returns:
    - énergie calculée.
    """
    remanence_step = int(remanence * 24 * 60 / 15) # 24 heures * 60 minutes / 15 minutes
    norm2 = np.nansum((temp_obs[:, remanence_step:] - temp_sim[:, remanence_step:]) ** 2)
    return 0.5 * norm2 / variance