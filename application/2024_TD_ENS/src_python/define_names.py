from enum import Enum
import json

# Temporal values
class TemporalValues:
    NSECINMIN = 60
    NSECINHOUR = 3600
    NSECINDAY = 86400
    NHOURINDAY = 24
    NDAYINYEAR = 365
    NDAYINMONTH = 30
    ABSURD_DATE = "1999/09/09 09:09:09"

# Default absurd values
class DefaultValues:
    CODE_TEMP = 959595
    CODE_SCALAR = -9999
    PARAMBOUND = 1e7


# Prior initialization
class PriorInitialization:
    permeability_range = (11, 15)
    permeability_SIGMA = 0.01

    N_INTERVAL = (0.01, 0.25)
    N_SIGMA = 0.01

    LAMBDA_S_INTERVAL = (1, 5)
    LAMBDA_S_SIGMA = 0.1

    RHOS_CS_INTERVAL = (1e6, 1e7)
    RHOS_CS_SIGMA = 1e5


# MCMC parametrization
NITMCMC = 1000

# Save parameter names to JSON files
def save_parameter_names():
    names_phy = ['lambda_s', 'rho_s', 'c_s', 'n', 'permeability']
    with open('names_phy.json', 'w') as f:
        json.dump(names_phy, f)

    names_red = ['alpha_e', 'kappa_e']
    with open('names_red.json', 'w') as f:
        json.dump(names_red, f)

    names_BC = ["period_T", "dH", "A", "T_mu", "sampling_period"]
    with open('names_BC.json', 'w') as f:
        json.dump(names_BC, f)

