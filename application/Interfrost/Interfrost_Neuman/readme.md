---
noteId: "d156165001be11f1b9103f30f863e303"
tags: []

---

# InterFrost TH1 Test Case - Heat Advection Benchmark

## Overview

Heat advection via constant pore water velocity is here added to heat exchange due to conduction and phase change. The analytical solution is obtained from the **Kurylyk et al. (2014, AWR)** paper based on a reassessment of solutions by **Lunardini (1998, Proc. 7th. Conf. Permf.)**. 

> Although not physically realistic (constant velocity independent of ice saturation) this solution can be used for benchmarking purposes: *"Lack of fidelity to physical processes does not limit ability to serve as a benchmark"*.

## Reference

The paper **Kurylyk et al. (2014, AWR)** details several alternative analytical solutions, the suggested benchmark cases, and the SUTRA code runs to compare with these solutions. 

**Full Reference:**
> "Analytical solutions for benchmarking cold regions subsurface water flow and energy transport models: One-dimensional soil thaw with conduction and advection" by B. Kurylyk, J. McKenzie, K. MacQuarrie, C. Voss in *Advances in Water Resources* 70 (2014) 172–184.

A presentation of the TH1 Case was made by **Barret Kurylyk** during the kick off meeting.

## Benchmark Cases

### TH1 Cases (Benchmarks 2 and 3)
Benchmark cases **2 and 3** recommended by Kurylyk et al. (2014) are included in the InterFrost project as two TH1 cases differing by the flow velocity considered:

- **Case 2:** Velocity of **10 m/yr**
- **Case 3:** Velocity of **100 m/yr**

These benchmarks are based on the following analytical solution where:
- **X** = depth to the thawing front
- **α** = thermal diffusivity  
- **v<sub>t</sub>** = thermal plume velocity due to advection
- **t** = time
- **S<sub>T</sub>** = dimensionless Stefan number

### Benchmark 1 (Neumann Solution)
**Benchmark 1** (Neumann solution) recommended by Kurylyk et al. (2014) is an option to complement the Lunardini case provided in the InterFrost project as **T1**. 

⚠️ **Note:** The initial temperature (T<sub>i</sub>) term in the Neumann solution is expressed as the number of degrees below 0°C (i.e., it is a positive number).

## TH1 Test Case Description

The initially uniform temperature is at **0°C**. This condition simplifies the energy balance at the thawing front by ensuring that there is no thermal gradient (or conductive flux) below the thawing front. 

### Key Features:
- Specified surface temperature is **above 0°C** and hence induces thaw
- Water is advected through the entire medium
- Divergence of the advective flux is **zero below the thawing front** due to uniform thermal conditions

## Initial and Boundary Conditions

### Initial Conditions
- **Temperature:** Uniform at 0°C throughout the domain
- **Phase state:** Initially frozen ground

### Boundary Conditions
- **Surface (z=0):** Constant temperature **T<sub>s</sub> > 0°C**
- **Bottom boundary:** Adiabatic or constant temperature at 0°C
- **Flow:** Constant pore water velocity throughout the domain

## Physical Parameters

### Thermal Properties
- **Thermal conductivity (frozen):** λ<sub>f</sub> = 2.619 W/(m·K)
- **Thermal conductivity (unfrozen):** λ<sub>u</sub> = 1.839 W/(m·K)
- **Heat capacity:** C = 2.0 × 10<sup>6</sup> J/(m³·K)
- **Latent heat of fusion:** L = 3.34 × 10<sup>8</sup> J/m³

### Hydraulic Properties
- **Porosity:** n = 0.4
- **Water density:** ρ<sub>w</sub> = 1000 kg/m³
- **Specific heat of water:** c<sub>w</sub> = 4180 J/(kg·K)

## Analytical Solution

The analytical solution provides the position of the thawing front as a function of time and can be used to validate numerical models. The solution accounts for:

1. **Conductive heat transfer** in both frozen and unfrozen zones
2. **Advective heat transport** due to groundwater flow
3. **Phase change** at the freezing front

## Implementation in Ginette

This benchmark is implemented in the Ginette code using:
- **NEUMA thermal conductivity model** for temperature-dependent conductivity
- **Gel/dégel physics** (icycle=1) for phase change modeling
- **Thermal transport equations** with conduction and advection terms

### Input Files

The TH1 test case requires the following input files for ginette_V2.f90:

#### Main Configuration
- `E_parametre.dat` - Main parameter file including:
  - `ytest=TH1` (line 54) - Activates TH1 benchmark mode
  - `ysolv=BIC` (line 79) - Solver selection (BIC or CGS)
  - Time stepping and hydraulic parameters
- `E_p_therm.dat` - Thermal parameters including:
  - `icycle=1` - Activates freeze/thaw physics
  - `ymoycondtherm=NEUMA` - Thermal conductivity model

#### Boundary Conditions and Geometry
- `E_cdt_aux_limites.dat` - Thermal boundary conditions (surface temperature +5°C) and Flow boundary conditions
- E_cdt_initiale.dat - Initial condition = T=-5°C

### Output Files

The TH1 benchmark generates specific output files:

#### Time Series Data
- `S_bound_permaf_1_t.dat` - Permafrost boundary positions:
  - Format: `time zs(1,1) zs(1,2) zl(1,1) zl(1,2)`
  - Tracks freeze/thaw front locations
- `S_Temp_id400_t.dat` - Temperature at specific node (400):
  - Format: `time temp(400)`
  - Monitoring point temperature evolution
- `S_bilan_therm_t.dat` - Thermal balance:
  - Format: `time qtherm(1)`
  - Heat flux at surface boundary

#### Spatial Profiles
- `S_pts` - Complete profile data (binary format):
  - Pressure, temperature, and saturation profiles
  - Updated at each output time step
- Binary output (unit 181818) - Compact format:
  - Format: `time gxl gxs`
  - Freeze/thaw front positions

#### Standard Output
- Console output with `OUT` statements showing:
  - Time step information
  - Temperature evolution at surface node
  - Record indicators

### TH1-Specific Features

The TH1 implementation includes several specialized features:

#### Physical Model Settings
- **Forced thaw mode**: `igel = 2` (always in thaw/degel state)
- **Water density**: Uses liquid water density for thaw calculations
- **Thermal conductivity**: NEUMA model with temperature dependence
  - λ<sub>frozen</sub> = 2.619 W/(m·K)
  - λ<sub>unfrozen</sub> = 1.839 W/(m·K)

#### Computational Features
- **Front tracking**: Automatic calculation of `gxl` and `gxs` parameters
  - `gxl = 2.0 - zl(1,1)` (liquid front position)
  - `gxs = 2.0 - zs(1,1)` (solid front position)
- **Thermal matrix**: Specialized assembly for freeze/thaw physics
- **Iterative solver**: BIC (BiCGSTAB) or CGS solver selection

#### Validation Outputs
- **Position tracking**: Continuous monitoring of freeze/thaw boundary
- **Temperature profiles**: Full spatial temperature distribution
- **Heat balance**: Surface heat flux for energy conservation check
- **Benchmark comparison**: Direct comparison with Neumann analytical solution

---

**Related InterFrost Documentation:** [Test Cases](https://wiki.lsce.ipsl.fr/interfrost/interfrost/doku.php?id=test_cases:seven.html)

*This benchmark is part of the InterFrost project for validating permafrost modeling codes.*