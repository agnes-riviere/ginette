#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 15:21:53 2025

@author: Maxime GAUTIER
"""


# %% IMPORTS:
import os
import numpy as np
import pandas as pd
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.pyplot as plt
cmap = matplotlib.colormaps['nipy_spectral'].resampled(256).reversed()
fontsize = 12
plt.rcParams["font.size"] = fontsize


# %% LOAD RESULTS:
results = pd.read_csv(os.path.join("results.txt"), delimiter=" ",
                      index_col=[0])


assess_var = "misfit_tot"
results.sort_values(by=assess_var, inplace=True, ascending=True)
best = results[results[assess_var] == results[assess_var].min()].ID.values[0]
sim_data = pd.read_csv(os.path.join("results", f"sim_temp_{best}.txt"),
                       delimiter=" ", index_col=[0])
obs_data = pd.read_csv(os.path.join("observed_data.txt"), delimiter=" ",
                       index_col=[0])

# Marginals 1D:
marginal_log_k = results.groupby("log_k")[assess_var].sum()
marginal_lam = results.groupby("lam")[assess_var].sum()


# %% PLOT:
fig = plt.figure(figsize=(16, 9), dpi=100)
gs = gridspec.GridSpec(2, 2, figure=fig)
axa = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[1, 0])
axc = fig.add_subplot(gs[0, 1])
axd = fig.add_subplot(gs[1, 1])

axa.plot(np.log10(marginal_log_k.index), marginal_log_k.values, lw=2, ls="-",
         color="black", label="Marginal 1D")
axa.axvline(-12, lw=4, ls="--", color="lime", label="True value")
axa.axvline(np.log10(results.loc[best, "log_k"]), lw=3, ls="--", color="gold",
            label="Best value")
axa.legend(loc="upper left")
axa.set_xlabel(r"$log_{10}(k)$")
axa.set_ylabel(r"$\Phi{}_{tot.}$")
axa.grid()

axd.plot(marginal_lam.index, marginal_lam.values, lw=2, ls="-",
         color="black", label="Marginal 1D")
axd.axvline(2.5, lw=4, ls="--", color="lime", label="True value")
axd.axvline(results.loc[best, "lam"], lw=3, ls="--", color="gold",
            label="Best value")
axd.legend(loc="upper left")
axd.set_xlabel(r"$\lambda\;(W/m°C)$")
axd.set_ylabel(r"$\Phi{}_{tot.}$")
axd.grid()

# Plot 2D:
gci = axb.scatter(np.log10(results.log_k), results.lam, c=results.misfit_tot,
                  cmap=cmap)
fig.colorbar(gci, ax=axb, label=r"$\Phi{}_{tot.}$")
axb.set_xlabel(r"$log_{10}(k)$")
axb.set_ylabel(r"$\lambda\;(W/m°C)$")
axb.grid()
axb.set_xlim(np.log10(results.log_k).min(), np.log10(results.log_k).max())
axb.set_ylim(results.lam.min(), results.lam.max())
axb.scatter(-12, 2.5, marker="o", edgecolors="lime", facecolors="None", s=100,
            lw=5, label="True values")
axb.scatter(np.log10(results.loc[best, "log_k"]),
            results.loc[best, "lam"], marker="*", edgecolors="black",
            facecolors="gold", s=80, lw=1, label="Best")
axb.legend(loc="upper left")

# Plot temperatures:
axc.plot(obs_data.Time, obs_data.Temp1, ls="-", lw=4, color="red",
         label="Temperature 1 (-10 cm)")
axc.plot(obs_data.Time, obs_data.Temp2, ls="-", lw=4, color="green",
         label="Temperature 2 (-20 cm)")
axc.plot(obs_data.Time, obs_data.Temp3, ls="-", lw=4, color="blue",
         label="Temperature 3 (-30 cm)")
axc.plot(sim_data.Time, sim_data.Temp1, ls="--", lw=2, color="orange",
         label="Sim. temp. 1")
axc.plot(sim_data.Time, sim_data.Temp2, ls="--", lw=2, color="lime",
         label="Sim. temp. 2")
axc.plot(sim_data.Time, sim_data.Temp3, ls="--", lw=2, color="dodgerblue",
         label="Sim. temp. 3")
axc.set_xlabel("Time (s)")
axc.set_ylabel("Temperature (°C)")
axc.grid()
axc.legend(loc="lower right")



