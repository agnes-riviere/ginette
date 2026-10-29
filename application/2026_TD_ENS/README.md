---
noteId: "3eecf830847a11f1aabba9711f319253"
tags: []

---

# 2026 TD ENS — Ginette : modele direct sur un cas synthetique

TD Ginette : mise en place d'un modele direct de transport de chaleur en zone hyporheique (echange nappe-riviere), sur un cas synthetique 1D.

- `TD_M2_direct_model.ipynb` : notebook du TD
- `SYNTHETIC_CASES/` : fichiers de configuration Ginette (`E_*.dat`). Les fichiers de sortie `Sim_*.dat` y sont ecrits a chaque simulation.
- `environment.yml` et `pyproject.toml` : installation de l'environnement (voir ci-dessous)

## Installation

Il faut **Python 3.12** et un compilateur Fortran (`gfortran`). Deux methodes sont possibles.

### Avec conda (recommande)

Conda installe Python 3.12 et gfortran dans un environnement isole :

```bash
cd application/2026_TD_ENS
conda env create -f environment.yml
conda activate ginette-td
```

### Avec pip

```bash
cd application/2026_TD_ENS
python3.12 -m venv .venv
# Linux / macOS :
source .venv/bin/activate
# Windows (PowerShell) :
# .venv\Scripts\activate
pip install -e .
```

Avec pip, `gfortran` doit etre installe et accessible depuis le terminal :

- Linux : `sudo apt install gfortran`
- macOS : `brew install gcc` (gfortran fait partie de gcc)
- Windows : installer [MSYS2](https://www.msys2.org/), puis `pacman -S mingw-w64-x86_64-gcc-fortran`, et ajouter le dossier `mingw64/bin` au `PATH`

Le notebook compile automatiquement l'executable `ginette` au premier lancement, s'il n'existe pas deja.

## Lancer le TD

Lancer Jupyter **depuis ce dossier** (`2026_TD_ENS`), pour que les chemins relatifs du notebook se resolvent correctement :

```bash
jupyter notebook TD_M2_direct_model.ipynb
```
