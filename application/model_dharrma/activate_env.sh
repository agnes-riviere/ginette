#!/bin/bash

# Script d'activation de l'environnement virtuel pour model_dharrma
# Usage: source activate_env.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_PATH="$SCRIPT_DIR/venv"

if [ -d "$VENV_PATH" ]; then
    echo "Activation de l'environnement virtuel model_dharrma..."
    export VIRTUAL_ENV="$VENV_PATH"
    export PATH="$VENV_PATH/bin:$PATH"
    export PYTHONPATH=""  # Nettoyer PYTHONPATH pour éviter les conflits
    
    # Modifier le prompt pour indiquer l'environnement actif
    if [ -z "${VIRTUAL_ENV_DISABLE_PROMPT:-}" ]; then
        PS1="(model_dharrma) $PS1"
        export PS1
    fi
    
    echo "Environnement virtuel activé !"
    echo "Python: $(which python)"
    echo "Pip: $(which pip)"
    python --version
    pip --version
else
    echo "Erreur: Environnement virtuel non trouvé dans $VENV_PATH"
    return 1
fi