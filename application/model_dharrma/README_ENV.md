# Environnement Virtuel Model DHARRMA

## Utilisation

### Activation de l'environnement
```bash
cd /home/ariviere/Programmes/ginette/application/model_dharrma
source activate_env.sh
```

### Installation de packages
Une fois l'environnement activé :
```bash
pip install <package_name>
```

### Désactivation
```bash
deactivate
```

### Utilisation directe sans activation
Vous pouvez aussi utiliser l'environnement directement sans l'activer :
```bash
# Python
./venv/bin/python script.py

# Pip
PYTHONPATH="" ./venv/bin/pip install <package>
```

## Structure
```
model_dharrma/
├── venv/                 # Environnement virtuel
├── activate_env.sh       # Script d'activation
└── README_ENV.md         # Ce fichier
```

## Notes importantes
- L'environnement utilise Python 3.11.14
- Pip 25.3 est installé et configuré correctement
- Le script d'activation nettoie automatiquement PYTHONPATH pour éviter les conflits avec les packages système