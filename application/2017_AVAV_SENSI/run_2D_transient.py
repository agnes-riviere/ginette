from __future__ import annotations

"""
Lance un cas AVA en cinq etapes:
1. un calcul permanent sans zone non saturee (ivg=0, yunconfined=UNC),
2. un calcul permanent avec ZNS et des parametres Van Genuchten tres adoucis (EASY),
3. un calcul permanent avec ZNS et des parametres Van Genuchten intermediaires (MEDIUM),
4. un calcul permanent avec ZNS et les vrais parametres Van Genuchten,
5. un calcul transitoire qui reutilise cette pression finale comme condition initiale.

Pourquoi ce script existe:
- sur ce cas AVA, on souhaite d'abord approcher un etat hydraulique stable,
- en commencant par une version plus simple sans ZNS,
- en reintroduisant la ZNS avec une courbe de retention tres douce (EASY, convergence Picard facile),
- puis avec une courbe intermediaire (MEDIUM, alpha/n a mi-chemin entre EASY et les vraies valeurs),
- puis en passant aux vrais parametres Van Genuchten de la geologie du site,
- puis demarrer le calcul transitoire depuis cet etat,
- sans activer la temperature, car elle n'est pas calculee pour ce cas.

Le script modifie donc quelques fichiers de controle Ginette entre les etapes:
- E_parametre.dat
- E_p_therm.dat
- E_cdt_initiale.dat
- E_cdt_aux_limites.dat
- E_pression_initiale.dat
- E_zone_parameter.dat  (sauvegarde + parametres adoucis + restauration)
"""

import subprocess
import sys
import importlib.util
import os
import time
from pathlib import Path

import pandas as pd


CASE_DIR = Path(__file__).resolve().parent
REPO_DIR = CASE_DIR.parents[1]

STEADY_TIMEOUT_SEC = int(os.environ.get("GINETTE_STEADY_TIMEOUT_SEC", "3600"))
TRANSIENT_TIMEOUT_SEC = int(os.environ.get("GINETTE_TRANSIENT_TIMEOUT_SEC", "3600"))
PROGRESS_EVERY_SEC = int(os.environ.get("GINETTE_PROGRESS_EVERY_SEC", "20"))
TRACE_LOG = CASE_DIR / "run_2D_transient.trace.log"

# Parametres Van Genuchten "adoucis" utilises pour l'etape intermediaire.
# L'objectif est d'aider la convergence de Picard en evitant les courbes de
# retention trop raides lors du premier calcul permanent ZNS.
#
# Colonnes (ordre impose par Ginette, cf. lecture unité 321 dans ginette_V2.f90) :
#   zone_id, ak [m/s], porosity, alpha [m-1], n, swres, lambda [W/m/K], cpm [J/kg/K], rho [kg/m3]
#
# Strategie : on garde ak, porosite et proprietes thermiques identiques aux
# vraies valeurs, mais on reduit alpha et n, et on augmente swres.
# Cela aplatit la courbe S(h) => derivee ds/dh plus faible => iterations Picard
# moins explosives quand la nappe est proche de la surface.
#
# MODIFIER CES VALEURS si un autre jeu de parametres doux convient mieux.
EASY_ZONE_PARAMS: dict[int, tuple] = {
    # id  ak          porosity   alpha    n      swres    lambda  cpm    rho
    1:  (1.00e-10,  0.040,  1.0,  1.5,  0.10,  2.24,  827.0,  2040.0),
    2:  (7.50e-14,  0.200,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    3:  (1.00e-09,  0.010,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    4:  (1.00e-11,  0.010,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    5:  (1.00e-11,  0.010,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    6:  (1.00e-11,  0.020,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    7:  (1.00e-10,  0.020,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
    8:  (1.00e-11,  0.200,  1.0,  1.5,  0.05,  0.23,  860.0,  2600.0),
}

# Parametres Van Genuchten intermediaires (a mi-chemin entre EASY et les vraies valeurs).
# Zone 1 : vraies valeurs alpha=14.5 n=2.68 => medium alpha~4, n~2.0
# Zones 2-8 : vraies valeurs alpha=6.2  n=2.0  => medium alpha~2.5, n~1.75
# Le swres et les proprietes thermiques sont laisses proches des vraies valeurs.
MEDIUM_ZONE_PARAMS: dict[int, tuple] = {
    # id  ak          porosity   alpha    n      swres    lambda  cpm    rho
    1:  (1.00e-10,  0.040,  4.0,  2.0,  0.12,  2.24,  827.0,  2040.0),
    2:  (7.50e-14,  0.200,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    3:  (1.00e-09,  0.010,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    4:  (1.00e-11,  0.010,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    5:  (1.00e-11,  0.010,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    6:  (1.00e-11,  0.020,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    7:  (1.00e-10,  0.020,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
    8:  (1.00e-11,  0.200,  2.5,  1.75, 0.025, 0.23,  860.0,  2600.0),
}


def log_trace(message: str) -> None:
	"""Ecrit un message horodate dans le journal du script."""
	timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
	line = f"[{timestamp}] {message}"
	print(line)
	with open(TRACE_LOG, "a", encoding="utf-8") as trace_file:
		trace_file.write(line + "\n")

def load_module(module_name: str, module_path: Path):
	"""
	Charge un module Python directement depuis son fichier.

	On evite ainsi `import src...`, car le package `src` importe aussi
	`src_gmsh`, qui depend de `pyvista`. Ce script n'a pas besoin de cette
	partie du projet.
	"""
	spec = importlib.util.spec_from_file_location(module_name, module_path)
	if spec is None or spec.loader is None:
		raise ImportError(f"Could not load module {module_name} from {module_path}")
	module = importlib.util.module_from_spec(spec)
	sys.modules[module_name] = module
	spec.loader.exec_module(module)
	return module


DIRECT_MODEL = load_module("ginette_direct_model", REPO_DIR / "src" / "src_python" / "Direct_model.py")
INIT_FOLDERS = load_module("ginette_init_folders", REPO_DIR / "src" / "src_python" / "Init_folders.py")

boundary_conditions_2D_tdirect = DIRECT_MODEL.boundary_conditions_2D_tdirect
boundary_conditions_perm_2D_tdirect = DIRECT_MODEL.boundary_conditions_perm_2D_tdirect
format_value = DIRECT_MODEL.format_value
compile_ginette_src = INIT_FOLDERS.compile_ginette_src


def replace_setting(file_path: Path, prefix: str, new_value: str) -> None:
	"""
	Remplace la valeur d'un parametre dans un fichier texte de configuration.

	Exemple:
	- prefix = "rp="
	- new_value = "0"
	transforme la ligne commencant par "rp=" en conservant le reste du commentaire.

	Cette fonction est pratique pour les fichiers Ginette, car ils contiennent
	une valeur suivie d'un commentaire explicatif sur la meme ligne.
	"""
	lines = file_path.read_text().splitlines()
	updated_lines = []
	replaced = False

	for line in lines:
		# On ne modifie que la ligne qui commence exactement par le prefixe demande.
		if line.startswith(prefix):
			value_start = len(prefix)
			value_end = value_start
			# La valeur s'arrete au premier espace ou tabulation.
			while value_end < len(line) and not line[value_end].isspace():
				value_end += 1
			suffix = line[value_end:]
			# On reconstruit la ligne en gardant intacte la partie commentaire.
			updated_lines.append(f"{prefix}{new_value}{suffix}")
			replaced = True
		else:
			updated_lines.append(line)

	if not replaced:
		raise ValueError(f"Setting '{prefix}' not found in {file_path}")

	file_path.write_text("\n".join(updated_lines) + "\n")


def read_first_scalar(file_path: Path) -> float:
	"""
	Lit la premiere valeur numerique d'un fichier colonne.

	On s'en sert ici pour recuperer les premieres charges imposees sur:
	- la rive droite,
	- la rive gauche,
	- la riviere.

	Ces valeurs servent a construire un run permanent coherent avec
	les conditions imposées au debut du transitoire.
	"""
	data = pd.read_csv(file_path, sep=r"\s+", header=None)
	if data.empty:
		raise ValueError(f"File is empty: {file_path}")
	return float(data.iloc[0, 0])


def read_nmi(file_path: Path) -> int:
	"""
	Extrait le nombre de mailles actives `nmi` depuis E_parametre.dat.

	Cette information est indispensable pour retrouver, dans S_pression_charge_temperature.dat,
	le dernier bloc complet de pression correspondant a toutes les mailles.
	"""
	for line in file_path.read_text().splitlines():
		if line.startswith("nmi="):
			raw_value = line.split()[0].split("=", 1)[1]
			return int(raw_value)
	raise ValueError(f"Could not find nmi in {file_path}")


def backup_zone_parameters() -> None:
	"""
	Sauvegarde E_zone_parameter.dat sous le nom E_zone_parameter.dat.orig.

	Cette copie sert a restaurer les vrais parametres de zone apres l'etape
	intermediaire qui utilise des parametres Van Genuchten adoucis.
	"""
	src = CASE_DIR / "E_zone_parameter.dat"
	dst = CASE_DIR / "E_zone_parameter.dat.orig"
	dst.write_text(src.read_text())


def restore_zone_parameters() -> None:
	"""
	Restaure E_zone_parameter.dat depuis la sauvegarde E_zone_parameter.dat.orig.

	Appele apres l'etape intermediaire pour que le permanent ZNS final
	utilise bien les vrais parametres geologiques du site.
	"""
	backup = CASE_DIR / "E_zone_parameter.dat.orig"
	if not backup.exists():
		raise FileNotFoundError(
			"E_zone_parameter.dat.orig introuvable — backup_zone_parameters() n'a pas ete appele"
		)
	(CASE_DIR / "E_zone_parameter.dat").write_text(backup.read_text())


def write_easy_zone_parameters() -> None:
	"""
	Ecrit les parametres Van Genuchten adoucis dans E_zone_parameter.dat.

	Les valeurs proviennent du dictionnaire EASY_ZONE_PARAMS defini en tete
	de script. On conserve le meme nombre de zones et le meme format de
	fichier (liste libre separee par des espaces) que le fichier original.

	Fortran lit ce fichier avec un `read(321, *)` (lecture libre),
	done les valeurs Python en notation scientifique standard sont acceptees.
	"""
	lines = []
	for zone_id, params in sorted(EASY_ZONE_PARAMS.items()):
		ak, om, asp, ans, swres, alambda, cpm, rho = params
		lines.append(
			f"{zone_id}\t{ak:.2E}\t{om:.2E}\t{asp:.2E}\t{ans:.2E}"
			f"\t{swres:.2E}\t{alambda:.2E}\t{cpm:.0f}\t{rho:.2E}"
		)
	(CASE_DIR / "E_zone_parameter.dat").write_text("\n".join(lines) + "\n")


def write_medium_zone_parameters() -> None:
	"""
	Ecrit les parametres Van Genuchten intermediaires dans E_zone_parameter.dat.

	Les valeurs proviennent du dictionnaire MEDIUM_ZONE_PARAMS defini en tete
	de script. Alpha et n sont a mi-chemin entre EASY_ZONE_PARAMS et les vraies
	valeurs geologiques, de facon a guider progressivement la convergence Picard.
	"""
	lines = []
	for zone_id, params in sorted(MEDIUM_ZONE_PARAMS.items()):
		ak, om, asp, ans, swres, alambda, cpm, rho = params
		lines.append(
			f"{zone_id}\t{ak:.2E}\t{om:.2E}\t{asp:.2E}\t{ans:.2E}"
			f"\t{swres:.2E}\t{alambda:.2E}\t{cpm:.0f}\t{rho:.2E}"
		)
	(CASE_DIR / "E_zone_parameter.dat").write_text("\n".join(lines) + "\n")


def clear_previous_outputs() -> None:
	"""
	Supprime quelques sorties de runs precedents.

	Objectif:
	- eviter de reutiliser par erreur un ancien S_pression_charge_temperature.dat,
	- forcer le script a travailler avec les fichiers generes pendant le run courant.
	"""
	for file_name in [
		"S_pression_charge_temperature.dat",
		"S_pression_t.dat",
		"S_temperature_t.dat",
		"S_saturation_t.dat",
		"S_flux_ZH_aq.dat",
	]:
		output_path = CASE_DIR / file_name
		if output_path.exists():
			output_path.unlink()


def set_steady_boundary_values() -> None:
	"""
	Ajuste E_cdt_aux_limites.dat pour le calcul permanent.

	Le helper `boundary_conditions_perm_2D_tdirect()` remet seulement `iclchgt=0`,
	mais ne fixe pas forcement les valeurs permanentes souhaitees pour ce cas.
	Ici, on impose donc explicitement les premieres charges lues dans:
	- E_chargeT_RD.dat
	- E_chargeT_RG.dat
	- E_chargeT_Riv.dat

	Cela permet de faire converger le permanent vers l'etat hydraulique
	correspondant au debut du forcage transitoire.
	"""
	charge_rd = read_first_scalar(CASE_DIR / "E_chargeT_RD.dat")
	charge_rg = read_first_scalar(CASE_DIR / "E_chargeT_RG.dat")
	charge_riv = read_first_scalar(CASE_DIR / "E_chargeT_Riv.dat")

	cdt_path = CASE_DIR / "E_cdt_aux_limites.dat"
	replace_setting(cdt_path, "valcl_gauche=", format_value(charge_rg))
	replace_setting(cdt_path, "valcl_droite=", format_value(charge_rd))


def configure_common_flags() -> None:
	"""
	Applique les drapeaux communs aux deux etapes du calcul.

	Pour AVA, on veut rester en hydraulique seule:
	- `ith=0` coupe le calcul thermique,
	- `rpth=0` evite tout transitoire thermique,
	- `itempi=0` indique qu'on n'utilise pas de temperature initiale variable,
	- `ichi2=1` indique qu'on utilise un champ de pression initiale variable,
	  lu dans E_pression_initiale.dat.
	"""
	replace_setting(CASE_DIR / "E_parametre.dat", "ith=", "0")
	replace_setting(CASE_DIR / "E_p_therm.dat", "rpth=", "0")
	replace_setting(CASE_DIR / "E_cdt_initiale.dat", "itempi=", "0")
	replace_setting(CASE_DIR / "E_cdt_initiale.dat", "ichi2=", "1")


def configure_saturated_steady_state() -> None:
	"""
	Prepare un permanent hydraulique sans zone non saturee.

	Ce premier run sert de preconditionnement numerique:
	- `rp=0` pour un calcul permanent,
	- `ivg=0` pour desactiver la ZNS,
	- `yunconfined=UNC` pour rester en mode captif/sature.
	"""
	configure_common_flags()
	replace_setting(CASE_DIR / "E_parametre.dat", "rp=", "0")
	replace_setting(CASE_DIR / "E_parametre.dat", "ivg=", "0")
	replace_setting(CASE_DIR / "E_parametre.dat", "yunconfined=", "UNC")
	boundary_conditions_perm_2D_tdirect()
	set_steady_boundary_values()


def configure_unsaturated_steady_state_easy() -> None:
	"""
	Prepare un permanent ZNS intermediaire avec des parametres Van Genuchten adoucis.

	C'est une etape de preconditionnement supplementaire entre le permanent sature
	et le permanent ZNS avec les vrais parametres :
	- `rp=0` calcul permanent,
	- `ivg=1` ZNS active,
	- `yunconfined=UNS` nappe libre,
	- E_zone_parameter.dat contient des alpha/n plus doux (EASY_ZONE_PARAMS).
	  Les vrais parametres sont sauvegardes dans E_zone_parameter.dat.orig.

	Cette etape est utile quand le saut direct du sature vers les vrais
	parametres ZNS bloque la convergence de Picard.
	"""
	configure_common_flags()
	replace_setting(CASE_DIR / "E_parametre.dat", "rp=", "0")
	replace_setting(CASE_DIR / "E_parametre.dat", "ivg=", "1")
	replace_setting(CASE_DIR / "E_parametre.dat", "yunconfined=", "UNS")
	boundary_conditions_perm_2D_tdirect()
	set_steady_boundary_values()
	write_easy_zone_parameters()


def configure_unsaturated_steady_state_medium() -> None:
	"""
	Prepare un permanent ZNS intermediaire avec des parametres Van Genuchten a mi-chemin.

	C'est la deuxieme etape de preconditionnement ZNS, entre EASY et les vrais parametres :
	- `rp=0` calcul permanent,
	- `ivg=1` ZNS active,
	- `yunconfined=UNS` nappe libre,
	- E_zone_parameter.dat contient MEDIUM_ZONE_PARAMS (alpha~4 pour zone 1, ~2.5 pour les autres).
	"""
	configure_common_flags()
	replace_setting(CASE_DIR / "E_parametre.dat", "rp=", "0")
	replace_setting(CASE_DIR / "E_parametre.dat", "ivg=", "1")
	replace_setting(CASE_DIR / "E_parametre.dat", "yunconfined=", "UNS")
	boundary_conditions_perm_2D_tdirect()
	set_steady_boundary_values()
	write_medium_zone_parameters()


def configure_unsaturated_steady_state() -> None:
	"""
	Prepare un permanent hydraulique avec zone non saturee et les vrais parametres.

	On repart du champ de pression fourni par le permanent ZNS intermediaire, puis:
	- `rp=0` garde un calcul permanent,
	- `ivg=1` reactive la courbe de retention ZNS,
	- `yunconfined=UNS` reactive le comportement nappe libre/ZNS,
	- E_zone_parameter.dat est restaure vers ses vrais parametres geologiques.
	"""
	configure_common_flags()
	replace_setting(CASE_DIR / "E_parametre.dat", "rp=", "0")
	replace_setting(CASE_DIR / "E_parametre.dat", "ivg=", "1")
	replace_setting(CASE_DIR / "E_parametre.dat", "yunconfined=", "UNS")
	boundary_conditions_perm_2D_tdirect()
	set_steady_boundary_values()


def configure_transient_state() -> None:
	"""
	Prepare les fichiers pour le run transitoire.

	- `rp=1` active le transitoire hydraulique,
	- `ivg=1` garde la zone non saturee active,
	- `yunconfined=UNS` conserve le meme mode hydraulique que le permanent ZNS,
	- `iclchgt=1` reactive les conditions aux limites variables dans le temps.
	"""
	configure_common_flags()
	replace_setting(CASE_DIR / "E_parametre.dat", "rp=", "1")
	replace_setting(CASE_DIR / "E_parametre.dat", "ivg=", "1")
	replace_setting(CASE_DIR / "E_parametre.dat", "yunconfined=", "UNS")
	boundary_conditions_2D_tdirect()


def run_ginette(step_name: str, timeout_sec: int) -> None:
	"""
	Execute l'executable `ginette` dans le dossier du cas.

	On passe par `cwd=CASE_DIR` parce que Ginette lit tous ses fichiers d'entree
	avec des chemins relatifs (`E_parametre.dat`, `E_p_therm.dat`, etc.).

	Comportement anti-blocage:
	- ecrit toute la sortie dans un log dedie,
	- affiche regulierement le temps ecoule,
	- arrete proprement le processus si le timeout est depasse.
	"""
	log_trace(f"Running Ginette in {step_name} mode...")
	log_path = CASE_DIR / f"ginette_{step_name.replace(' ', '_')}.log"
	start = time.time()
	next_progress = start + PROGRESS_EVERY_SEC

	with open(log_path, "w", encoding="utf-8") as log_file:
		proc = subprocess.Popen(
			["./ginette"],
			cwd=CASE_DIR,
			stdout=log_file,
			stderr=subprocess.STDOUT,
			text=True,
		)

		while proc.poll() is None:
			now = time.time()
			elapsed = int(now - start)
			if now >= next_progress:
				log_trace(f"[{step_name}] still running... {elapsed}s elapsed")
				next_progress = now + PROGRESS_EVERY_SEC

			if elapsed >= timeout_sec:
				proc.terminate()
				try:
					proc.wait(timeout=10)
				except subprocess.TimeoutExpired:
					proc.kill()
				log_trace(f"[{step_name}] timeout after {timeout_sec}s")
				raise TimeoutError(
					f"Ginette timeout during {step_name} after {timeout_sec}s. "
					f"See log: {log_path}"
				)

			time.sleep(1)

		if proc.returncode != 0:
			if log_path.exists() and log_path.stat().st_size == 0:
				log_trace(f"[{step_name}] ginette returned {proc.returncode} with empty stdout/stderr log")
			else:
				log_trace(f"[{step_name}] ginette returned {proc.returncode}, see {log_path}")
			raise RuntimeError(
				f"Ginette failed during {step_name} run with exit code {proc.returncode}. "
				f"See log: {log_path}"
			)

	log_trace(f"[{step_name}] completed successfully. Log: {log_path}")


def extract_last_pressure_snapshot() -> None:
	"""
	Construit un nouveau E_pression_initiale.dat a partir du permanent.

	Le fichier S_pression_charge_temperature.dat contient une ou plusieurs sorties
	successives de pression/charge/temperature.
	Chaque sortie contient `nmi` lignes, une par maille.
	
	Strategie:
	- lire tout le fichier,
	- recuperer le dernier bloc complet de taille `nmi`,
	- extraire la 4e colonne, qui contient la pression finale,
	- l'ecrire dans E_pression_initiale.dat.

	Ce fichier devient ensuite la condition initiale du run transitoire.
	"""
	source_path = CASE_DIR / "S_pression_charge_temperature.dat"
	if not source_path.exists() or source_path.stat().st_size == 0:
		raise FileNotFoundError(
			"S_pression_charge_temperature.dat was not generated by the steady-state run"
		)

	data = pd.read_csv(source_path, sep=r"\s+", header=None)
	if data.shape[1] < 4:
		raise ValueError("S_pression_charge_temperature.dat does not contain the pressure column")

	n_cells = read_nmi(CASE_DIR / "E_parametre.dat")
	if len(data) < n_cells:
		raise ValueError(
			f"S_pression_charge_temperature.dat has only {len(data)} rows, expected at least {n_cells}"
		)

	last_block = data.iloc[-n_cells:, 3]
	last_block.to_csv(CASE_DIR / "E_pression_initiale.dat", index=False, header=False)


def main() -> None:
	"""
	Orchestre la sequence complete:
	1.  compilation si besoin,
	2.  nettoyage des anciennes sorties critiques,
	3.  permanent sature (ivg=0, yunconfined=UNC),
	4.  extraction du champ de pression sature,
	5.  sauvegarde de E_zone_parameter.dat,
	6.  permanent ZNS EASY avec parametres VG tres adoucis (alpha~1, n~1.5),
	7.  extraction du champ de pression,
	8.  permanent ZNS MEDIUM avec parametres VG intermediaires (alpha~4, n~2),
	9.  extraction du champ de pression,
	10. restauration des vrais parametres VG,
	11. permanent ZNS avec les vrais parametres geologiques,
	12. extraction du champ de pression,
	13. reconfiguration en transitoire,
	14. run transitoire.
	"""
	# Le script peut etre lance depuis n'importe quel repertoire.
	# On se place donc explicitement dans le dossier du cas pour que:
	# - les helpers de src/src_python trouvent les fichiers E_*.dat,
	# - la compilation produise l'executable `ginette` au bon endroit,
	# - les sorties soient ecrites dans le dossier du cas AVA.
	previous_cwd = Path.cwd()
	os.chdir(CASE_DIR)
	if TRACE_LOG.exists():
		TRACE_LOG.unlink()
	log_trace(f"Starting run script from {previous_cwd}")
	log_trace(f"Case directory: {CASE_DIR}")

	try:
		# Compile `ginette` si l'executable n'existe pas encore dans le dossier du cas.
		log_trace("Compiling ginette if needed")
		compile_ginette_src(str(REPO_DIR))

		# On efface les sorties dont on depend pour la reprise afin d'eviter
		# toute confusion avec un ancien calcul.
		log_trace("Clearing previous outputs")
		clear_previous_outputs()

		# Etape 1: calcul permanent sature, sans ZNS.
		log_trace("Configuring saturated steady-state")
		configure_saturated_steady_state()
		run_ginette("steady-state-saturated", timeout_sec=STEADY_TIMEOUT_SEC)

		# Etape 2: on recopie la pression finale du premier permanent
		# vers le fichier d'initialisation du permanent ZNS intermediaire.
		log_trace("Extracting final pressure field from saturated steady-state")
		extract_last_pressure_snapshot()

		# Etape 3: on sauvegarde les vrais parametres de zone avant de les adoucir.
		log_trace("Backing up E_zone_parameter.dat")
		backup_zone_parameters()

		# Etape 4: calcul permanent ZNS intermediaire avec des parametres VG adoucis.
		# Ce preconditionneur aide la convergence de Picard avant les vrais parametres.
		log_trace("Configuring easy unsaturated steady-state (soft VG parameters)")
		configure_unsaturated_steady_state_easy()
		run_ginette("steady-state-unsaturated-easy", timeout_sec=STEADY_TIMEOUT_SEC)

		# Etape 5: on recopie le champ de pression du permanent EASY.
		log_trace("Extracting final pressure field from easy unsaturated steady-state")
		extract_last_pressure_snapshot()

		# Etape 6: calcul permanent ZNS MEDIUM avec parametres a mi-chemin.
		log_trace("Configuring medium unsaturated steady-state (intermediate VG parameters)")
		configure_unsaturated_steady_state_medium()
		run_ginette("steady-state-unsaturated-medium", timeout_sec=STEADY_TIMEOUT_SEC)

		# Etape 7: on recopie le champ de pression du permanent MEDIUM
		# et on restaure les vrais parametres de zone.
		log_trace("Extracting final pressure field from medium unsaturated steady-state")
		extract_last_pressure_snapshot()
		log_trace("Restoring original E_zone_parameter.dat")
		restore_zone_parameters()

		# Etape 8: calcul permanent ZNS avec les vrais parametres geologiques.
		log_trace("Configuring unsaturated steady-state (real VG parameters)")
		configure_unsaturated_steady_state()
		run_ginette("steady-state-unsaturated", timeout_sec=STEADY_TIMEOUT_SEC)

		# Etape 9: on recopie la pression finale du permanent ZNS reel
		# vers le fichier d'initialisation du transitoire.
		log_trace("Extracting final pressure field from unsaturated steady-state")
		extract_last_pressure_snapshot()

		# Etape 10: calcul transitoire en repartant de cette pression.
		log_trace("Configuring transient state")
		configure_transient_state()
		run_ginette("transient", timeout_sec=TRANSIENT_TIMEOUT_SEC)
	except Exception as exc:
		log_trace(f"Script failed: {exc}")
		raise
	finally:
		# On restaure le repertoire initial pour ne pas surprendre l'utilisateur
		# si le script est appele depuis un autre emplacement.
		os.chdir(previous_cwd)
		log_trace(f"Returned to initial directory: {previous_cwd}")

	log_trace("Sequence complete: saturated perm -> ZNS easy perm -> ZNS medium perm -> ZNS real perm -> transient.")
	log_trace("Temperature stays disabled for the AVA case (ith=0, itempi=0, rpth=0).")


if __name__ == "__main__":
	# Point d'entree standard du script.
	main()
