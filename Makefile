PROG = ginette

RELEASE = v2

EXECUTABLE = $(PROG)

F90 = gfortran
F90FLAGS = -cpp
GIN_F90_DIR = src/ginette_V2.f90

GCC = gcovr
GCChtmlFLAGS = --html --html-details
GCCcompileFLAGS = -O0 --coverage

# Peut-être modifié pour être plus pratique
TARGET_PATH = application/2017_AVAV_SENSI

# target

compile : $(GIN_F90_DIR)
	$(F90) $(F90FLAGS) $(GIN_F90_DIR) -o $(TARGET_PATH)/$(EXECUTABLE)

compile_debug :
	$(F90) $(F90FLAGS) -DDEBUG $(GIN_F90_DIR) -o $(TARGET_PATH)/$(EXECUTABLE)

run : 
	cd $(TARGET_PATH) && ./$(EXECUTABLE)

clean :
	cd $(TARGET_PATH) && rm -f $(EXECUTABLE) S_* Sim* *.gcno *.gcda *.html *.css

sup_run : 
	$(MAKE) clean 
	$(MAKE) compile
	$(MAKE) run

debug :
	$(MAKE) clean
	$(MAKE) compile_debug
	$(MAKE) run


GCC : $(GIN_F90_DIR)
	$(MAKE) clean
	$(F90) $(F90FLAGS) $(GCCcompileFLAGS) $(GIN_F90_DIR) -o $(TARGET_PATH)/$(EXECUTABLE)
	$(MAKE) run
	$(GCC) --root src --object-directory . $(GCChtmlFLAGS) -o $(TARGET_PATH)/$(EXECUTABLE).html
	firefox $(TARGET_PATH)/$(EXECUTABLE).html

#Spécialement pour l'application Dharrma
DHARRMA_PATH = application/model_dharrma
DHARRMA_VENV = $(DHARRMA_PATH)/venv
DHARRMA_PYTHON = venv/bin/python3

# Crée le venv et installe les paquets nécessaires (build + exécution) seulement s'il n'existe pas déjà
$(DHARRMA_VENV)/bin/python3 :
	python3 -m venv $(DHARRMA_VENV)
	$(DHARRMA_VENV)/bin/pip install --upgrade pip
	$(DHARRMA_VENV)/bin/pip install "Cython==3.3.0" "numpy==2.2.6" "pandas==2.3.3" "matplotlib==3.10.9" "pgcore==1.5.5" "pygimli==1.5.5"

init_dharrma : $(DHARRMA_VENV)/bin/python3
	cd $(DHARRMA_PATH) && $(DHARRMA_PYTHON) setup.py build_ext --inplace && touch lib/__init__.py

run_dharrma : $(DHARRMA_VENV)/bin/python3
	cd $(DHARRMA_PATH) && $(DHARRMA_PYTHON) main_DHARRMA.py

sup_run_dharrma : $(DHARRMA_VENV)/bin/python3
	cd $(DHARRMA_PATH)/input_ginette && rm -f $(EXECUTABLE) S_* Sim* *.gcno *.gcda *.html *.css
	cd $(DHARRMA_PATH) && $(DHARRMA_PYTHON) main_DHARRMA.py


