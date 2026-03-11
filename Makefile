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
TARGET_PATH = application/2017_AVA_SENSI

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
init_dharrma : 
	cd $(DHARRMA_PATH) && python3 setup.py build_ext --inplace && touch lib/__init__.py

run_dharrma :
	cd $(DHARRMA_PATH) && python3 main_DHARRMA.py

sup_run_dharrma :
	cd $(DHARRMA_PATH)/input_ginette && rm -f $(EXECUTABLE) S_* Sim* *.gcno *.gcda *.html *.css
	cd $(DHARRMA_PATH) && python3 main_DHARRMA.py

