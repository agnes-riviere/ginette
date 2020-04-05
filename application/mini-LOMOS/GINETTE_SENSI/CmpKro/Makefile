PROGRAMME = $(NAME)$(VERSION)

NAME=CmpKro

SHELL = /bin/bash
VERSION=0.01

SOURCES.c = CmpKro.c\

SOURCES.l =

SOURCES.y =

SOURCES.h =

SOURCES.f =


#CC=g++
#CC=/usr/local/gcc-4.7.2/bin/g++ #for numerix which has gcc3.3 by default
#CC=/usr/local/gcc-4.7.2/bin/gcc #for numerix which has gcc3.3 by default
#CC=${MPI_DIR}/bin/mpicc
CC=gcc
#test2=/usr/fin

#VERSION of GCC is needed to compile
GCC_VER=$(shell $(CC) -v 2>&1|tail -1|awk '{print $$3}'|awk -F "." '{print $$1$$2$$3}')
#FLAG_CC=$(shell echo $$(CC)|awk -F "/" '{for(i=2;i<NF+1;i++) print $$i}')
#FLAG_CC=$(shell echo $$test2|awk -F "/" '{for(i=2;i<NF+1;i++) print $$i}')
#FLAG_CC=$(shell echo $$test2)

COMPILO=$(shell echo $(CC) | awk -F "/" '{print $$NF}')

#Configure dependencies
#libprint and libts
#Versions
V_LP=1.24
V_TS=1.69
#Paths for Include
INCL_TS=-I$(LIB_HYDROSYSTEM_PATH)/libts/trunk/src
INCL_LP=-I$(LIB_HYDROSYSTEM_PATH)/libprint/trunk/src

LIB_PATH=  -L$(LIB_HYDROSYSTEM_PATH)
LIBS=   -lts$(V_TS)_$(COMPILO) -lprint$(V_LP)_$(COMPILO) -lm -lc
# Attention, l'ordre est important!

OPTCOMP = -g
#OPTW= -Wall
OPTW= -Wno-write-strings #-fpermissive
OPTOMP= -fopenmp
OPTD = -DOMP $(OPTOMP)

CFLAGS_GLOBAUX= $(OPTCOMP) $(OPTW)
CFLAGS= $(CFLAGS_GLOBAUX) $(OPTD)


INCLUDES = -I/usr/local/gcc-4.7.2/include  -I/usr/local/gcc-4.7.2/lib/gcc/x86_64-unknown-linux-gnu/4.7.2/include/#for system using a native gcc older than 4.x. For those systems, an installation of another gcc in /usr/local is needed, as well as the definition of the paths here

#CPP : Preprocessor, gere les inclure.
#C'est ici qu'il faut indiquer les chemins d'acces aux .h
CPPFLAGS= -I$(INCLUDES) $(INCL_LP) $(INCL_TS)

#CPP : Preprocessor, gere les inclure.
#C'est ici qu'il faut indiquer les chemins d'acces aux .h
YACC= bison -y
#LEX= /share/apps/flex/bin/flex
LEX= flex
YFLAGS= -dv
LFLAGS= -v


# Variables déduites des précédentes:
OBJ1=$(SOURCES.y:.y=.o) $(SOURCES.l:.l=.o)
OBJ2=$(SOURCES.c:.c=.o)
OBJ=  $(OBJ1) $(OBJ2)

%.o : %.c
	@echo 'Creating $(@)...'
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $(basename $@).c -o $@

$(PROGRAMME) :  $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LIB_PATH) -o $(PROGRAMME)  $(OBJ)  $(LIBS)


#all :  $(OBJ)
#	$(LINK.c) -o $(PROGRAMME)  $(OBJ) $(LIBS)

all : clean $(PROGRAMME)


install : all
	rm ~/executables/$(NAME)
	ln -s ~/Programmes/$(NAME)/trunk/src/$(NAME)$(VERSION) ~/executables/$(NAME)
	@echo "You can start" $(NAME)$(VERSION) "from everywhere with the command >" $(NAME) "*.COMM *.log"


# Pour faire un lint des programmes:
lint: $(LINTFILES)
	$(LINT.c) $(LINTFILES)

# Pour nettoyer quand on change de machine:
clean:
	rm -f $(OBJ) $(LINTFILES) $(PROGRAMME) core y.* apoub*


# Suppression du $(RM) par défaut :
.y.o :
	$(YACC.y) $<
	$(COMPILE.c)  -o $@ y.tab.c


print :
	@echo "Settings for the Compilation of version $(PROGRAMME) "
	 @echo "$(CC) version : $(GCC_VER)"
	@echo "Compiler extension : $(COMPILO)"
	 @echo "INCLUDES sat to $(INCLUDES)"
	 @echo "compiling options $(CFLAGS)"
	 @echo "with inclusion of  $(CPPFLAGS)"
	@echo "PETSc folder is $(PETSC_FOLD) for the inclusion of rules and variables files"
	@echo "and static libraries $(LIBS)"
