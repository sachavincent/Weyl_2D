# L3 Informatique - Informatique graphique, traitement et analyse d'images
# Bases du traitement et de l'analyse d'images
#
# Dans ce qui suit, remplacez monprog par le nom de votre programme
#
# Pour lancer les compilations, tapez :
# make
# Pour effacer tous les fichiers objets et executables, tapez :
# make clean
#
# Pour prendre en compte un nouveau programme, ajoutez-le a tous endroits designes par -->
# ci-dessous.

# Ne modifiez pas cette partie 
# Exception : en cas de refus de l'option -Wpedantic, la remplacer par -pedantic
CC=gcc
CFLAGS=-Wall -Wextra -std=c99 -pedantic -fopenmp
LDFLAGS= $(CFLAGS) 
LDLIBS=-lm 
RM=rm -f

# --> ci-dessous, ajoutez a la suite les noms des fichiers objets 
# separes par un espace
OBJECTS=limace.o main.o Stereo.o Suivi.o Utils.o

# --> ci-dessous, ajoutez a la suite les noms des fichiers executables 
# separes par un espace
EXE=main Stereo Suivi 

# Ne modifiez pas cette partie
.PHONY: all
all: $(EXE)

# Dependances non implicites des executables
OBJ=limace.o
main: $(OBJ) main.o
stereo: $(OBJ) Stereo.o
suivi: $(OBJ) Suivi.o
utils: $(OBJ) Utils.o
# --> ajoutez ici une ligne par programme selon le meme modele


# Dependances non implicites des objets
HEAD=limace.h
limace.o: $(HEAD)
main.o: $(HEAD)
stereo.o: $(HEAD)
suivi.o: $(HEAD)
utils.o: $(HEAD)
# --> ajoutez ici une ligne par programme selon le meme modele


# Ne modifiez pas cette partie
clean:
	$(RM) $(OBJECTS) $(EXE)
