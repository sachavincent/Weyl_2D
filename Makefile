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
OBJECTS=limace.o utils.o Stereo.o Track.o FindBestMatch.o

# --> ci-dessous, ajoutez a la suite les noms des fichiers executables 
# separes par un espace
EXE=Stereo Track FindBestMatch

# Ne modifiez pas cette partie
.PHONY: all
all: $(EXE)

# Dependances non implicites des executables
OBJ=limace.o utils.o
stereo: $(OBJ) Stereo.o
track: $(OBJ) Track.o
findbestmatch: $(OBJ) FindBestMatch.o
# --> ajoutez ici une ligne par programme selon le meme modele


# Dependances non implicites des objets
HEAD=limace.h utils.h
limace.o: $(HEAD)
utils.o: $(HEAD)
stereo.o: $(HEAD)
track.o: $(HEAD)
findbestmatch.o: $(HEAD)
# --> ajoutez ici une ligne par programme selon le meme modele


# Ne modifiez pas cette partie
clean:
	$(RM) $(OBJECTS) $(EXE)
