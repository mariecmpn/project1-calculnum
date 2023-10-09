# Compiler used
CC = gcc-11

# Regles pour la partie 1 du projet
partie1: main.o methodesnum.o donnees.o fonctions.o
	$(CC) -o partie1 main.o methodesnum.o donnees.o fonctions.o

methodenum.o: methodesnum.c
	$(CC) -c methodesnum.c

donnees.o: donnees.c
	$(CC) -c donnees.c

fonctions.o: fonctions.c donnees.h methodesnum.h
	$(CC) -c fonctions.c

main.o: main.c fonctions.c donnees.h methodesnum.h fonctions.h
	$(CC) -c main.c

# Regles pour la partie 2 du projet

# Regle clean
clean: 
	rm *.o partie1