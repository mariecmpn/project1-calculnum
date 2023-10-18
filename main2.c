#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"
#include "fonctions2.c"

int main() {

    /***************************
     definition des variables 
    ****************************/ 

    int i, j; // entiers pour les boucles for
    int N = 50; // nombre de points de maillage en x
    double pas = 1./(N-1); // pas de maillage en x
    double *Points = malloc(sizeof(int[N])); // tableau qui contient les valeurs des x_i
    double eps = 1.E-6; // tolerance epsilon pour la methode de Newton

    // recuperation des donnees du probleme definies dans donnees.c
    double L = recup_L(L), H = recup_H(H);
    int M = recup_M(M), n = recup_n(n);

    double Alpha[22]; // tableau qui contient les differentes valeurs de alpha

    double *Gamma_app = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_alpha calculees numeriquement
    double *Gamma_exact = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_ex calculees numeriquement

    /***************************
     remplissage des tableaux
    ****************************/ 

    // remplissage de Points
    double x_i = 0.;
    for (i = 0; i<N; i++) {
        Points[i] = x_i;
        x_i  = x_i+pas;
    }

    // remplissage de Gamma_exact
    for (i = 0; i<N; i++) {
        Gamma_exact[i] = Gamma_ex(Points[i]);
    }

    // remplissage de Alpha
    Alpha[0] = 10.;
    Alpha[1] = 1.;
    for (i = 2; i<22; i++) {
        Alpha[i] = 1./pow(10.,i+1);
    }

    /***************************
     Calcul de la frontiere libre
     pour plusieurs alpha
    ****************************/ 

   



    // on retourne 0
    return 0;
}