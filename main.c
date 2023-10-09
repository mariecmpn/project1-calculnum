#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"


int main() {
    // definition des variables
    float L, H;
    int M, n;

    // recuperation des donnees du probleme definies dans donnees.c
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);


    // test methode num gauss pour une fonction dont on connait l'integrale
    float err;
    err = fabs(gauss(n,fctn_test,-1,1,0,0)-8./3.);
    printf("%s%f", "Erreur = ", err);

    // on retourne 0
    return 0;
}