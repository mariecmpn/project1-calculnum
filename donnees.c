#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "donnees.h"


float recup_L(float L) {
    /* fonction permettant de recuperer la valeur de la longueur (horizontale) du domaine */
    return 1.;
}

float recup_H(float H) {
    /* fonction permettant de recuperer la valeur de la longueur (verticale) du domaine */
    return 1.;
}

int recup_M(float M) {
    /* fonction permettant de recuperer le nombre de termes pour les sommes infinies que l'on a tronque */
    return 100;
}

int recup_n(float n) {
    /* fonction permettant de recuperer le nombre de points de quadrature utilises pour la methode de Gauss-Legendre */
    return 2;
}