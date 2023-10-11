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


    // test f_3
    float x = 0.5;
    //float y = 0.5;
    /*float app, ex;

    app = f_3(x, 1.E-1);
    ex = f_3_ex(x);

    printf("%s%f\n", "Valeur approchee ", app);
    printf("%s%f\n", "Valeur exacte ", ex);*/

    x = gauss(n,inte_h,0,L,1,0);
    printf("%f", x);

    // on retourne 0
    return 0;
}