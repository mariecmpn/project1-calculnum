#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"


int main() {
    // definition des variables
    double L, H;
    int M, n;

    // recuperation des donnees du probleme definies dans donnees.c
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);


    // test f_3
    /*double x = 0.1;
    double y = 0.5;
    float app, ex;

    app = f_3(x, 1.E-10);
    ex = f_3_ex(x);

    printf("%s%f\n", "Valeur approchee ", app);
    printf("%s%f\n", "Valeur exacte ", ex);*/

    // test newton

    double err;
    double newt = newton(0., fctn_test, derivee_test, 0.00001, 0.);

    err = fabs(newt-1./3.);
    printf("%f", err);

    /*double alpha = 1.E-1;
    //y = (1./alpha)*h(x);
    y = h(x);
    printf("%f\n", y);
    y = (1./alpha)*(H/(alpha*H + 1.));
    printf("%f\n", y);*/


    // on retourne 0
    return 0;
}