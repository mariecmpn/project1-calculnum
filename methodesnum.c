#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "methodesnum.h"



float gauss(int n, float (*f)(float,int,float), float a, float b, int m, float alpha) {
    /* fonction pour calculer une integrale numeriquement par la methode de Gauss-Legendre
    n: nombre de points de quadrature: 1, 2 ou 3
    f: la fonction dont on souhaite approcher l'integrale
    a: borne inf du domaine sur lequel on integre
    b: borne sup du domaine sur lequel on integre
    retourne: la valeur approchee de l'integrale
    m: pour les sommes (on recupere l'indice m)
    alpha: parametre de Lavrentier qu'on recupere pour calculer f_3 */

    int i;
    float Inte;
    float *T = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les points de quadrature
    float *Tp = malloc(sizeof(int[n])); // allocation dynamique du tableau T qui contient les poids de gauss

    if (n!=1 && n!=2 && n!=3) { 
        printf("%s", "Nombre de points de quadrature non defini (n = 1, 2 ou 3)");
    }
    else {
        switch (n) { // on remplit les tableaux pour les poids et points de quadrature en fonction de n
            case 1:
                T[0] = 0.;
                Tp[0] = 2.;
                break;
            case 2:
                T[0] = -1./sqrt(3.);
                T[1] = 1./sqrt(3.);
                Tp[0] = 1.;
                Tp[1] = 1.;
                break;
            case 3:
                T[0] = -sqrt(3.)/sqrt(5.);
                T[1] = 0.;
                T[2] = sqrt(3.)/sqrt(5.);
                Tp[0] = 5./9.;
                Tp[1] = 8./9.;
                Tp[2] = 5./9.;
        }
        Inte = 0;
        for (i = 0; i < n; i++) { // on calcule notre integrale
            Inte = Inte + ((b-a)/2.)*(Tp[i] * f(((b-a)*T[i]/2.)+(a+b)/2.,m,alpha));
        }
    }

    // on desalloue la place memoire des tableaux
    free(T);
    free(Tp);

    return Inte; // on retourne la valeur de l'integrale
}
