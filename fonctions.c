#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" // besoin des donnees du probleme
#include "methodesnum.h" // besoin de la methode de quadrature de gauss-legendre


/************************
 fonctions conditions 
 aux limites et exactes
*************************/

double T_ex(double x, double y) {
    /* fonction solution exacte 
    x: reel dont on veut calculer l'image */
    return cosh(M_PI*y)*cos(M_PI*x);
}

double f_0(double x, int m, double alpha) {
    /* fonction f_0 condition limite sur Gamma_0 
    x: reel dont on veut calculer l'image */
    return cos(M_PI*x);
}

double f_3_ex(double x) {
    /* fonction f_3 exacte condition limite sur Gamma_3 
    x: reel dont on veut calculer l'image */

    //on recupere d'abord les donnees du probleme
    double H;
    H = recup_H(H);
    // puis on retourne la fonction souhaitee
    return cosh(M_PI*H)*cos(M_PI*x);
}

double q_0(double x) {
    /* fonction q_0 exacte condition limite sur Gamma_0 */
    return 0.;
}

/************************
 fonctions a approcher
*************************/

double inte_h(double x, int m, double alpha) {
    /* fonction a integrer dans l'expression de h */
    //on recupere d'abord les donnees du probleme
    double L,H;
    L = recup_L(L);
    H = recup_H(H);
    // puis on retourne la fonction souhaitee
    return cos(m*M_PI*x/L);
}

double h(double x){
    /* fonction qui approche (par une somme de M termes) la fonction h
    x: reel dont on veut calculer l'image */
    double S;
    int i;
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    S = -q_0(x) + (1/(L*H))*gauss(n,f_0,0,L,0,0);
    for (i = 1; i<=M; i++) {
        //S = S + (2*M_PI/pow(L,2)) * (i*cos(i*M_PI*x/L)/sinh(i*M_PI*x/L)) * f_0(x,i,0.) * cosh(i*M_PI*H/L) * gauss(n,inte_h,0,L,i,0);
        S = S + (2*M_PI/pow(L,2)) * ((i*cos(i*M_PI*x/L))/sinh(i*M_PI*x/L))*f_0(x, i, 0.)*cosh(i*M_PI*H/L) * gauss(n,inte_h,0.,L,i,0.);
    }
    return S;
}

double inte_f3(double x, int m, double alpha) {
    /* fonction que l'on integre dans l'expression de f_3 */
    //on recupere les donnees du probleme
    float L;
    L = recup_L(L);
    return cos(m*M_PI*x/L)*h(x);
}

double f_3(double x, double alpha){
    /* fonction qui approche f_3 
    x: reel dont on veut calculer l'image par f_3
    alpha: parametre de Lavrentier */
    int i;
    double S;
    //on recupere les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);

    S = (1./alpha)*h(x) - (1./alpha)*(H/(alpha*H + 1.));
    printf("%f\n",S);
    for (i=1; i <= M; i++) {
        S = S - (1./alpha)*((pow(L,2)*sinh((i*M_PI*H)/L)/(i*M_PI+alpha*pow(L,2)*sinh((i*M_PI*H)/L)))*cos((i*M_PI*x)/L)*gauss(n,inte_f3,0,L,i,alpha));
        printf("%f\n",S);
    }
    return S;
}

/*************************
 fonctions a integrer pour
 les coefficients
**************************/

double B0_f3_f0(double x, int m, double alpha){
    /* fonction a integer pour le coefficient B_0
    car en argument de gauss() il faut une fonction
    x: reel dont on veut calculer l'image */
    m = 0; // ici on a le terme 0 de la somme
    return f_3(x,alpha)- f_0(x,m,alpha);
}

double Am_f3_f0(double x, int m, double alpha){
    /* fonction a integer pour le coefficient A_m
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule A_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (f_3(x,alpha) - exp((-m*M_PI*H)/L)*f_0(x,m,alpha))*cos((m*M_PI*x)/L);
}

double Bm_f3_f0(double x, int m, double alpha) {
    /* fonction a integer pour le coefficient B_m
    x: reel dont on veut calculer l'image
    m: indice m pour lequel on calcule B_m */
    
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (exp((m*M_PI*H)/L)*f_3(x,alpha) - f_0(x,m,alpha))*cos((m*M_PI*x)/L);
}

/************************
 fonctions coefficients
*************************/

double A_0() {
    /* fonction qui calcule le coefficient A_0 pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (1/L)*gauss(n,f_0,0,L,0,0);
}

double B_0(double alpha) {
    /* fonction qui calcule le coefficient B_0 pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (1/(L*H))*gauss(n,B0_f3_f0,0,L,0,alpha);
}

double A_m(int m, double alpha) {
    /* fonction qui calcule le coefficient A_m pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    float L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (1/L*sinh((m*M_PI*H)/L))*gauss(n,Am_f3_f0,0,L,m,alpha);
}

double B_m(int m, double alpha) {
    /* fonction qui calcule le coefficient B_m pour l'approximation de T_tilde */
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    return (1/L*sinh((m*M_PI*H)/L))*gauss(n,Bm_f3_f0,0,L,m,alpha);
}

/************************
 fonction T tilde approchee
*************************/

double T_tilde(double x, double y, double alpha) {
    /* fonction qui donne l'approximation de T_tilde */
    double T;
    int i;
    //on recupere d'abord les donnees du probleme
    double L,H;
    int M,n;
    L = recup_L(L);
    H = recup_H(H);
    M = recup_M(M);
    n = recup_n(n);
    // puis on retourne la fonction souhaitee
    T = A_0()+ B_0(alpha)*y;
    for (i = 1; i<=M; i++) {
        T = T + (A_m(i,alpha)*exp(i*M_PI*y/L) + B_m(i,alpha)*exp(-i*M_PI*y/L))*cos(i*M_PI*x/L);
    }
    return T;
}

double fctn_test(double x, double y, double alpha) {
    return 3.*y+x-1.;
}

double derivee_test(double x, double y, double alpha) {
    return 3.;
}