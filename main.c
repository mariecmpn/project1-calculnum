#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonctions.h"
#include "donnees.h" 
#include "methodesnum.h"


int main() {
    /***************************
     definition des variables 
    ****************************/ 

    int i, j; // entier pour les boucles for
    int N = 50; // nombre de points de maillage en x
    double h = 1./(N-1); // pas de maillage en x
    double *Points = malloc(sizeof(int[N])); // tableau qui contient les valeurs des x_i

    double Alpha[20]; // tableau qui contient les differentes valeurs de alpha

    // recuperation des donnees du probleme definies dans donnees.c
    double L = recup_L(L), H = recup_H(H);
    int M = recup_M(M), n = recup_n(n);

    double *T_app = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_alpha calculees numeriquement
    double *T_exact = malloc(sizeof(int[N])); // tableau qui contient les images des x_i par T_ex calculees numeriquement

    double y = H; // on se place sur Gamma_3

    /***************************
     remplissage des tableaux
    ****************************/ 

    // remplissage de Points
    double x_i = 0.;
    for (i = 0; i<N; i++) {
        Points[i] = x_i;
        x_i  = x_i+h;
    }

    // remplissage de Alpha
    for (i = 0; i<20; i++) {
        Alpha[i] = 1./pow(10.,i);
    }



    /***************************
     question 4
    ****************************/ 


    // test f_3
    /*double x = 0.5;
    float app, ex;

    app = T_tilde(x, y, 1.E-10);
    ex = T_ex(x,y);

    printf("%s%f\n", "Valeur approchee ", app);
    printf("%s%f\n", "Valeur exacte ", ex);*/


    /*double alpha = 1.E-1;
    //y = (1./alpha)*h(x);
    y = h(x);
    printf("%f\n", y);
    y = (1./alpha)*(H/(alpha*H + 1.));
    printf("%f\n", y);*/

    y = H; // on se place sur Gamma_3
    FILE *approche;
    approche = fopen("approche_1.txt", "w");
    for (i = 0; i < 20; i++) {
        for (j = 0; j < N; j++) {
            T_app[j] = T_tilde(Points[j], y, Alpha[i]);
            fprintf(approche, "%f", T_app[j]);
            fputs(" ", approche);
        }
        fputs("\n", approche);
    }

    /***************************
    enregistrement dans un 
    fichier texte de T_ex
    ****************************/ 

    for (i = 0; i<N; i++) {
        T_exact[i] = T_ex(Points[i], y); // calcul de T sur Gamma_3
    }

    FILE *exact;
    exact = fopen("exact_1.txt", "w"); // on ouvre le fichier

    for (i = 0; i<N; i++) { // on ecrit dans le fichier
        fprintf(exact, "%f", T_exact[i]);
        fputs(" ", exact);
    }

    fclose(exact); // on ferme le fichier


    // on desalloue l'espace memoire des tableaux alloues dynamiquement
    free(Points);
    free(T_exact);
    free(T_app);

    // on retourne 0
    return 0;
}