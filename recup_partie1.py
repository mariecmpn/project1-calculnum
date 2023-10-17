# recuperation des donnees pour le projet 1 - calcul numerique scientifique
# Partie 1 -resolution du probleme de Cauchy

#%% importation packages

import numpy as np
import matplotlib.pyplot as plt
import os

#%% definition des variables

L = 1.
H = 1.

Alpha = np.zeros(22)

Alpha[0] = 10.
Alpha[1] = 1.
for i in range(2,22):
    Alpha[i] = 10**(-(i-1))
    
N = 50 #nombre de points de maillage pour x
h = L/(N-1) #pas de maillage

X = np.zeros(N)
x = 0.
for i in range(N):
    X[i] = x
    x = x+h
  
#%% execution du programme partie1

#os.system('make clean')
#os.system('make partie1')
os.system('./partie1') # on execute le fichier executable partie1 depuis notre console python
# penser a le compiler d'abord
    

#%% lecture des fichiers

# solution exacte sur Gamma_3
file = open('exact_1.txt', 'r')

data = file.read()
Exact = data.split()
Ex = [float(i) for i in Exact]

#print(Ex)


# solution approchee sur Gamma_3
file = open('approche_1.txt', 'r')

data = file.read()
Approche = data.split()
Approche = [float(i) for i in Approche]
App = np.array(Approche)

#print(App)

# On divise le tableau en 20 tableaux (1 tableau pour chaque alpha)
T_alpha = np.split(App, 22)

#%% calcul des erreurs

Err = np.zeros(22)
for i in range(22):
    Err[i] = np.linalg.norm(T_alpha[i] - Ex)



#%% graphiques

plt.loglog(Alpha, Err)
plt.title('Erreur entre les solutions exacte et approchée sur Gamma_3')
plt.show()


plt.plot(X, Ex)
plt.title('Solution exacte T_ex sur Gamma_3')
plt.show()

plt.plot(X, T_alpha[0])
plt.title('Solution approchee T_app sur Gamma_3')
plt.show()


#%% 

#def f_3(x):
    #return np.cosh(np.pi*H)*np.cos(np.pi*x)

#plt.plot(X, f_3(X))
#plt.show()


#%% lecture des fichiers pour un alpha donne = 1.

file = open('solapp_omega.txt', 'r')
data = file.read()
SolApp = data.split()
SolApp = [float(i) for i in SolApp]
SolApp = np.array(SolApp)
SolApp = np.split(SolApp, 50)
