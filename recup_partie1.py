# recuperation des donnees pour le projet 1 - calcul numerique scientifique
# Partie 1 -resolution du probleme de Cauchy

#%% importation packages

import numpy as np
import matplotlib.pyplot as plt

#%% definition des variables

L = 1.
H = 1.

Alpha = np.zeros(20)

for i in range(20):
    Alpha[i] = 10**(-(i+1))
    
N = 50 #nombre de points de maillage pour x
h = L/(N-1) #pas de maillage

X = np.zeros(N)
x = 0.
for i in range(N):
    X[i] = x
    x = x+h

#%% lecture des fichiers

file = open('exact_1.txt', 'r')

data = file.read()
Exact = data.split()
Ex = [float(i) for i in Exact]

print(Ex)



#%% graphiques

plt.plot(X, Ex)
plt.title('Solution exacte T_ex sur Gamma_3')
plt.show()



