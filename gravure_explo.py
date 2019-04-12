# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 09:57:53 2017

@author: godef
"""
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
import matplotlib as mpl

#-----Fonctions_------------------------------
def Balayage(ij, i, spin, n, H0) :
    # On boucle
    for j in range(ij, N + 1, 2) :
        # On calcul E
        E = Energie(i, j, spin, n, H0)
        # On test la probabilité
        if np.random.binomial(1, Metropolis(E, i, j, n)) == 1  :
            # On met à jour le réseau
            spin[i, j] *= -1
    return spin
# Fin du balayage

def Energie(i, j, spin, n, H0) :
    # On somme l'énergie sur les spin des voisins
#    sumSpin = reverse*spin[i, j]*spin[i + dx[:], j + dy[:]].sum()
    sumSpin = 0
    for k in [-1,1]: #boucle sur les voisins
        sumSpin -= J * spin[i, j] * ((spin[i+k, j])+(spin[i, j+k]))
    
    sumSpin -= J * spin[i, j] * ((spin[i-1, j-1])+(spin[i-1, j+1]) + (spin[i+1, j-1])+(spin[i+1, j+1]))
    # On retourne l'énergie
    return sumSpin - H0*spin[i, j]

def Metropolis(E, i, j, t) :
    return np.min((1, np.exp(2*E/(k*temperature(i, j, t)))))
#---------------------------------------------

#---Function Ising Modelisation---------------
def Ising(spin, H0) :
    magne = np.zeros(nIter)
    # Boucle temporelle
    for n in range(0, nIter) :
        print('Temperature : ' + str(T), 'Iteration : ' + str(n))
        # On balaye les noeud "blanc" du reseau
        for i in range(1, N + 1, 2) :
            # On alterne le noeud de départ
            ij = (i % 2) + 1
            Balayage(ij, i, spin, n, H0)
            
        # On balaye les noeud "noir" du reseau
        for i in range(1, N + 1, 2) :
            # On alterne le noeud de départ
            ij = ((i + 1) % 2) + 1
            Balayage(ij, i, spin, n, H0)
        # On balaye les noeud "blanc" du reseau
        for i in range(2, N + 1, 2) :
            # On alterne le noeud de départ
            ij = (i % 2) + 1
            Balayage(ij, i, spin, n, H0)
            
        # On balaye les noeud "noir" du reseau
        for i in range(2, N + 1, 2) :
            # On alterne le noeud de départ
            ij = ((i + 1) % 2) + 1
            Balayage(ij, i, spin, n, H0)
            
        # Condition périodique sur les spins
        spin[0, :] = spin[-1, :]
        spin[:, 0] = spin[:, -1]
        spin[0, 0] = spin[-1, -1]
        spin[0, -1] = spin[-1, 0]

        # On calcul la magnétisation
        # On calcul l'énergie moyenne
        magne[n] = spin.sum() /(N**2)
    
    return spin, magne
    
def spaceProfil(i, j):
    if (i - 32)**2 + (j - 32)**2 <= R**2 :
        return 1
    else:
        return 0
        
def temperature(i, j, t) :
    return T0 + deltaT*((i - i0)**2 + (j - i0)**2 <= R**2)*np.exp(-((t-t0)**2)/sigma**2)

def MapSpin(spin, T) :
    plt.figure(facecolor='white', figsize=(12,8), dpi=72)
    grid=spin[1:-1,1:-1]
    color = mpl.colors.ListedColormap(['k','w'])
    bounds=[-2,0,2]
    norme = mpl.colors.BoundaryNorm(bounds,color.N)
    plt.imshow(grid,interpolation='nearest',cmap=color,norm=norme)
    plt.title('Tempterature ' + str(T))
    plt.savefig('gravure_exploration_T-{}.png'.format(T))
#---------------------------------------------

#-----Constantes------------------------------
# Taille du réseau
N = 64
# Nombre d'itérations temporelle
nIter = 500
# Constante d'interraction spin-spin
J = 1/3
# Constante de perméabilité du vide
mu0 = constants.mu_0
# Constance pour Metropolis (Boltzmann)
k = 1
# Température ambiante
T = 1.5
# Stencil de voisinage
dx = np.array([-1, 0, 1, 0])
dy = np.array([0, -1, 0, 1])
# Tableau de l'énergie moyenne
energie = np.zeros(nIter)
# Tableau de magnétisation moyenne
magn = np.zeros(nIter)
# Représente la position i du centre du faisceau
i0 = 32
# Représente la position j du centre du faisceau
j0 = 32
# Rayon de la région chauffée
R = 10
# Hausse maximale de temprérature
deltaT = 0.1
# Iteration à intensité maximale
t0 = 100
# La durée du pulse
sigma = 10
# Champ magnétique externe
H0 = -0.2
#---------------------------------------------

#---Ising spin network initialisation---------
# On initialise le reseau de spin à +1 partout
spin = np.ones([N + 2, N + 2])

# Conditions limites
spin[0, :] = spin[-1, :]
spin[:, 0] = spin[:, -1]
spin[0, 0] = spin[-1, -1]
spin[0, -1] = spin[-1, 0]

# On conserve le modèle initial de spin
spinInit = np.array(spin)
#---------------------------------------------

#---Animation Ising Model----------------------
space = np.arange(0, nIter)
deltaT = 15
T0 = 0.05

s, m = Ising(spin, H0)

MapSpin(s, T)