#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 15:51:06 2026

@author: skraburskia
"""

# Bibliothèque
import numpy as np 
import matplotlib.pyplot as plt

#=========================================
# Equation Spatio-temporelle
#=========================================

#=========================================
# Fonction initiale
#=========================================

def F_spatiotemp_continuite(rho : list, v : list, R : list, xi : list, t : list, rho0 : float, v0 : float, R0 : float, t0 : float) -> list : 
    """ Fonction totale de l'équation de continuité, variation temporelle de rho """
    
    # Initialisation des tableaux
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    # paramètres élémentaires
    # dxi = xi[1]-xi[0]
    # dt = t[1] -t[0]
    
    dxi = xi[-1]/len(xi)
    dt = t[-1]/len(t)

    # Obtention des termes tel que L = F + S 
    S[1:-1] = xi[1:-1]*rho[1:-1]*((R[2:]-R[:-2])/(2*dt)) - (2*R0*rho[1:-1]*v[1:-1])/(xi[1:-1]*R[1:-1]*t0) - (rho[1:-1]*R0)/(t0*R[1:-1])*((v[2:]-v[:-2])/(2*dxi))
    F[1:-1] = -(v[1:-1]*R0)/(t0*R[1:-1]) * (rho[2:]-rho[:-2])/(2*dxi)
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L

def F_spatiotemp_euler(rho : list, v : list, R : list, xi : list, t : list, rho0 : float, v0 : float, R0 : float, U0, t0 : float, K : float, gamma : float) -> list : 
    """ Fonction totale de l'équation d'euler, variation temporelle de v """
    
    # Initialisation des tableaux
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    # paramètres élémentaires

    # dxi = xi[1]-xi[0]
    # dt = t[1] -t[0]

    dxi = xi[-1]/len(xi)
    dt = t[-1]/len(t)
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*dxi)) - (rho[1:-1]*rho0)/R[1:-1] ((U0[2:] - U0[:-2])/(2*dxi)) + (rho[1:-1]*rho0*R0*xi[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt)) # définition de U0 problématique
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*dxi))
    
    L[1:-1] = (F[1:-1] + S[1:-1])*(rho[1:-1]*rho0*(R0/t0))
    
    return L

def F_geometrique(rho,t,rho0) : # R doit évoluer donc on doit lui appliquer une fonction spécifique ? 
    F = np.zeros(N)

    dt = t[-1]/len(t)

    F[1:-1] = (M_ns_moy/((4/3)*np.pi*rho0))**(1/3)*(-3/rho[1:-1]**4)*((rho[2:]-rho[:-2])/(2*dt))

    return F

def F_spatiotemp_poisson(): # Faire plus tard
    pass

#=========================================
# Fonction Globale
#=========================================

def F_global_euler(t, Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R = Y
    return F_spatiotemp_euler(rho, v, R, xi, t, h, rho0, v0, R0, U0, t0, K, gamma)

def F_global_continuite(t, Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R = Y 
    return F_spatiotemp_continuite(rho, v, R, xi, t, rho0, v0, R0, t0)

def F_global_geometrique(t,Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """

    rho, _ , _ = Y
    return F_geometrique(rho,t,rho0)

def F_global_poisson(t, Y):
    pass

#=========================================
# Initialisation Runge-Kutta 3
#=========================================

def RK4_spatiotemp_base(f,Y,t,h):
    """ Patern RK4 1 unité général """
    
    k1 = f(t,Y)
    k2 = f(t+h/2,Y+(h*k1)/2)
    k3 = f(t+h/2,Y +(h*k2)/2)
    k4 = f(t+h,Y + h*k3)
    
    Y_next = Y + (h/6)*(k1 +2*k2+ 2*k3 + k4)
    t_next = t+h
    return Y_next, t_next

def RK4_spatiotemp_boucle(f, Y_init, h):
    """ Patern RK4 boucle général """
    
    Y_approx = []
    t = 0
    
    Y=Y_init.copy()
    
    for _ in range(3):
        Y_approx.append(Y.copy())
        Y,t = RK4_spatiotemp_base(f,Y,t,h)
    
    return Y_approx, t 
    
#=========================================
# Adam-Bashforth ordre 3
#=========================================   
    
def AB3_spatiotemp_base(f, Y_init , h : float, N : int):
    """ Patern AB3 général, permet d'obtenir les valeurs d'une grandeur particulière via la fonction choisi et pas toutes """
        
    Y_rk, t = RK4_spatiotemp_boucle(f,Y_init,h)
    Y0, Y1, Y2 = Y_rk
    t0, t1, t2 = 0*h, 1*h, 2*h
    
    # Obtention des primaires
    F0, F1, F2 = f(t0,Y0), f(t1,Y1), f(t2,Y2)
    
    # Démarrage à n = 2, et paramètres 
    Y = Y2.copy()
    t = 2*h
    
    res = [Y0,Y1,Y2]
    
    # récupération des valeurs 
    for _ in range(2,N):
        Y_next = Y + h*(23*F2 - 16*F1 + 5*F0)/12
        t_next = t + h

        res.append(Y_next.copy())

        # mise à jour pour le prochain pas
        F0, F1, F2 = F1, F2, f(t_next, Y_next)
        Y, t = Y_next, t_next
    
    return res[-1] # ne renvois que la dernière valeur

#=========================================
# Fonction décriptive de la résolution
#=========================================

def solveur(Y,t,h):
    # Initialisation des variables à évoluer 
    rho_tab = np.zeros(N)
    v_tab = np.zeros(N)
    R_tab = np.zeros(N)

    # Initialisation des dictionnaires pour l'étude couplée des grandeurs
    Ys = {"rho" : rho_tab,
         "v" : v_tab,
         "R" : R_tab}

    Ys[all][0] = Y

    pass 

#=========================================
# Utilisation du polynôme Tchebychev (à voir plus tard)
#=========================================
    

if __name__ == "__main__" :
    
    # Détermination de grandeur modulable
    
    # Grandeurs
    c = 3e8 # m.s-1
    gamma = 2 # 1
    K = 4.25e-3 # m5.kg^-1.s^-2
    G = 6.67e-11 # m3.kg^-1.s^-2
    Ms = 1.98e30 # kg
    M_ns_moy = 1.4*Ms # kg Que l'on suppose constante
    
    # variables initiales
    v0 = c
    rho0 = 1.66e17 # kg.m^-3
    p0 = K*rho0**gamma
    R0 = 10e3 # m
    t0 = R0/v0 # s

    rho_c = 1.9891e18 # kg/m-3
    r0 = np.sqrt((np.pi*K)/(2*G))
    U0 = -2*K*rho_c*(1+np.sinc(np.pi*((R*xi)/r0))) # s2 # utilisation de variable qui ne sont pas encore défini
    xi0 = 1/R0 # à vérifier (m^-1) # initialiser correctement
    
    # Nombre d'itération
    N = 1000
    t = np.linspace(1,100,N) # s

    xi = np.linspace(0,1,N) # discrétisation à revoir
    
    # Variables d'algorithme 
    h = 1 # durée donc cohérence avec l'oscillation des étoiles minutes, heures... ? 
    Y = np.array([rho0,v0,R0,U0]) # RK4





    
    
    
    
    
    
    
    
    
    
    