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

def F_spatiotemp_continuite(rho : list, v : list, R : list, epsilon : list, t : list, rho0 : float, v0 : float, R0 : float, t0 : float) -> list : 
    """ Fonction totale de l'équation de continuité """
    
    # Initialisation des tableaux
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    # paramètres élémentaires
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = epsilon[1:-1]*rho[1:-1]*((R[2:]-R[:-2])/(2*dt)) - (2*R0*rho[1:-1]*v[1:-1])/(epsilon[1:-1]*R[1:-1]*t0) - (rho[1:-1]*R0)/(t0*R[1:-1])*((v[2:]-v[:-2])/(2*depsilon))
    F[1:-1] = -(v[1:-1]*R0)/(t0*R[1:-1]) * (rho[2:]-rho[:-2])/(2*depsilon)
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L

def F_spatiotemp_euler(rho : list, v : list, R : list, epsilon : list, U : list, t : list, rho0 : float, v0 : float, R0 : float, U0 : float, t0 : float, K : float, gamma : float) -> list : 
    """ Fonction totale de l'équation d'euler """
    
    # Initialisation des tableaux
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    # paramètres élémentaires
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*depsilon)) - (rho[1:-1]*rho0)/R[1:-1] * ((U[2:] - U[:-2])/(2*depsilon)) + (rho[1:-1]*rho0*R0*epsilon[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt))
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*depsilon))
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L
    
    
def F_spatiotemp_poisson(rho : list, epsilon : list, U : list, rho0 : float, U0 : float):
    """ Fonction totale de l'équation de poisson """
    pass
    

#=========================================
# Fonction Globale
#=========================================

def F_global_euler(t, Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R, epsilon, U = Y
    return F_spatiotemp_euler(rho, v, R, epsilon, U, t, h, rho0, v0, R0, U0, t0, K, gamma)

def F_global_continuite(t, Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R, epsilon, U = Y 
    return F_spatiotemp_continuite(rho, v, R, epsilon, t, rho0, v0, R0, t0)

def F_global_poisson(t, Y):
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R, epsilon, U = Y
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
    
def AB3_spatiotemp_base(f, Y_init :list[list] , h : float, N : int):
    """ Patern AB3 général """
        
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
    for n in range(2,N):
        Y_next = Y + h*(23*F2 - 16*F1 + 5*F0)/12
        t_next = t + h

        res.append(Y_next.copy())

        # mise à jour pour le prochain pas
        F0, F1, F2 = F1, F2, f(t_next, Y_next)
        Y, t = Y_next, t_next
    
    return res 

#=========================================
# Utilisation du polynôme Tchebychev
#=========================================
    

if __name__ == "__main__" :
    
    # Détermination de grandeur modulable
    
    # Grandeurs
    c = 3e8 #m.s-1
    K = 1 # à changer plus tard
    gamma = 2
    G = 6.67e-11 # m3.kg^-1.s^-2
    
    # variables initiales
    v0 = c
    rho0 = 1.66e17 # kg.m^-3
    p0 = K*rho0**gamma
    R0 = 10 # km
    t0 = R0/v0
    U0 = v0**2 # même unité 
    epsilon0 = 1/R0 # à vérifier
    
    # Nombre d'itération
    N = 500
    t = np.linspace(1,N)
    
    # Variables amovibles 
    rho = np.random.uniform(N)
    v = np.random.uniform(N)
    p = np.random.uniform(N)
    R = np.random.uniform(N)
    U = np.random.uniform(N)
    epsilon = np.random.uniform(N)
    
    # Variables d'algorithme 
    h = 1 # à changer, RK4
    Y = np.array([rho0,v0,R0,epsilon0,U0]) # RK4
    
    # Déroulement opérationnel
    pass
    
    
    
    
    
    
    
    
    
    
    