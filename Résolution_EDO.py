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
    """ Fonction totale de l'équation de continuité, variation temporelle de rho """
    
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
    """ Fonction totale de l'équation d'euler, variation temporelle de v """
    
    # Initialisation des tableaux
    S = np.zeros(len(t))
    F = np.zeros(len(t))
    L = np.zeros(len(t))
    
    # paramètres élémentaires
    depsilon = epsilon[1]-epsilon[0]
    dt = t[1] -t[0]
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*depsilon)) - (rho[1:-1]*rho0)/R[1:-1] * U0 * ((U[2:] - U[:-2])/(2*depsilon)) + (rho[1:-1]*rho0*R0*epsilon[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt))
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*depsilon))
    
    L[1:-1] = (F[1:-1] + S[1:-1])*(rho[1:-1]*rho0*(R0/t0))
    
    return L
    
    
def F_spatiotemp_poisson(rho : list, R : list,  epsilon : list, U : list, rho0 : float, U0 : float): # Problématique sur la mise en forme
    """ Fonction totale de l'équation de poisson, variation spatiale de U """
    depsilon = epsilon[1]-epsilon[0]

    Up = np.zeros(len(t))
    F = np.zeros(len(t))

    Up[1:-1] = (U[2:] - U[:-2])/(2*depsilon)

    Um = np.array([Up[1:-1], U[1:-1]])

    # du/depsilon = au + b = f(u,epsilon,t,F) tq F = 4np.pi*G*rho*(rho0/u0)*R^2
    A = np.array([[-2/epsilon[1:-1], 0],[1,0]])
    B = np.array([[4*np.pi*G*rho[1:-1]*(rho0/U0)*R[1:-1]**2]])

    F[1:-1] = A@Um + B

    return F  
    

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

def F_global_poisson(t, Y): # comme il n'y a pas de dépendance de t directement, doit-on mettre epsilon ou laisser t ?
    """ Fonction intermédiaire pour calculer la fonction voulue """
    
    rho, v, R, epsilon, U = Y
    return F_spatiotemp_poisson(rho,R, epsilon, U, rho0, U0)

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
    for n in range(2,N):
        Y_next = Y + h*(23*F2 - 16*F1 + 5*F0)/12
        t_next = t + h

        res.append(Y_next.copy())

        # mise à jour pour le prochain pas
        F0, F1, F2 = F1, F2, f(t_next, Y_next)
        Y, t = Y_next, t_next
    
    return res[-1]

#=========================================
# Fonction décriptive de la résolution
#=========================================

def solveur(Y,t,h):

    # Obtention de l'évolution des paramètres tout au long de l'expérience

    # Paramètre évolutifs rho, v, R, U - pour moi R est relié à rho
    # Les paramètres sont calculés indépendamment alors que chaque variable dépend des autres, comment faire
    
    rho_star = np.zeros(N)
    v_star = np.zeros(N)
    R_star = np.zeros(N)
    U_star = np.zeros(N)
    Epsilon_star = np.zeros(N) # problématique dans l'évolution, constante ?

    # je ne comprends pas l'utilité de ceci
    F_euler = np.zeros(N)
    F_continuite = np.zeros(N)
    F_poisson = np.zeros(N)

    rho0, v0, R0, epsilon0, R0 = Y

    rho_star[0] = rho0
    v_star[0] = v0
    R_star[0] = R0
    U_star[0] = U0
    Epsilon_star[0] = epsilon0

    grandeur = [rho_star, v_star, U_star]
    L_fonc = [0,1,2]
    for f,n,g in zip(F,L_fonc,grandeur) :
        Y_rk,_= RK4_spatiotemp_boucle(f,Y_init,h)
        _, Y1, Y2 = Y_rk
        # a finir
    

    # remplir les 3 premières valeurs des tableaux, epsilon constant ? 

    for n in range(2,N-1):

        rho_star[n+1] = AB3_spatiotemp_base(F_global_continuite,Y_init, h,n)
        v_star[n+1] = AB3_spatiotemp_base(F_global_euler,Y_init, h, n)
        R_star[n+1] = (M_ns_moy/((4/3)*np.pi*rho_star))**(1/3)
        U_star[n+1] = AB3_spatiotemp_base(F_global_poisson,Y_init,h,n)

        Y_next = [rho_star[n+1], v_star[n+1], R_star[n+1], Epsilon_star[n+1], U_star[n+1]]

        Y_init = Y_next

    return rho_star, v_star, R_star, U_star, Epsilon_star

#=========================================
# Utilisation du polynôme Tchebychev (à voir)
#=========================================
    

if __name__ == "__main__" :
    
    # Détermination de grandeur modulable
    
    # Grandeurs
    c = 3e8 # m.s-1
    gamma = 2 # 1
    K = 4.25e-3 # m5.kg^-1.s^-2
    G = 6.67e-11 # m3.kg^-1.s^-2
    Ms = 1.98e30 # kg
    M_ns_moy = 1.4*Ms # kg
    
    # variables initiales
    v0 = c
    rho0 = 1.66e17 # kg.m^-3
    p0 = K*rho0**gamma
    R0 = 10e3 # m
    t0 = R0/v0 # s
    U0 = v0**2 # s2
    epsilon0 = 1/R0 # à vérifier (m^-1)
    
    # Nombre d'itération
    N = 1000
    t = np.linspace(1,100,N) # s
    
    # Variables d'algorithme 
    h = 1 # durée donc cohérence avec l'oscillation des étoiles minutes, heures... ? 
    Y = np.array([rho0,v0,R0,epsilon0,U0]) # RK4

    F = [F_global_continuite, F_global_euler, F_global_poisson]




    
    
    
    
    
    
    
    
    
    
    