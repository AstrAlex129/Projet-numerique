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
    S = np.zeros(len(rho))
    F = np.zeros(len(rho))
    L = np.zeros(len(rho))
    
    # paramètres élémentaires
    # dxi = xi[1]-xi[0]
    # dt = t[1] -t[0]
    
    dxi = xi[-1]/len(xi)
    dt = t[-1]/len(t) # le temps est réellement une liste ? à initialiser ou pas ?

    # Obtention de la source selon le fichier de référence fourni


    # Obtention des termes tel que L = F + S 
    S[1:-1] = xi[1:-1]*rho[1:-1]*((R[2:]-R[:-2])/(2*dt)) - (2*R0*rho[1:-1]*v[1:-1])/(xi[1:-1]*R[1:-1]*t0) - (rho[1:-1]*R0)/(t0*R[1:-1])*((v[2:]-v[:-2])/(2*dxi))
    F[1:-1] = -(v[1:-1]*R0)/(t0*R[1:-1]) * (rho[2:]-rho[:-2])/(2*dxi)
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L

def F_spatiotemp_euler(rho : list, v : list, R : list, xi : list, t : list, rho0 : float, v0 : float, R0 : float, U0, t0 : float, K : float, gamma : float) -> list : 
    """ Fonction totale de l'équation d'euler, variation temporelle de v """
    
    # Initialisation des tableaux
    S = np.zeros(len(v))
    F = np.zeros(len(v))
    L = np.zeros(len(v))
    
    # paramètres élémentaires

    # dxi = xi[1]-xi[0]
    # dt = t[1] -t[0]

    dxi = xi[-1]/len(xi)
    dt = t[-1]/len(t) # même questionnement que la fonction de continuité

    # Idem que la fonction de continuité
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*dxi)) - (rho[1:-1]*rho0)/R[1:-1] ((U0[2:] - U0[:-2])/(2*dxi)) + (rho[1:-1]*rho0*R0*xi[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt)) # définition de U0 problématique
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*dxi))
    
    L[1:-1] = (F[1:-1] + S[1:-1])*(rho[1:-1]*rho0*(R0/t0))
    
    return L

def F_geometrique(v) : 
    return v[-1] # R[t+dt] = R[t] + v[-1]*dt

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

    _ , v , _ = Y
    return F_geometrique(v)

def F_global_poisson(t, Y):
    pass

#=========================================
# Initialisation Runge-Kutta 3
#=========================================

def RK4_spatiotemp_base(f,Y,t,h): # pour chaque grandeur (faire manuellement)
    """ Patern RK4 1 unité général """

    # Ici Y : list[list]

    # Evolution simultanée

    # ordre 1
    k1_rho = f[0](t,Y)
    k1_v = f[1](t,Y)
    k1_R = f[2](t,Y)

    # ordre 2
    Y_temp_1 = [Y[0] + (h*k1_rho)/2, Y[1] + (h*k1_v)/2, Y[2] + (h*k1_R)/2]

    k2_rho = f[0](t+h/2,Y_temp_1)
    k2_v = f[1](t+h/2,Y_temp_1)
    k2_R = f[2](t+h/2,Y_temp_1)

    # ordre 3
    Y_temp_2 = [Y_temp_1[0] + (h*k2_rho)/2, Y_temp_1[1] + (h*k2_v)/2, Y_temp_1[2] + (h*k2_R)/2]

    k3_rho = f[0](t+h/2,Y_temp_2)
    k3_v = f[1](t+h/2,Y_temp_2)
    k3_R = f[2](t+h/2,Y_temp_2)

    # ordre 4
    Y_temp_3 = [Y_temp_2[0] + (h*k3_rho), Y_temp_2[1] + (h*k3_v), Y_temp_2[2] + (h*k3_R)]

    k4_rho = f[0](t+h,Y_temp_3)
    k4_v = f[1](t+h,Y_temp_3)
    k4_R = f[2](t+h,Y_temp_3) 
    
    Y_next_rho = Y[0] + (h/6)*(k1_rho +2*k2_rho+ 2*k3_rho + k4_rho)
    Y_next_v = Y[1] + (h/6)*(k1_v +2*k2_v+ 2*k3_v + k4_v)
    Y_next_R = Y[2] + (h/6)*(k1_R +2*k2_R+ 2*k3_R + k4_R)

    Y_next = [Y_next_rho, Y_next_v, Y_next_R]
    t_next = t+h

    return Y_next, t_next # revoir l'implémentation de Y et l'utilisation future de Y_next comportant trois Y. Dois-je utiliser chaque Y pour chaque grandeur ou le tout pour chaque ?

def RK4_spatiotemp_boucle(f, Y_init, h): # pour une grandeur 
    """ Patern RK4 boucle général """

    Y_approx_l = [] # liste simple

    t = 0
    
    Y=Y_init.copy()
    
    for _ in range(3):
        Y_approx_l.append([Y[0].copy(),Y[1].copy(),Y[2].copy()]) # ajouter chaque composante de variable individuellement mais pour "un seul" élément
        Y,t = RK4_spatiotemp_base(f,Y,t,h)
    
    return Y_approx_l, t 
    
#=========================================
# Adam-Bashforth ordre 3
#=========================================   
    
def AB3_spatiotemp_base(f, Y_init , h : float, N : int): # Pour une grandeur 
    """ Patern AB3 général, permet d'obtenir les valeurs d'une grandeur particulière via la fonction choisi et pas toutes """
        
    Y_rk, t = RK4_spatiotemp_boucle(f,Y_init,h)
    Y0, Y1, Y2 = Y_rk
    t0, t1, t2 = 0*h, 1*h, 2*h
    
    # Obtention des primaires

    #F0
    F0_rho = f[0](t0, Y0)
    F0_v = f[1](t0, Y0)
    F0_R = f[2](t0, Y0)

    #F1
    F1_rho = f[0](t1, Y1)
    F1_v = f[1](t1, Y1)
    F1_R = f[2](t1, Y1)

    #F2
    F2_rho = f[0](t2, Y2)
    F2_v = f[1](t2, Y2)
    F2_R = f[2](t2, Y2)

     
    
    # Démarrage à n = 2, et paramètres 
    Y = [Y2[0].copy(), Y2[1].copy(), Y2[2].copy()]
    t = 2*h
    
    res = [Y0,Y1,Y2]
    
    # récupération des valeurs 
    for _ in range(2,N):
        Y_next_rho = Y[0] + h*(23*F2_rho - 16*F1_rho + 5*F0_rho)/12
        Y_next_v   = Y[1] + h*(23*F2_v   - 16*F1_v   + 5*F0_v)/12
        Y_next_R   = Y[2] + h*(23*F2_R   - 16*F1_R   + 5*F0_R)/12

        Y_next = [Y_next_rho, Y_next_v, Y_next_R]

        res.append([Y_next_rho.copy(), Y_next_v.copy(), Y_next_R.copy()]) # .copy() uniquement des listes pas des tableaux
        t_next = t + h

        # mise à jour pour le prochain pas
        F0_rho, F0_v, F0_R = F1_rho, F1_v, F1_R
        F1_rho, F1_v, F1_R= F2_rho, F2_v, F2_R 
        F2_rho, F2_v, F2_R= f[0](t_next, Y_next), f[1](t_next, Y_next), f[2](t_next, Y_next)
        Y, t = Y_next, t_next
    
    return res

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
    U0 = v0**2
    p0 = K*rho0**gamma
    R0 = 10e3 # m
    t0 = R0/v0 # s
    xi0 = 1/R0 # à vérifier (m^-1) # initialiser correctement

    # Implémentation de sigma, tau, w et r dans la suite, ainsi que dans les fonctions précédentes
    # Relecture de le document de référence

    # variables à l'équilibre
    rho_c = 1.9891e18 # kg/m-3
    r0 = np.sqrt((np.pi*K)/(2*G))
    rho0_eq = rho_c*np.sinc(np.pi*(r/r0))
    U0_eq = -2*K*rho_c*(1+np.sinc(np.pi*((R*xi)/r0))) # s2 # utilisation de variable qui ne sont pas encore défini
    w_eq = 0 

    # Nombre d'itération
    N = 1000
    t = np.linspace(1,100,N) # s

    xi = np.linspace(0,1,N) # discrétisation à revoir
    
    # Variables d'algorithme 
    h = 1 # durée donc cohérence avec l'oscillation des étoiles minutes, heures... ? 
    Y = np.array([rho0,v0,R0]) # RK4
    f = [F_global_continuite, F_global_euler, F_global_geometrique]

    # Obtention des grandeurs à étudier 

    res_grandeur = AB3_spatiotemp_base(f,Y,h,N)

    rho_evolution = [Yg[0] for Yg in res_grandeur] 
    v_evolution = [Yg[1] for Yg in res_grandeur] 
    R_evolution = [Yg[2] for Yg in res_grandeur] 






    
    
    
    
    
    
    
    
    
    
    