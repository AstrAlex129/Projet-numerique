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
# Fonctions initiales
#=========================================

def F_spatiotemp_continuite(rho, v, w , R, xi, t): 
    """ 
    Implémente la fonction totale de l'équation de continuité suite à la transition physique -> numérique.
    Permet de déterminer l'évolution de rho.

    Args : 
        rho (list) : Densité de l'étoile à neutron.
        v (list) : Vitesse d'oscillation (ou du déplacement interne).
        w (list) : Paramètre numérique v(r,t)/r
        R (list) : Rayon évolutif physique de l'étoile à neutron.
        xi (list) : [Définition] représente r/R(t).
        t (list) : Temps d'étude.
    
    Returns :
        list : (L) - L est l'expression numérique de la fonction destinées
    """    
    # Initialisation des tableaux
    S = np.zeros(len(rho)) # terme source de la fonction - initialisation 
    F = np.zeros(len(rho)) # terme principal de la fonction - initialisation
    L = np.zeros(len(rho)) # terme total souhaité - initialisation
   
    ### [A faire : Modification de fonction via w, tau et sigma] ###

    # Obtention des termes tel que L = F + S 
    S[1:-1] = xi[1:-1]*rho[1:-1]*((R[2:]-R[:-2])/(2*dt)) - (2*R0*rho[1:-1]*v[1:-1])/(xi[1:-1]*R[1:-1]*t0) - (rho[1:-1]*R0)/(t0*R[1:-1])*((v[2:]-v[:-2])/(2*dxi))
    F[1:-1] = -(v[1:-1]*R0)/(t0*R[1:-1]) * (rho[2:]-rho[:-2])/(2*dxi)
    
    L[1:-1] = F[1:-1] + S[1:-1]
    
    return L

def F_spatiotemp_euler(rho, v, w, R, xi, t): 
    """ 
    Implémente la fonction totale de l'équation d'Euler suite à la transition physique -> numérique.
    Permet de déterminer l'évolution de v.
    
    Args : 
        rho (list) : Densité de l'étoile à neutron.
        v (list) : Vitesse d'oscillation (ou du déplacement interne).
        w (list) : Paramètre numérique v(r,t)/r
        R (list) : Rayon évolutif physique de l'étoile à neutron.
        xi (list) : [Définition] représente r/R(t).
        t (list) : Temps d'étude.
        
    Returns :
        list : (L) - L est l'expression numérique de la fonction destinées
    """     
    # Initialisation des tableaux
    S = np.zeros(len(v))    # terme source de la fonction - initialisation
    F = np.zeros(len(v))    # terme principal de la fonction - initialisation
    L = np.zeros(len(v))    # terme finale souhaité - initialisation

    ### [A faire : Modification de fonction via w, tau et sigma] ###
    
    # Obtention des termes tel que L = F + S 
    S[1:-1] = -(gamma*rho[1:-1]*K*rho0**(gamma-1))/(rho[1:-1]*R[1:-1]) * ((rho[2:] - rho[:-2])/(2*dxi)) - (rho[1:-1]*rho0)/R[1:-1] ((U0_eq[2:] - U0_eq[:-2])/(2*dxi)) + (rho[1:-1]*rho0*R0*xi[1:-1]*v[1:-1])/t0 * ((R[2:] - R[:-2])/(2*dt)) 
    F[1:-1] = -(rho[1:-1]*rho0*R0**2*v[1:-1])/(t0**2*R[1:-1])*((v[2:] - v[:-2])/(2*dxi))
    
    L[1:-1] = (F[1:-1] + S[1:-1])*(rho[1:-1]*rho0*(R0/t0))
    
    return L

def F_geometrique(v_inst, R_prec) : # correspond à la dernière valeur de v tel que v va être un float
    """
    Implémente la fonction géométrique nécessaire pour l'obtention de R(t).

    Args :
        v_inst (float) : Vitesse instantannée 
        R_prec (float) : Rayon instantannée de l'étoile à neutron à l'instant précédent

    Returns :
        float (R_inst) - R_inst est le rayon instantannée obtenu à une vitesse donnée à l'instant t
    """
    R_inst = R_prec + v_inst*dt
    return R_inst

def F_spatiotemp_poisson(rho, R, U): # Faire plus tard. Pour l'instant U = U_eq
    """
    Implémente la fonction totale de Poisson suite à la transition physique -> numérique.
    Permet d'obtenir l'évolution de U.

    Args :
        rho (list) : Densité de l'étoile à neutron.
        R (list) : Rayon physique de l'étoile à neutron.
        U (list) : Potentiel associé
     
    Returns :
        list : (L) - L est l'expression numérique de la fonction destinées
    """
    pass

#=========================================
# Fonction Globale
#=========================================

def F_global_euler(t, Y):
    """ 
    Fonction intermédiaire pour calculer la fonction globale d'Euler.

    Args :
        t (list) : Temps de simulation
        Y(list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné 
    
    Returns :
        function : (F_spatiotemp_euler) - Permet de calculer la fonction global et d'obtenir son return 
    """
    rho, v, R = Y
    return F_spatiotemp_euler(rho, v, R, xi, t)

def F_global_continuite(t, Y):
    """ 
    Fonction intermédiaire pour calculer la fonction globale de continuité.
    
    Args :
        t (list) : Temps de simulation.
        Y(list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné.
        
    Returns :
        function : (F_spatiotemp_continuite) - Permet de calculer la fonction global et d'obtenir son return 
    """
    rho, v, R = Y 
    return F_spatiotemp_continuite(rho, v, R, xi, t)

def F_global_geometrique(t,Y,R_prec):
    """ 
    Fonction intermédiaire pour calculer la fonction globale géométrique.
        
    Args :
        t (list) : Temps de simulation
        Y(list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné
        R_prec (float) : Rayon instantannée du 'palier' précédent 
            
    Returns :
        function : (F_geometrique) - Permet de calculer la fonction global et d'obtenir son return 
    """
    _ , v , _ = Y
    return F_geometrique(v,R_prec)

def F_global_poisson(t, Y):
    """ 
    Fonction intermédiaire pour calculer la fonction globale géométrique.
            
    Args :
        t (list) : Temps de simulation
        Y(list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné
                
    Returns :
        func : (F_spatiotemp_poisson) - Permet de calculer la fonction global et d'obtenir son return 
    """
    pass

#=========================================
# Initialisation Runge-Kutta 3
#=========================================

def RK4_spatiotemp_base(f,Y,t,h): # pour chaque grandeur (faire manuellement)
    """ 
    Evolution des différents paramètres via un algorithme RK4 d'une fonction vectorisée.
    Permet d'obtenir la valeur suivante, de chaque grandeur.

    Args :
        f (func) : Fonction associant la valeur étudiée souhaitée
        Y (list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné
        t (list) : Temps utilisé pour la simulation 
        h (float) : Pas de temps 
    
    Returns :
        tuples : (Y_next, t_next) - Y_next est le tableau des CI pour l'étape d'après de chaque grandeur, idem pour t_next 
    """
    ## Evolution simultanée
    # Expression manuelle de chaque constante selon chaque ordre.
    # Permet d'obtenir une étude simultannée via les équations couplées.
    
    # Etapes principales pour le RK4
    # Ordre 1
    k1_rho = f[0](t,Y)
    k1_v = f[1](t,Y)
    k1_R = f[2](t,Y)

    # Ordre 2
    Y_temp_1 = [Y[0] + (h*k1_rho)/2, Y[1] + (h*k1_v)/2, Y[2] + (h*k1_R)/2] # Mise au point nécessaire pour chaque ordre supérieur

    k2_rho = f[0](t+h/2,Y_temp_1)
    k2_v = f[1](t+h/2,Y_temp_1)
    k2_R = f[2](t+h/2,Y_temp_1)

    # Ordre 3
    Y_temp_2 = [Y_temp_1[0] + (h*k2_rho)/2, Y_temp_1[1] + (h*k2_v)/2, Y_temp_1[2] + (h*k2_R)/2]

    k3_rho = f[0](t+h/2,Y_temp_2)
    k3_v = f[1](t+h/2,Y_temp_2)
    k3_R = f[2](t+h/2,Y_temp_2)

    # Ordre 4
    Y_temp_3 = [Y_temp_2[0] + (h*k3_rho), Y_temp_2[1] + (h*k3_v), Y_temp_2[2] + (h*k3_R)]

    k4_rho = f[0](t+h,Y_temp_3)
    k4_v = f[1](t+h,Y_temp_3)
    k4_R = f[2](t+h,Y_temp_3) 

    # Obtention de l'évolution couplé de chaque variable selon l'évolution
    Y_next_rho = Y[0] + (h/6)*(k1_rho +2*k2_rho+ 2*k3_rho + k4_rho)
    Y_next_v = Y[1] + (h/6)*(k1_v +2*k2_v+ 2*k3_v + k4_v)
    Y_next_R = Y[2] + (h/6)*(k1_R +2*k2_R+ 2*k3_R + k4_R)

    Y_next = [Y_next_rho, Y_next_v, Y_next_R]
    t_next = t+h

    return Y_next, t_next ## revoir l'implémentation de Y et l'utilisation future de Y_next comportant trois Y. Dois-je utiliser chaque Y pour chaque grandeur ou le tout pour chaque ?

def RK4_spatiotemp_boucle(f, Y_init, h): # pour une grandeur 
    """ 
    Evolution des différents paramètres via un algorithme RK4 d'une fonction vectorisée.
    Permet d'obtenir les 3 premières valeurs, de chaque grandeur, selon le paramètre souhaité 
    
    Args :
        f (func) : Fonction associant la valeur étudiée souhaitée
        Y (list[list]) : Tableau contenant les valeurs d'études - initiales ou à un instant donné
        h (float) : Pas de temps 
        
    Returns :
        tuples : (Y_approx_l, t) - Y_approx_l est le tableau globale de chaque grandeur évolutive, t est le temps associé à la simulation 
    """
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
    
def AB3_spatiotemp_base(f, Y_init , h, N): # Pour une grandeur 
    """
    Evolution de chaque grandeur et obtention des valeurs souhaitées selon la fonction ou grandeur de référence (pour obtenir sa propre évolution)

    Args :
        f (func) : Fonction de référence.
        Y_init (list) : Liste des paramètres initiaux de référence.
        h (float) : Pas de temps.
        N (int) : Nombre d'itération de simulation

    Returns :
        list[list] : (res) - res est l'ensemble des valeurs de grandeurs obtenu pour l'étude d'une valeur de référence
    """  
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
    
    # Récupération des valeurs | Problématique sur l'utilisation du temps (ou pas)
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

#=========================================
# Mise en place des sphères représentatives
#=========================================

def sphere_numerique(R_fixe,theta,phi):
    """
    Représentation d'une sphère fixe numérique (sans le plot ou avec, à voir) via les grandeurs utilisées.
    Permet de mettre en évidence l'étude futur de l'oscillation.

    Args :
        R_fixe (float) : Rayon de la sphère fixe.
        theta (list) : Paramètre descriptif (et non calculé) pour la création de la sphère [0,pi]
        phi (list) : idem mais [0,2pi]

    Returns :
        tuple(float) : (x,y,z) - permet de fournir les éléments nécessaire pour la création de la sphère fixe (si pas de plot), ou de l'afficher directement. 
    """
    # Variable d'affichage initialisation [à voir si c'est utile]
    x_tab_fixe = np.zeros(N)
    y_tab_fixe = np.zeros(N)
    z_tab_fixe = np.zeros(N)

    # Variable d'affichage évolution
    x_tab_fixe = R_fixe*np.sin(theta)*np.cos(phi)
    y_tab_fixe = R_fixe*np.sin(theta)*np.sin(phi)
    z_tab_fixe = R_fixe*np.cos(theta)

    return (x_tab_fixe, y_tab_fixe, z_tab_fixe)

def sphere_physique(R,theta,phi):
    """
    Représentation d'une sphère oscillante physique (sans le plot ou avec, à voir) via les grandeurs utilisées.
    Permet de mettre en évidence l'étude futur de l'oscillation en 'temps réel'.

    Args :
        R (list) : Rayon de la sphère au cours de l'évolution. Représentation itérative.
        theta (list) : Paramètre descriptif (et non calculé) pour la création de la sphère [0,pi].
        phi (list) : idem mais [0,2pi].

    Returns :
        list : () - permet de fournir les éléments nécessaire pour la création de la sphère (si pas de plot), ou de l'afficher directement.
    """
    # Variable d'affichage initialisation [à voir si c'est utile]
    x_mat = np.zeros((N,N))
    y_mat = np.zeros((N,N))
    z_mat = np.zeros((N,N))

    # Variable d'affichage évolution 
    for n in range (N):
        x_mat[n,:] = R[n]*np.sin(theta)*np.cos(phi)
        y_mat[n,:] = R[n]*np.sin(theta)*np.sin(phi)
        z_mat[n,:] = R[n]*np.cos(theta)

        # ou 
        # x_mat[n,:] = R[n]*np.sin(theta[:])*np.cos(phi[:])
        # y_mat[n,:] = R[n]*np.sin(theta[:])*np.sin(phi[:])
        # z_mat[n,:] = R[n]*np.cos(theta[:])

    return (x_mat, y_mat, z_mat)

#=========================================
# Etude du système
#=========================================
    

if __name__ == "__main__" :
    
    # Détermination de grandeur modulable
    
    # Grandeurs
    c = 3e8                         # m.s-1         - vitesse de la lumière
    gamma = 2                       # 1             - [facteur à définir]
    K = 4.25e-3                     # m5.kg^-1.s^-2 - [facteur à défini]
    G = 6.67e-11                    # m3.kg^-1.s^-2 - constante de gravitation
    Ms = 1.98e30                    # kg            - masse du soleil 
    M_ns_moy = 1.4*Ms               # kg            - masse moyenne de l'étoile à neutron suposée constante
    
    # variables initiales
    v0 = c                          # m.s^-1        - vitesse initiale de l'oscillation
    rho0 = 1.66e17                  # kg.m^-3       - densité initiale de l'étoile à neutron 
    # U0 = v0**2                      # m2^s-1        - Potentiel initial 
    p0 = K*rho0**gamma              # kg.m^-1.s^-2  - préssion initiale
    R0 = 10e3                       # m             - rayon initial de l'étoile à neutron 
    t0 = R0/v0                      # s             - temps initial d'étude
    xi0 = 1/R0                      # à vérifier (m^-1) # initialiser correctement

    # Nombre d'itération
    N = 1000
    t = np.linspace(0,100,N) # s
    xi = np.linspace(0,1,N) # 1 
        
    # Paramètres élémentaires    
    dxi = xi[-1]/len(xi)    # Pas de [Paramètre]
    dt = t[-1]/len(t)       # Pas de temps

    # grandeur numérique remplacement
    tau = np.random.uniform() # 1
    r = xi*R0

    # Implémentation de sigma, tau, w et r dans la suite, ainsi que dans les fonctions précédentes
    # Relecture de le document de référence

    # variables à l'équilibre
    rho_c = 1.9891e18                                       # kg/m-3 - densité de l'étoile à neutron à l'équilibre                     
    r0 = np.sqrt((np.pi*K)/(2*G))                           # m - rayon numérique initiale    
    rho0_eq = rho_c*np.sinc(np.pi*(r/r0))
    U0_eq = -2*K*rho_c*(1+np.sinc(np.pi*(r/r0)))
    w_eq = 0 

    # Variables d'algorithme 
    h = tau*dt # durée donc cohérence avec l'oscillation des étoiles minutes, heures... ? 
    Y = np.array([rho0,v0,R0]) # RK4
    f = [F_global_continuite, F_global_euler, F_global_geometrique]

    # Obtention des grandeurs à étudier 

    res_grandeur = AB3_spatiotemp_base(f,Y,h,N)

    rho_evolution = [Yg[0] for Yg in res_grandeur] 
    v_evolution = [Yg[1] for Yg in res_grandeur] 
    R_evolution = [Yg[2] for Yg in res_grandeur] 

    # Paramètre pour les sphères représentatives
    
    # Variable repère sphérique
    theta = np.arange(0,np.pi,N) # rad
    phi = np.arange(0,2*np.pi,N) # rad






    
    
    
    
    
    
    
    
    
    
    