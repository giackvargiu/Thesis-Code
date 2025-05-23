#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 15:42:53 2025

@author: giacomovargiu
"""

import math
import numpy as np
import matplotlib.pyplot as plt

def planet_OP(number):
    # [a, a_dot, e, e_dot, i, i_dot, Omega, Omega_dot, omega, omega_dot, L, L_dot]
    planets = {
        1: [0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749,
            48.33076593, -0.12534081, 77.45779628, 0.16047689, 252.25032350, 149472.67411175],
        2: [0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890,
            76.67984255, -0.27769418, 131.60246718, 0.00268329, 181.97909950, 58517.81538729],
        3: [1.00000261, 0.00000562, 0.01671123, -0.00004391, -0.00001531, -0.01294668,
            0.0, 0.0, 102.93768193, 0.32327364, 100.46457166, 35999.37244981],
        4: [1.52371034, 0.0001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131,
            49.55953891, -0.29257343, -23.94362959, 0.44441088, -4.55343205, 19140.30268499],
        5: [5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714,
            100.47390909, 0.20469106, 14.72847983, 0.21252668, 34.39644501, 3034.74612775],
        6: [9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609,
            113.66242448, -0.28867794, 92.59887831, -0.41897216, 49.95424423, 1222.49362201],
        7: [19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939,
            74.01692503, 0.04240589, 170.95427630, 0.40805281, 313.23810451, 428.48202785],
        8: [30.06992276, 0.00026291, 0.00859048, 0.00005105, 1.77004347, 0.00035372,
            131.78422574, -0.00508664, 44.96476227, -0.32241464, -55.12002969, 218.45945325]
    }
    return planets.get(number, None)

def to_julian(y, m, d, utch, utcm, utcs):
    if m <= 2:
        y -= 1
        m += 12
    UT = utch + (utcm / 60.0) + (utcs / 3600.0)
    JD = int(365.25 * y) + int(30.6001 * (m + 1)) + d + UT / 24 + 1720981.5
    JC = (JD - 2451545) / 36525
    return JD, JC

def propagation(E, E_dot, jul_cent):
    return E + E_dot * jul_cent

def angular_momentum(mu, a, e):
    return math.sqrt(mu * a * (1 - e**2))

def mean_anomaly(omega_prop, L_prop):
    return L_prop - omega_prop

def true_anomaly(e, E):
    return 2.0 * math.atan(math.tan(E / 2.0) * math.sqrt((1 + e) / (1 - e)))

def rotational_matrix(Omega, omega, inclination):
    R = np.zeros((3, 3))
    R[0, 0] = np.cos(Omega)*np.cos(omega) - np.sin(Omega)*np.cos(inclination)*np.sin(omega)
    R[0, 1] = -np.cos(Omega)*np.sin(omega) - np.sin(Omega)*np.cos(inclination)*np.cos(omega)
    R[0, 2] = np.sin(Omega)*np.sin(inclination)
    R[1, 0] = np.sin(Omega)*np.cos(omega) + np.cos(Omega)*np.cos(inclination)*np.sin(omega)
    R[1, 1] = -np.sin(Omega)*np.sin(omega) + np.cos(Omega)*np.cos(inclination)*np.cos(omega)
    R[1, 2] = -np.cos(Omega)*np.sin(inclination)
    R[2, 0] = np.sin(inclination)*np.sin(omega)
    R[2, 1] = np.sin(inclination)*np.cos(omega)
    R[2, 2] = np.cos(inclination)
    return R

def arg_per(omega_prop, Omega_prop):
    return omega_prop - Omega_prop

def kepler(M, E, e):
    return E - e * np.sin(E) - M

def kepler_der(E, e):
    return 1.0 - e * np.cos(E)

def newton_raphson(M, e, N, eps):
    x0 = M
    for _ in range(N):
        x1 = x0 - kepler(M, x0, e) / kepler_der(x0, e)
        if abs(x1 - x0) < eps:
            return x1
        x0 = x1
    return x0

def state_vector(planet_OP, jul_cent):
    # Constantes
    mu = 1.32712440018e20
    Au = 149597870700

    # Parámetros orbitales
    a      = planet_OP[0]
    a_dot  = planet_OP[1]
    e      = planet_OP[2]
    e_dot  = planet_OP[3]
    i      = planet_OP[4]
    i_dot  = planet_OP[5]
    Omega  = planet_OP[6]
    Omega_dot = planet_OP[7]
    omega  = planet_OP[8]
    omega_dot = planet_OP[9]
    L      = planet_OP[10]
    L_dot  = planet_OP[11]

    # Propagar los parámetros orbitales
    a_prop = propagation(a, a_dot, jul_cent) * Au
    e_prop = propagation(e, e_dot, jul_cent)
    i_prop = math.radians(propagation(i, i_dot, jul_cent))
    Omega_prop = np.mod(math.radians(propagation(Omega, Omega_dot, jul_cent)), 2*math.pi)
    omega_prop = np.mod(math.radians(propagation(omega, omega_dot, jul_cent)), 2*math.pi)
    L_prop = np.mod(math.radians(propagation(L, L_dot, jul_cent)), 2*math.pi)

    # Cálculo del momento angular, anomalía media y argumento del perihelio
    h = angular_momentum(mu, a_prop, e_prop)
    omega_2 = np.mod(arg_per(omega_prop, Omega_prop), 2*math.pi)
    M = mean_anomaly(omega_prop, L_prop)

    # Resolver la ecuación de Kepler usando el método de Newton-Raphson
    N_iter = 10000000
    eps = 1e-12
    E = newton_raphson(M, e_prop, N_iter, eps)

    # Anomalía verdadera
    theta = true_anomaly(e_prop, E)

    # Calcular el vector posición
    r = np.zeros(3)
    r0 = (h**2 / mu) / (1.0 + e_prop * np.cos(theta))
    r[0] = r0 * np.cos(theta)
    r[1] = r0 * np.sin(theta)
    r[2] = 0.0

    # Calcular el vector velocidad
    v = np.zeros(3)
    v[0] = -(mu / h) * np.sin(theta)
    v[1] = (mu / h) * (e_prop + np.cos(theta))
    v[2] = 0.0

    # Aplicar la matriz de rotación para obtener el estado en el sistema de referencia inercial
    R = rotational_matrix(Omega_prop, omega_2, i_prop)
    r = np.dot(R, r)
    v = np.dot(R, v)

    return r, v, theta
    