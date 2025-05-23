#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 16:02:43 2025

@author: giacomovargiu
"""
import numpy as np
from functions_plan import orbit_3D, state_vector, rotational_matrix, orbital_elements
#from Lambert_Bate import Lambert_Bate_1971_Z_BI_FG_SIMPLE
from lambert_izzo_2015_x_hh_rt_norel import lambert_izzo_2015_x_hh_rt_norel
from ANGLE import angle_from_coord


def orbital_parameters(r, v, mu):
    """
    Calcula los elementos orbitales dados el vector posición r1, el vector velocidad v1 y el parámetro gravitacional mu.
    """
    # Momento angular
    h = np.cross(r, v)
    
    # Vector de excentricidad y excentricidad
    e_vector = np.cross(v, h) / mu - r / np.linalg.norm(r)
    e = np.linalg.norm(e_vector)
    
    # Nodo (producto cruz entre k y h1, donde k es el eje z)
    k = np.array([0, 0, 1])
    N = np.cross(k, h)
    
    # Inclinación
    i = np.arccos(h[2] / np.linalg.norm(h))
    
    # Longitud del nodo ascendente
    if N[1] >= 0:
        Omega = np.arccos(N[0] / np.linalg.norm(N))
    else:
        Omega = 2 * np.pi - np.arccos(N[0] / np.linalg.norm(N))
    
    # Argumento del periapsis
    if e_vector[2] >= 0:
        omega = np.arccos(np.dot(N, e_vector) / (np.linalg.norm(N) * e))
    else:
        omega = 2 * np.pi - np.arccos(np.dot(N, e_vector) / (np.linalg.norm(N) * e))
    
    # Semieje mayor
    a = 1 / ((2 / np.linalg.norm(r)) - (np.linalg.norm(v) ** 2) / mu)
    
    # True anomaly at departure
    if np.dot(r,v) >= 0:
        theta = np.arccos(np.dot(e_vector,r)/(e * np.linalg.norm(r)));
    else:
        theta = 2*np.pi - np.arccos(np.dot(e_vector,r)/(e * np.linalg.norm(r)));
    
    return a, e, Omega, omega, i, theta, h

def data_plots():
    AU = 149597870700  # 1 UA en metros
    mu = 1.327e20      # m^3/s^2

    # Diccionarios para nombres y colores de planetas
    planet_names = {
        1: "Mercury", 2: "Venus", 3: "Earth", 4: "Mars",
        5: "Jupiter", 6: "Saturn", 7: "Uranus", 8: "Neptune"
    }
    colors = {
    1: "lightgray",             # Mercurio
    2: "lightcoral",            # Venus
    3: "lightskyblue",  # Tierra
    4: "lemonchiffon",          # Marte
    5: "peachpuff",             # Júpiter
    6: "powderblue",            # Saturno
    7: "lemonchiffon",          # Urano
    8: "plum"                   # Neptuno
}
    # Masas aproximadas de planetas (kg) y del Sol (kg)
    planet_masses = {1: 3.3011e23,  # Mercury 
                     2: 4.8675e24,  # Venus
                     3: 5.972e24,   # Earth
                     4: 6.4171e23,  # Mars
                     5: 1.8982e27,  # Jupiter
                     6: 5.6834e26,  # Saturn
                     7: 8.6810e25,  # Uranus
                     8: 1.02413e26  # Neptune
                     }
    sun_mass = 1.989e30  # kg
    return AU, mu, planet_names, colors, planet_masses, sun_mass



def Lambert_3D(JC_dep, JC_arr, departure_planet, arrival_planet, tof_sec):
    
    # Data
    AU, mu, planet_names, colors, planet_masses, sun_mass = data_plots()
    
    # Obtener las coordenadas 3D para la órbita de departure_planet
    x_dep, y_dep, z_dep = orbit_3D(JC_dep, departure_planet)
    # Obtener las coordenadas 3D para la órbita de arrival_planet
    x_arr, y_arr, z_arr = orbit_3D(JC_arr, arrival_planet)

    # Obtener la posición y estado del planeta de partida
    r_dep, v_dep, theta_dep = state_vector(departure_planet, JC_dep)
    # Obtener la posición y estado del planeta de llegada
    r_arr, v_arr, theta_arr = state_vector(arrival_planet, JC_arr)
    

    # ================================================= #
    # Cálculo y graficado de la órbita de transferencia #
    # ================================================= #
    a_dep, e_dep, i_dep, omega_dep, Omega_dep, theta_dep = orbital_elements(departure_planet, JC_dep)
    a_arr, e_arr, i_arr, omega_arr, Omega_arr, theta_arr = orbital_elements(arrival_planet, JC_arr)
    
    transfer_angle = angle_from_coord(r_dep, r_arr)
    k = 0 if transfer_angle <= np.pi else 1
    
    #v1_transfer, v2_transfer = Lambert_Bate_1971_Z_BI_FG_SIMPLE(r_dep, r_arr, tof_sec, mu, 0, 1E-6 )
    v1_transfer, v2_transfer, flag, iter_count = lambert_izzo_2015_x_hh_rt_norel(r_dep, r_arr, tof_sec, mu, k, 0, 0, 1E-6)

    # Calcular los parámetros orbitales de la órbita de transferencia usando el estado de salida modificado
    a_transfer1, e_transfer1, Omega_transfer1, omega_transfer1, i_transfer1, theta_transfer1, h_transfer1 = orbital_parameters(r_dep, v1_transfer, mu)
    a_transfer2, e_transfer2, Omega_transfer2, omega_transfer2, i_transfer2, theta_transfer2, h_transfer2 = orbital_parameters(r_arr, v2_transfer, mu)
    # Generar puntos para la órbita de transferencia
    # Hay que tener en cuenta la rotación que se aplicará después según omega, Omega e i
    if theta_transfer1 > theta_transfer2:
        theta_transfer2 += 2 * np.pi
    f_vals_transfer = np.linspace(theta_transfer1, theta_transfer2, 500)  # radianes
    # Ecuación de la órbita: r = a(1 - e²) / (1 + e*cos(f))
    r_vals_transfer = a_transfer1 * (1 - e_transfer1**2) / (1 + e_transfer1 * np.cos(f_vals_transfer))
  
    x_orb = r_vals_transfer * np.cos(f_vals_transfer)
    y_orb = r_vals_transfer * np.sin(f_vals_transfer)
    z_orb = np.zeros_like(x_orb)  # todo en el plano orbital
                
    # Construir la matriz de rotación en 3D
    R = rotational_matrix(Omega_transfer1, omega_transfer1, i_transfer1)
    
    # Apilar las coordenadas de la órbita en un array 3xN
    coords_orb = np.vstack((x_orb, y_orb, z_orb))
    
    # Aplicar la matriz de rotación para pasar a coordenadas inerciales
    coords_rot = R @ coords_orb
    
    # Separar las componentes X, Y, Z ya rotadas
    x_transfer_rot = coords_rot[0, :]
    y_transfer_rot = coords_rot[1, :]
    z_transfer_rot = coords_rot[2, :]
    
    # ============================================================= #
    # Cálculo del punto de intersección con la esfera de influencia #
    # ============================================================= #
    
    # DEPARTURE SOI
    
    # Sphere of Influence (SOI)
    a_AU1 = np.linalg.norm(r_dep) / AU
    m1 = planet_masses[departure_planet]
    r_SOI1 = a_AU1 * (m1 / sun_mass) ** (2/5)  # Radio de la SOI en UA
    
    # Convertir las coordenadas de la órbita de transferencia a UA
    x_transfer_AU1 = x_transfer_rot / AU
    y_transfer_AU1 = y_transfer_rot / AU
    z_transfer_AU1 = z_transfer_rot / AU

    # Posición del planeta de llegada en UA
    r_dep_AU = r_dep / AU

    # Calcular la distancia de cada punto de la transferencia al planeta de llegada
    distances1 = np.sqrt((x_transfer_AU1 - r_dep_AU[0])**2 +
                        (y_transfer_AU1 - r_dep_AU[1])**2 +
                        (z_transfer_AU1 - r_dep_AU[2])**2)
    
    # Buscar el índice donde la distancia es la más cercana al radio de la SOI
    idx1 = np.argmin(np.abs(distances1 - r_SOI1))
    intersection_SOI1 = (x_transfer_AU1[idx1], y_transfer_AU1[idx1], z_transfer_AU1[idx1])
    
    
    # ARRIVAL SOI
    
    # Sphere of Influence (SOI)
    a_AU2 = np.linalg.norm(r_arr) / AU
    m2 = planet_masses[arrival_planet]
    r_SOI2 = a_AU2 * (m2 / sun_mass) ** (2/5)  # Radio de la SOI en UA
    
    # Convertir las coordenadas de la órbita de transferencia a UA
    x_transfer_AU2 = x_transfer_rot / AU
    y_transfer_AU2 = y_transfer_rot / AU
    z_transfer_AU2 = z_transfer_rot / AU

    # Posición del planeta de llegada en UA
    r_arr_AU = r_arr / AU

    # Calcular la distancia de cada punto de la transferencia al planeta de llegada
    distances2 = np.sqrt((x_transfer_AU2 - r_arr_AU[0])**2 +
                        (y_transfer_AU2 - r_arr_AU[1])**2 +
                        (z_transfer_AU2 - r_arr_AU[2])**2)
    
    # Buscar el índice donde la distancia es la más cercana al radio de la SOI
    idx2 = np.argmin(np.abs(distances2 - r_SOI2))
    intersection_SOI2 = (x_transfer_AU2[idx2], y_transfer_AU2[idx2], z_transfer_AU2[idx2])
    
    
    #print((f"DELTA_THETA: {transfer_angle:.2f} rad"))
    
    return h_transfer1, intersection_SOI1, intersection_SOI2, v1_transfer, v2_transfer, r_arr, r_dep, x_dep, y_dep, z_dep, x_arr, y_arr, z_arr, x_transfer_rot, y_transfer_rot, z_transfer_rot

