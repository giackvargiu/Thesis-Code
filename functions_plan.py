import math
import numpy as np
import matplotlib.pyplot as plt

def planet_OP(number):
    # [a, a_dot, e, e_dot, i, i_dot, Omega, Omega_dot, omega, omega_dot, L, L_dot]
    planets = {
        # Mercury
        1: [0.38709927, 0.00000037, 0.20563593, 0.00001906, 7.00497902, -0.00594749,
            48.33076593, -0.12534081, 77.45779628, 0.16047689, 252.25032350, 149472.67411175],
        # Venus
        2: [0.72333566, 0.00000390, 0.00677672, -0.00004107, 3.39467605, -0.00078890,
            76.67984255, -0.27769418, 131.60246718, 0.00268329, 181.97909950, 58517.81538729],
        # Earth
        3: [1.00000261, 0.00000562, 0.01671123, -0.00004391, -0.00001531, -0.01294668,
            0.0, 0.0, 102.93768193, 0.32327364, 100.46457166, 35999.37244981],
        # Mars
        4: [1.52371034, 0.0001847, 0.09339410, 0.00007882, 1.84969142, -0.00813131,
            49.55953891, -0.29257343, -23.94362959, 0.44441088, -4.55343205, 19140.30268499],
        # Jupiter
        5: [5.20288700, -0.00011607, 0.04838624, -0.00013253, 1.30439695, -0.00183714,
            100.47390909, 0.20469106, 14.72847983, 0.21252668, 34.39644501, 3034.74612775],
        # Saturn
        6: [9.53667594, -0.00125060, 0.05386179, -0.00050991, 2.48599187, 0.00193609,
            113.66242448, -0.28867794, 92.59887831, -0.41897216, 49.95424423, 1222.49362201],
        # Uranus
        7: [19.18916464, -0.00196176, 0.04725744, -0.00004397, 0.77263783, -0.00242939,
            74.01692503, 0.04240589, 170.95427630, 0.40805281, 313.23810451, 428.48202785],
        # Neptune
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

def julian_day(JC):
    JD = JC * 36525 + 2451545
    return JD

def propagation(E, E_dot, jul_cent):
    return E + E_dot * jul_cent

def angular_momentum(mu, a, e):
    return math.sqrt(mu * a * (1 - e**2))

def mean_anomaly(omega_prop, L_prop):
    return L_prop - omega_prop

def true_anomaly(e, E):
    return 2.0 * math.atan(math.tan(E / 2.0) * math.sqrt((1 + e) / (1 - e)))

def eccentric_anomaly(M, e, N, tolerance, E0):
    """
    Resuelve la ecuación de Kepler mediante el método de Newton-Raphson.
    
    Parámetros:
        M         : Mean anomaly (valor escalar)
        e         : Excentricidad (valor escalar)
        N         : Número máximo de iteraciones
        tolerance : Tolerancia para la convergencia
        E0        : Valor inicial (guess)
        
    Devuelve:
        E_solution : Valor final de E que resuelve la ecuación
        iterations : Número de iteraciones realizadas
        E_values   : Vector con los valores de E en cada iteración
    """
    E = E0
    
    for i in range(1, N+1):
        # Evaluamos la ecuación de Kepler y su derivada
        f = E - e * np.sin(E) - M
        dot_f = 1 - e * np.cos(E)
        
        # Actualización de Newton-Raphson
        E_next = E - f / dot_f
        
        if np.abs(E_next - E) < tolerance:
            return E_next
        else:
            E = E_next
    
    # Si no se converge, se emite una advertencia y se retorna el último valor
    print("Advertencia: El método de Newton-Raphson no ha convergido.")
    return E

def rotational_matrix(Omega, omega, inclination):
    R = np.zeros((3, 3))
    
    R[0, 0] = np.cos(Omega) * np.cos(omega) - np.sin(Omega) * np.cos(inclination) * np.sin(omega)
    R[0, 1] = -np.cos(Omega) * np.sin(omega) - np.sin(Omega) * np.cos(inclination) * np.cos(omega)
    R[0, 2] = np.sin(Omega) * np.sin(inclination)

    R[1, 0] = np.sin(Omega) * np.cos(omega) + np.cos(Omega) * np.cos(inclination) * np.sin(omega)
    R[1, 1] = -np.sin(Omega) * np.sin(omega) + np.cos(Omega) * np.cos(inclination) * np.cos(omega)
    R[1, 2] = -np.cos(Omega) * np.sin(inclination)

    R[2, 0] = np.sin(inclination) * np.sin(omega)
    R[2, 1] = np.sin(inclination) * np.cos(omega)
    R[2, 2] = np.cos(inclination)

    return R


def arg_per(omega_prop, Omega_prop):
    return omega_prop - Omega_prop

def kepler(M, E, e):
    return E - e * np.sin(E) - M

def kepler_der(E, e):
    return 1.0 - e * np.cos(E)
    
def orbital_elements(number, JC):
    op = planet_OP(number)
    a, e, i, Omega, varpi = propagation_OP(number, JC)

    # Longitud media propagada
    L = propagation(op[10], op[11], JC)        # en grados
    L = math.radians(L) % (2*math.pi)

    # Argumento de perihelio = ϖ − Ω
    omega = (varpi - Omega) % (2*math.pi)

    # Anomalía media = L − ϖ
    M = (L - varpi) % (2*math.pi)

    # Resolver Kepler para E y luego obtener θ
    E = eccentric_anomaly(M, e, 100, 1e-8, M)
    theta = true_anomaly(e, E)

    return a, e, i, omega, Omega, theta

def position(julian_century):
    planets = {}
    # Cálculo de las posiciones de cada planeta usando la función `state_vector`
    for i in range(1, 9):
        r, v, theta = state_vector(i, julian_century)
        planets[i] = r
    return planets

def propagation_OP(i, JC):
    # Dado un planeta y una epoc, devuelve los parámetros orbitales de su órbita:
    # a, e, i, omega, Omega, true anomaly
    
    AU = 149597870700  # 1 AU en metros
    op = planet_OP(i)
    a_prop = propagation(op[0], op[1], JC) * AU  # semieje mayor en metros
    e_prop = propagation(op[2], op[3], JC)
    # Convertir a radianes la inclinación, ascensión del nodo y argumento del perihelio
    i_prop = math.radians(propagation(op[4], op[5], JC))
    omega_prop = math.radians(propagation(op[8], op[9], JC)) % (2 * math.pi) # normalizarlo en el rango [0, 2*pi]
    Omega_prop = math.radians(propagation(op[6], op[7], JC)) % (2 * math.pi) # normalizarlo en el rango [0, 2*pi]
    # falta true anomaly
    return a_prop, e_prop, i_prop, omega_prop, Omega_prop

def orbit_xyz(a_prop, e_prop, f_vals):
    # f_vals puede ser un vector (para una órbita entera) o una valor único (para un punto exacto en la órbita)
    
    # Ecuación de la órbita: r = a*(1 - e^2) / (1 + e*cos(f))
    r_vals = a_prop * (1 - e_prop**2) / (1 + e_prop * np.cos(f_vals))

    # Coordenadas en el plano orbital (el perihelio se alinea con el eje x)
    x_orb = r_vals * np.cos(f_vals)
    y_orb = r_vals * np.sin(f_vals)
    z_orb = np.zeros_like(x_orb)  # plano orbital: z=0
    
    return x_orb, y_orb, z_orb

def state_vector(i, jul_cent):
    # Constantes
    mu = 1.32712440018e20

    a_prop, e_prop, i_prop, omega_prop, Omega_prop, theta = orbital_elements(i, jul_cent)

    # Cálculo del momento angular
    h = angular_momentum(mu, a_prop, e_prop)

    # Calcular el vector posición
    x_orb, y_orb, z_orb = orbit_xyz(a_prop, e_prop, theta)
    
    r = np.zeros(3)
    r[0] = x_orb
    r[1] = y_orb
    r[2] = z_orb
    
    
    # Calcular el vector velocidad
    v = np.zeros(3)
    v[0] = -(mu / h) * np.sin(theta)
    v[1] = (mu / h) * (e_prop + np.cos(theta))
    v[2] = 0.0

    # Aplicar la matriz de rotación para obtener el estado en el sistema de referencia inercial
    R = rotational_matrix(Omega_prop, omega_prop, i_prop)
    r = np.dot(R, r)
    v = np.dot(R, v)

    return r, v, theta


def orbit_3D(JC, i):
    AU = 149597870700  # 1 AU en metros

    # Obtener parámetros orbitales propagados para el planeta
    a_prop, e_prop, i_prop, omega_prop, Omega_prop, theta = orbital_elements(i, JC)
    
    # Generar ángulos verdaderos de 0 a 2pi para definir la órbita en el plano orbital
    f_vals = np.linspace(0, 2 * math.pi, 500)
    
    x_orb, y_orb, z_orb = orbit_xyz(a_prop, e_prop, f_vals)
    
    # Matriz de rotación para transformar al sistema inercial
    R = rotational_matrix(Omega_prop, omega_prop, i_prop)
    
    # Aplicar la matriz de rotación a cada punto de la órbita
    coords_inercial = np.array([x_orb, y_orb, z_orb])
    coords_rot = R @ coords_inercial  # producto matricial
    
    # Extraer las coordenadas y convertir de metros a AU
    x_orbita_AU = coords_rot[0, :] / AU
    y_orbita_AU = coords_rot[1, :] / AU
    z_orbita_AU = coords_rot[2, :] / AU
    
    return x_orbita_AU, y_orbita_AU, z_orbita_AU
    

def plot_orbits(planets, day, month, year, utch=0, utcm=0, utcs=0):
    # planets es un array con las posiciones r de cada planeta

    # Diccionarios de nombres y colores para los planetas
    planet_names = {
        1: "Mercury",
        2: "Venus",
        3: "Earth",
        4: "Mars",
        5: "Jupiter",
        6: "Saturn",
        7: "Uranus",
        8: "Neptune"
    }
    colors = {
        1: "gray",
        2: "yellow",
        3: "blue",
        4: "red",
        5: "orange",
        6: "lightblue",
        7: "gold",
        8: "purple"
    }
    
    Au = 149597870700  # 1 AU en metros
    #mu = 1.32712440018e20

    JD, JC = to_julian(year, month, day, utch, utcm, utcs)

    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Representar el Sol en el origen
    ax.plot(0, 0, marker='*', color='yellow', markersize=12,
            markerfacecolor='none', linestyle='None', label='Sun')

    # Para cada planeta, dibujar su órbita elíptica y su posición actual
    for i in sorted(planets.keys()):
       r_current = planets[i]  # Posición actual en metros
       color = colors[i]
    
       # Obtener los parámetros orbitales del planeta y propagarlos con JC
       op = planet_OP(i)
       a_prop = propagation(op[0], op[1], JC) * Au  # semieje mayor en metros
       e_prop = propagation(op[2], op[3], JC)
       # Propagar el argumento del perihelio y convertir a radianes
       omega_prop = np.mod(math.radians(propagation(op[8], op[9], JC)), 2 * math.pi)
    
    
       # Generar puntos de la órbita en el plano orbital (perihelio en el eje x)
       f_vals = np.linspace(0, 2 * math.pi, 360)
       # Ecuación de la órbita: r = a(1 - e²) / (1 + e*cos(f))
       r_vals = a_prop * (1 - e_prop**2) / (1 + e_prop * np.cos(f_vals))
       x_orb = r_vals * np.cos(f_vals)
       y_orb = r_vals * np.sin(f_vals)
    
       # Rotar la órbita según el argumento del perihelio
       x_orbita = x_orb * np.cos(omega_prop) - y_orb * np.sin(omega_prop)
       y_orbita = x_orb * np.sin(omega_prop) + y_orb * np.cos(omega_prop)
    
       # Convertir la órbita de metros a AU
       x_orbita_AU = x_orbita / Au
       y_orbita_AU = y_orbita / Au
    
       # Dibujar la órbita elíptica
       ax.plot(x_orbita_AU, y_orbita_AU, '--', color=color, alpha=0.6)
    
       # Dibujar la posición actual del planeta
       r_current_AU = r_current / Au
       ax.plot(r_current_AU[0], r_current_AU[1], 'o', color=color, label=planet_names[i])
    
       # Calcular y dibujar la SOI
       # Se utiliza la norma de la posición actual para aproximar el semieje en AU
       a_AU = np.linalg.norm(r_current) / Au
       # Escalar la SOI para una mejor visualización

    ax.set_aspect('equal', 'box')
    ax.set_xlabel('X [AU]')
    ax.set_ylabel('Y [AU]')

    # Determinar el título y los límites de los ejes según los planetas mostrados
    keys = set(planets.keys())
    if keys == {1, 2, 3, 4}:
       title_prefix = "Inner Planets of the Solar System"
       ax.set_xlim(-5, 5)
       ax.set_ylim(-5, 5)
    elif keys == {5, 6, 7, 8}:
       title_prefix = "Outer Planets of the Solar System"
       ax.set_xlim(-35, 35)
       ax.set_ylim(-35, 35)
    else:
       title_prefix = "Solar System"
       ax.set_xlim(-35, 35)
       ax.set_ylim(-35, 35)

    ax.set_title(f'{title_prefix} on the {day}/{month}/{year}')
    ax.legend(loc='upper left', fontsize='small')
    plt.show()
       
       
       