# Main code for plotting the planets of the Solar System
# Irati Cadenato
# Interplanetary trajectory Design and Optimization 2024/2025

import matplotlib.pyplot as plt
from functions_plan import to_julian, position, plot_orbits

# Date
Y = 2002
M = 6
D = 28
utch = 0.0
utcm = 0.0
utc = 0.0
julian_day, julian_century = to_julian(Y, M, D, utch, utcm, utc)

plt.style.use('dark_background')

# Positions of the planets
planets = position(julian_century)

# 1. All the planets 
plot_orbits(planets, D, M, Y)

# 2. Inner planets: Mercury (1), Venus (2), Earth (3), Mars (4)
inner_planets = {k: planets[k] for k in (1, 2, 3, 4)}
plot_orbits(inner_planets, D, M, Y)

# 3. Outer planets: Jupiter (5), Saturn (6), Uranus (7), Neptune (8)
outer_planets = {k: planets[k] for k in (5, 6, 7, 8)}
plot_orbits(outer_planets, D, M, Y)

