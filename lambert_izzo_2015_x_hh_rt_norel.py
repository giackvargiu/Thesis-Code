#!/usr/bin/env python3
"""
lambert_izzo_2015_x_hh_rt_norel.py

Translation of the MATLAB implementation of the Lambert problem solver
(Lambert_Izzo_2015_X_HH_RT_NOREL) into Python3.

Algorithm:
    D. Izzo, Revisiting Lambert's problem, Celest. Mech. Dyn. Astron., vol. 121,
    no. 1, pp. 1–15, Jan. 2015.

The code uses the Householder iterative method and includes subroutines
for computing time-of-flight from the universal variable x and its derivatives.

Functions:
    lambert_izzo_2015_x_hh_rt_norel(r1, r2, tof, mu, lw, mr, lp, tol)
        Main function to solve the Lambert problem.

    lambert_izzo_2015_x_tmin(r1, r2, mu, lw, mr)
        Computes the minimum time-of-flight for a multirevolution transfer.

    dtdx(x, t, lam2, lam3, lam5)
        Computes the 1st–3rd derivatives for the Householder iterations.

    x2tof(x, n, lam, lam2)
        Computes the time-of-flight corresponding to the universal variable x.

    hypergeo2f1(x)
        Computes the hypergeometric series 2F1(3,1,5/2,x).
"""

import numpy as np


def lambert_izzo_2015_x_hh_rt_norel(r1, r2, tof, mu, lw, mr, lp, tol):
    """
    Solves the Lambert problem.

    Parameters:
        r1  : Position vector at departure [list or array of 3 floats].
        r2  : Position vector at arrival [list or array of 3 floats].
        tof : Time of flight.
        mu  : Gravitational parameter.
        lw  : Transfer type (0: type I (short-way), 1: type II (long-way)).
        mr  : Number of revolutions (0, 1, 2, ...).
        lp  : Multi-revolution orbit type (0: short-period, 1: long-period).
        tol : Tolerance for iteration (absolute tolerance on time error).

    Returns:
        v1    : Final velocity vector at departure.
        v2    : Initial velocity vector at arrival.
        flag  : Status flag (0=OK, 201,202,203,... see inline comments).
        i     : Number of iterations performed.
    """
    # Convert to numpy arrays
    r1 = np.array(r1, dtype=float)
    r2 = np.array(r2, dtype=float)

    # Initial internal parameters
    flag = 0  # Status flag
    iter_count = 0  # Iteration counter
    ni = 100  # Maximum number of iterations for Householder method
    v1 = np.array([np.nan, np.nan, np.nan])  # Preallocate output vectors
    v2 = np.array([np.nan, np.nan, np.nan])

    # Auxiliary magnitudes
    r1n = np.linalg.norm(r1)
    r2n = np.linalg.norm(r2)
    lws = -(2 * lw - 1)  # Long-way sign (lw=0 -> 1, lw=1 -> -1)

    # Check validity of inputs
    if r1n == 0.0 or r2n == 0.0:
        flag = 201  # Bad input: position vector with zero norm
        return v1, v2, flag, iter_count
    if tof <= 0.0:
        flag = 202  # Negative time of flight
        return v1, v2, flag, iter_count
    if mu <= 0.0:
        flag = 203  # Negative gravitational parameter
        return v1, v2, flag, iter_count
    if tol <= 0.0:
        flag = 206  # Negative tolerance
        return v1, v2, flag, iter_count

    # Compute geometrical parameters
    h = np.cross(r1, r2)  # Orbit normal vector
    c = np.linalg.norm(r2 - r1)  # Length of chord P1_P2
    s = 0.5 * (r1n + r2n + c)  # Semiperimeter of triangle P1_P2_F
    lam2 = 1 - c / s  # Battin's Lambda parameter, squared
    lam = lws * np.sqrt(lam2)  # Lambda parameter
    lam3 = lam2 * lam  # Powers of Lambda
    lam5 = lam3 * lam2  # Powers of Lambda
    # Non-dimensional time-of-flight
    t_non_dim = np.sqrt(2 * mu / (s ** 3)) * tof

    # Account for pi transfers
    phi = np.arccos(np.dot(r1, r2) / (r1n * r2n))
    if np.isclose(phi % np.pi, 0, atol=1e-12):
        h = np.array([0.0, 0.0, 1.0])  # Arbitrary plane on pi transfers

    # Initial guess for the universal variable x
    mmax = int(np.floor(t_non_dim / np.pi))
    t00 = np.arccos(lam) + lam * np.sqrt(1 - lam2)  # Eq. 19 in [1]
    t0m = t00 + mmax * np.pi  # Eq. 19 in [1]
    t1 = 2.0 / 3.0 * (1 - lam3)  # TOF for a parabolic orbit

    # For multi-revolution transfers: check if transfer is possible and compute Tmin
    if mr > 0:
        # Get Tmin via Halley root solver from a separate function
        x_tmin, tmin, flag_tmin, iter1 = lambert_izzo_2015_x_tmin(r1, r2, mu, lw, mr)
        # Normalize Tmin magnitude
        tmin = tmin * np.sqrt(2 * mu / (s ** 3))
        if tmin > t_non_dim:
            flag = 700  # No solution exists (time-of-flight too low)
            return v1, v2, flag, iter_count

    # Compute starter (initial guess)
    if mr == 0:  # Single-revolution starter
        if t_non_dim >= t0m:
            x0 = -(t_non_dim - t00) / (t_non_dim - t00 + 4)
        elif t_non_dim <= t1:
            x0 = 2.5 * (t1 * (t1 - t_non_dim)) / (t_non_dim * (1 - lam5)) + 1
        else:
            exponent = np.log(2) / np.log(t1 / t00)
            x0 = (t_non_dim / t00) ** exponent - 1
    else:  # Multi-revolution starter, from inverting Eq. 29 in [1]
        if not lp:  # Left branch (short-period)
            term = ((mr * np.pi + np.pi) / (8 * t_non_dim)) ** (2 / 3)
            x0 = (term - 1) / (term + 1)
        else:  # Right branch (long-period)
            term = ((8 * t_non_dim) / (mr * np.pi)) ** (2 / 3)
            x0 = (term - 1) / (term + 1)
    x = x0

    # Householder iterative solver for x
    for i in range(ni):
        # Safety check: x out-of-bounds
        if x <= -1 or (mr > 1 and x >= 1):
            flag = 303
            return v1, v2, flag, iter_count

        # Compute time-of-flight for current x
        ti = x2tof(x, mr, lam, lam2)
        # Compute derivatives
        dt, d2t, d3t = dtdx(x, ti, lam2, lam3, lam5)

        d = ti - t_non_dim
        d2 = d * d
        dt2 = dt * dt
        denom = dt * (dt2 - d * d2t) + d3t * d2 / 6.0
        term = d * (dt2 - d * d2t / 2.0) / (denom if denom != 0 else 1)
        x = x - term
        iter_count = i + 1
        # Check convergence conditions
        if abs(d) <= tol:
            break
        if abs(term) <= np.finfo(float).eps:
            flag = 308
            break
        if i == ni - 1:  # Max iterations exceeded
            flag = 302
            return v1, v2, flag, iter_count

    # Compute velocity vectors from the converged solution
    r1u = r1 / r1n  # Radial unit vector at r1
    r2u = r2 / r2n  # Radial unit vector at r2
    hxr1 = lws * np.cross(h, r1)
    norm_hxr1 = np.linalg.norm(hxr1)
    hxr1u = hxr1 / norm_hxr1 if norm_hxr1 != 0 else hxr1
    hxr2 = lws * np.cross(h, r2)
    norm_hxr2 = np.linalg.norm(hxr2)
    hxr2u = hxr2 / norm_hxr2 if norm_hxr2 != 0 else hxr2

    gamma = np.sqrt(mu * s / 2.0)
    rho = (r1n - r2n) / c
    sigma = np.sqrt(1 - rho ** 2)
    omx2 = 1 - x ** 2
    y_val = np.sqrt(1 - lam2 * omx2)
    term1 = lam * y_val - x
    term2 = rho * (lam * y_val + x)
    vr1n = +gamma * (term1 - term2) / r1n
    vr2n = -gamma * (term1 + term2) / r2n
    vt = gamma * sigma * (y_val + lam * x)
    vt1n = vt / r1n
    vt2n = vt / r2n

    v1 = vr1n * r1u + vt1n * hxr1u
    v2 = vr2n * r2u + vt2n * hxr2u

    return v1, v2, flag, iter_count


def dtdx(x, t_value, lam2, lam3, lam5):
    """
    Computes the 1st to 3rd derivatives for the Householder iterations.

    Parameters:
        x      : Universal variable.
        t_value: Time-of-flight for current x.
        lam2   : lambda^2 (Battin's parameter squared).
        lam3   : lambda^3.
        lam5   : lambda^5.

    Returns:
        dt   : First derivative.
        d2t  : Second derivative.
        d3t  : Third derivative.
    """
    omx2 = 1 - x * x
    omx2inv = 1.0 / omx2 if omx2 != 0 else np.inf
    y = np.sqrt(1 - lam2 * omx2)  # [1] p. 6
    if x == 0:
        dt = -2.0  # Minimum energy ellipse case
    elif x == 1:
        dt = 0.4 * (lam5 - 1)  # Parabolic orbit case (x=1)
    else:
        dt = omx2inv * (3 * t_value * x - 2 + 2 * lam3 * x / y)
    d2t = omx2inv * (3 * t_value + 5 * x * dt + 2 * (1 - lam2) * lam3 / (y ** 3))
    d3t = omx2inv * (7 * x * d2t + 8 * dt - 6 * (1 - lam2) * lam5 * x / (y ** 5))
    return dt, d2t, d3t


def x2tof(x, n, lam, lam2):
    """
    Computes the time-of-flight from the universal variable x.

    Parameters:
        x    : Universal variable.
        n    : Number of revolutions.
        lam  : Lambda parameter.
        lam2 : lambda squared.

    Returns:
        t_val: Time-of-flight corresponding to x.
    """
    battin = 0.01
    lagrange = 0.2
    dist = abs(x - 1)

    # Use Lagrange TOF expression if within the specified interval
    if lagrange > dist > battin:
        a = 1 / (1 - x * x)  # Semi-major axis
        if a > 0:  # Elliptic orbit
            alpha = 2 * np.arccos(x)  # Eq. 10 in [1]
            beta = 2 * np.arcsin(np.sqrt(lam2 / a))  # Eq. 11 in [1]
            if lam < 0:
                beta = -beta
            t_val = (a * np.sqrt(a) * ((alpha - np.sin(alpha)) - (beta - np.sin(beta)) + 2 * np.pi * n)) / 2
        else:  # Hyperbolic orbit
            alpha = 2 * np.arccosh(x)
            beta = 2 * np.arcsinh(np.sqrt(-lam2 / a))
            if lam < 0:
                beta = -beta
            t_val = -a * np.sqrt(-a) * ((beta - np.sinh(beta)) - (alpha - np.sinh(alpha))) / 2
    else:
        # Use Battin's / Lancaster TOF expressions
        e = x * x - 1
        rho_val = abs(e)
        z = np.sqrt(1 + lam2 * e)
        if dist < battin:
            eta = z - lam * x
            s1 = 0.5 * (1 - lam - x * eta)
            q = (4 / 3) * hypergeo2f1(s1)
            t_val = (eta ** 3 * q + 4 * lam * eta) / 2 + n * np.pi / (rho_val * np.sqrt(rho_val))
        else:
            y_val = np.sqrt(rho_val)
            g = x * z - lam * e
            if e < 0:
                l_val = np.arccos(g)
                d_val = n * np.pi + l_val
            else:
                f_val = y_val * (z - lam * x)
                d_val = np.log(f_val + g)
            t_val = (x - lam * z - d_val / y_val) / e
    return t_val


def hypergeo2f1(x):
    """
    Computes the hypergeometric series 2F1(3,1,5/2,x).

    Parameters:
        x : Argument for the hypergeometric series.

    Returns:
        f_val : Value of the hypergeometric series.
    """
    tol = 1e-9
    kmax = 1000
    f_val = 1.0
    term = 1.0
    for k in range(kmax + 1):
        term = term * (3 + k) * (1 + k) / (2.5 + k) * x / (k + 1)
        f_val = f_val + term
        if abs(term) <= tol:
            break
    return f_val


def lambert_izzo_2015_x_tmin(r1, r2, mu, lw, mr):
    """
    Computes the minimum time-of-flight (Tmin) for a multi-revolution transfer.

    Parameters:
        r1 : Position vector at departure.
        r2 : Position vector at arrival.
        mu : Gravitational parameter.
        lw : Transfer type (0: type I, 1: type II).
        mr : Number of revolutions.

    Returns:
        x       : Universal variable value for Tmin.
        tofb    : Minimum time-of-flight value.
        flag    : Status flag (0 if OK, nonzero otherwise).
        iter_count : Number of iterations used in the Halley solver.
    """
    r1 = np.array(r1, dtype=float)
    r2 = np.array(r2, dtype=float)

    flag = 0
    tol_halley = 1e-13
    itermax_halley = 100
    tofb = np.nan
    r1n = np.linalg.norm(r1)
    r2n = np.linalg.norm(r2)
    lws = -(2 * lw - 1)

    c = np.linalg.norm(r2 - r1)
    s = 0.5 * (r1n + r2n + c)
    lam2 = 1.0 - c / s
    lam = lws * np.sqrt(lam2)
    lam3 = lam2 * lam
    lam5 = lam3 * lam2

    t00 = np.arccos(lam) + lam * np.sqrt(1 - lam2)  # Eq. 19 in [1]
    t0 = t00 + mr * np.pi

    if mr > 0:
        x_old = 0.0
        t_value = t0
        x_new = np.nan
        tmin = np.nan
        iter_count = 0
        for i in range(itermax_halley):
            dt, d2t, d3t = dtdx(x_old, t_value, lam2, lam3, lam5)
            if dt != 0:
                x_new = x_old - dt * d2t / (d2t * d2t - dt * d3t / 2.0)
            else:
                x_new = x_old
            if abs(dt) < tol_halley:
                break
            if i >= itermax_halley - 1:
                flag = 3  # Maximum iterations exceeded
                return x_new, tofb, flag, i + 1
            tmin = x2tof(x_new, mr, lam, lam2)
            x_old = x_new
            iter_count = i + 1
        x_val = x_new
        tofb = tmin / np.sqrt(2.0 * mu / (s ** 3))
        return x_val, tofb, flag, iter_count
    else:
        # For mr = 0, Tmin is not computed here.
        return 0, 0, flag, 0
