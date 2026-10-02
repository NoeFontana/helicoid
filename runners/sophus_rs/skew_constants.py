# /// script
# dependencies = ["mpmath"]
# ///
"""The `SKEW_*` constants of `src/convert.rs`, at 60 digits, from definitions.

`uv run skew_constants.py` prints them. They are the matrix exponential of the 4x4 hat, the series
`sum(M^n / (n+1)!)` of the 3x3 hat and of `ad` (the commutator of 4x4 hats), and `mp.inverse`: no
closed form of the functions the runner asks of sophus-rs. The angle is the binary64 value the
tests use, `pi / 2 * k / 7` in that order of operations.
"""

import math

from mpmath import eye, expm, inverse, matrix, mp

mp.dps = 60

PHI = [mp.mpf(math.pi / 2 * k / 7) for k in (2, 3, 6)]
RHO = [1, 2, 3]


def hat3(v):
    return matrix([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])


def hat4(tau):
    """`[[[phi]x, rho], [0, 0]]` of the rotation-first tangent `[phi; rho]`."""
    m = matrix(4, 4)
    m[0:3, 0:3] = hat3(tau[:3])
    for i in range(3):
        m[i, 3] = tau[3 + i]
    return m


def vee4(m):
    return [m[2, 1], m[0, 2], m[1, 0], m[0, 3], m[1, 3], m[2, 3]]


def ad(tau):
    """Column `j` is `vee([hat4(tau), hat4(e_j)])`."""
    h = hat4(tau)
    columns = []
    for j in range(6):
        e = hat4([1 if i == j else 0 for i in range(6)])
        columns.append(vee4(h * e - e * h))
    return matrix(columns).T


def series(m):
    """`sum(m^n / (n+1)!)`; 80 terms, far below 1e-60 for these norms."""
    total, power, factorial = matrix(m.rows, m.rows), eye(m.rows), mp.mpf(1)
    for n in range(80):
        total += power / factorial
        power *= m
        factorial *= n + 2
    return total


def show(name, m, rows, cols):
    print(name, [[float(m[r, c]) for c in cols] for r in rows])


TAU = PHI + RHO
exp_tau = expm(hat4(TAU))
print("SKEW_X", [float(exp_tau[i, 3]) for i in range(3)])
jl = series(hat3(PHI))
show("SKEW_JL", jl, range(3), range(3))
show("SKEW_JL_INV", inverse(jl), range(3), range(3))
jl6 = series(ad(TAU))
show("SKEW_Q", jl6, range(3, 6), range(3))
show("SKEW_Q_INV", inverse(jl6), range(3, 6), range(3))
