#!/usr/bin/env python3
"""
Morton-Taylor-Turner (1956) integral entrainment model for a turbulent round jet that
moves DOWN into linearly stratified (stable) air. Top-hat variables, s = distance below
the source:
  Q = pi b^2 w,  M = pi b^2 w^2,  F = Q g'        (volume, momentum, buoyancy fluxes)
  dQ/ds = 2 sqrt(pi) alpha sqrt(M)
  dM/ds = F Q / M                                   (F < 0 decelerates the jet)
  dF/ds = -N^2 Q
Source: the rotor wake, M0 = T/rho, Q0 = A v_i (actuator-disk momentum theory), F0 = 0
(the drone pushes down air from its own level). The jet stops where M -> 0.
For a pure jet the result scales as z_m = C(alpha) (M0/N^2)^(1/4).
"""
import math


def penetration(M0, Q0, N, alpha, ds=0.01):
    Q, M, F, s = Q0, M0, 0.0, 0.0
    while M > 1e-6 * M0:
        dQ = 2 * math.sqrt(math.pi) * alpha * math.sqrt(M)
        dM = F * Q / M
        dF = -N * N * Q
        Q, M, F, s = Q + dQ * ds, M + dM * ds, F + dF * ds, s + ds
        if s > 1e4:
            break
    return s


if __name__ == "__main__":
    rho, g, th0 = 1.25, 9.81, 285.0
    N = math.sqrt(g / th0 * 0.05)
    for name, T, A in (("heavy 50 kg", 50 * 9.81, 4 * 1.474), ("quad 2 kg", 2 * 9.81, 4 * 0.0507)):
        M0 = T / rho
        vi = math.sqrt(T / (2 * rho * A))
        Q0 = A * vi
        L = (M0 / N ** 2) ** 0.25
        print(name, f"M0={M0:.1f} m4/s2  N={N:.4f} 1/s  (M0/N^2)^1/4={L:.1f} m")
        for a in (0.05, 0.065, 0.076, 0.09, 0.1):
            z = penetration(M0, Q0, N, a)
            print(f"   alpha={a:.3f}  z_m={z:6.1f} m  C={z / L:.2f}")
