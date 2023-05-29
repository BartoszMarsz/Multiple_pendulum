from typing import Optional, Any

import numpy as np
import mpmath

mpmath.mp.dps = 50

g = mpmath.mpf(9.81)

N = 3
L = np.zeros(N)
M = np.zeros(N)
Th = np.zeros(N)
Om = np.zeros(N)

for i in range(N):
    # input
    L[i] = mpmath.mpf(1)
    M[i] = mpmath.mpf(1 / 3)
    Th[i] = mpmath.mpf(0)
    Om[i] = mpmath.mpf(0)


def F(Th, Om):
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    G = np.zeros(N)
    Y = Om ** 2

    for k in range(N):
        G[k] = -g * np.sin(Th[k]) * np.sum(M[k:])
        for i in range(N):
            if i > k:
                mass_sum = np.sum(M[i:])
            else:
                mass_sum = np.sum((M[k:]))
            A[k][i] = L[i] * np.cos(Th[i] - Th[k]) * mass_sum
            B[k][i] = L[i] * np.sin(Th[i] - Th[k]) * mass_sum

    return np.linalg.inv(A) * (np.matmul(B, Y) + G)
