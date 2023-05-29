import numpy as np
import mpmath

mpmath.mp.dps = 50

g = mpmath.mpf(9.81)

N = 3 #Number of stages
L = np.ones(N) #Vector of pendulum lengths
M = np.ones(N) #Vector of masses
Th = np.zeros((N,1)) #Vector of angles
Om = np.zeros((N,1)) #Vector of angular velocities

for i in range(N):
    # input
    L[i] = mpmath.mpf(1)
    M[i] = mpmath.mpf(1)
    Th[i][0] = mpmath.mpf(np.pi/2)
    Om[i][0] = mpmath.mpf(0)


def F(Th, Om):
    '''
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Vector of functions F(Th, Om)
    '''
    A = np.zeros((N, N))
    B = np.zeros((N, N))
    G = np.zeros((N,1))
    Y = Om ** 2

    for k in range(N):
        G[k][0] = -g * np.sin(Th[k][0]) * np.sum(M[k:])
        for i in range(N):
            if i > k:
                mass_sum = np.sum(M[i:])
            else:
                mass_sum = np.sum((M[k:]))
            A[k][i] = L[i] * np.cos(Th[i][0] - Th[k][0]) * mass_sum
            B[k][i] = L[i] * np.sin(Th[i][0] - Th[k][0]) * mass_sum

    return np.linalg.inv(A) * (np.matmul(B, Y) + G)

def Runge_Kutta(Th, Om, h):
    '''
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :param h: Step
    :return: Vectors of angles and velocities after step
    '''

    K1_Th = Om
    K1_Om = F(Th, Om)
    K2_Th = Om + 0.5 * h * K1_Om
    K2_Om = F(Th + 0.5 * h * K1_Th, Om + 0.5 * h * K1_Om)
    K3_Th = Om + 0.5 * h * K2_Om
    K3_Om = F(Th + 0.5 * h * K2_Th, Om + 0.5 * h * K2_Om)
    K4_Th = Om + h * K3_Om
    K4_Om = F(Th + h * K3_Th, Om + h * K3_Om)

    return Th + mpmath.mpf(1/6) * h * (K1_Th + 0.5 * K2_Th + 0.5 * K3_Th + K4_Th), \
           Om + mpmath.mpf(1/6) * h * (K1_Om + 0.5 * K2_Om + 0.5 * K3_Om + K4_Om)
