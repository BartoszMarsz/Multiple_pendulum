import numpy as np
import mpmath as mp

mp.mp.dps = 50

g = mp.mpf(9.81)

h = mp.mpf(0.1) #jump
N = 3 #Number of stages
L = mp.zeros(N, 1) #Vector of pendulum lengths
M = mp.zeros(N, 1) #Vector of masses
Th = mp.zeros(N, 1) #Vector of angles
Om = mp.zeros(N, 1) #Vector of angular velocities

for i in range(N):
    # input
    L[i] = mp.mpf(1.0)
    M[i] = mp.mpf(1.0)
    Th[i, 0] = mp.mpf(np.pi/2)
    Om[i, 0] = mp.mpf(0.0)


def F(Th, Om):
    '''
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Vector of functions F(Th, Om)
    '''
    A = mp.zeros(N)
    B = mp.zeros(N)
    G = mp.zeros(N, 1)
    Y = mp.zeros(N, 1)
    for i in range(N):
        Y[i] = Om[i] ** 2

    for k in range(N):
        mass_sum = mp.mpf(0)
        for n in range(k,N):
            mass_sum = mass_sum + M[n]
        G[k, 0] = -g * mp.sin(Th[k, 0]) * mass_sum
        for i in range(N):
            if i > k:
                mass_sum = mp.mpf(0)
                for n in range(i, N):
                    mass_sum = mass_sum + M[n]
            A[k, i] = L[i] * mp.cos(Th[i, 0] - Th[k, 0]) * mass_sum
            B[k, i] = L[i] * mp.sin(Th[i, 0] - Th[k, 0]) * mass_sum

    return (A ** -1) * ((B * Y) + G)

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

    return Th + mp.mpf(1/6) * h * (K1_Th + 0.5 * K2_Th + 0.5 * K3_Th + K4_Th), \
           Om + mp.mpf(1/6) * h * (K1_Om + 0.5 * K2_Om + 0.5 * K3_Om + K4_Om)

mp.nprint(Th, 10)
mp.nprint(Om, 10)
Th, Om = Runge_Kutta(Th, Om, h)
Th, Om = Runge_Kutta(Th, Om, h)
Th, Om = Runge_Kutta(Th, Om, h)
Th, Om = Runge_Kutta(Th, Om, h)
mp.nprint(Th, 10)
mp.nprint(Om, 10)
