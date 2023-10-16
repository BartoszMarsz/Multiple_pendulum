from numpy import vectorize
import mpmath as mp
import matplotlib.pyplot as plt

mp.mp.dps = 50

g = mp.mpf(9.81)

time = 5
h = mp.mpf(1)/mp.mpf(50) #step
N = 10 #Number of stages
L = mp.matrix(N, 1) #Vector of pendulum lengths
M = mp.matrix(N, 1) #Vector of masses
Th = mp.matrix(N, 1) #Vector of angles
Om = mp.matrix(N, 1) #Vector of angular velocities

for i in range(N):
    # input
    L[i] = mp.mpf(1.0/mp.mpf(N))
    M[i] = mp.mpf(1.0/mp.mpf(N))
    Th[i, 0] = mp.mpf(mp.radians(120))
    Om[i, 0] = mp.mpf(0.0)


def sum_to(V, start, stop):
    '''
    :param V: Vector to sum
    :param start: index start (included)
    :param stop: index of stop (included)
    :return: Sum of elements from start to stop
    '''
    sum = 0
    for i in range(start, stop+1):
        sum = sum + V[i]
    return sum


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


DATA = open('data.pdb', 'w')
DATA.write(str(0))
for j in range(N):
    DATA.write(' ' + str(round(Th[j],5)))
for j in range(N):
    DATA.write(' ' + str(round(Om[j],5)))
DATA.write('\n')
for i in range(1, int(time/h+1)):
    DATA.write(str(i))
    Th, Om = Runge_Kutta(Th, Om, h)
    for j in range(N):
        DATA.write(' ' + str(round(Th[j],5)))
    for j in range(N):
        DATA.write(' ' + str(round(Om[j],5)))
    DATA.write('\n')
DATA.close()

DATA = open('params.pdb', 'w')
DATA.write('Time= ' + str(time) + '\n')
DATA.write('Step= ' + str(h) + '\n')
DATA.write('Number_of_stages= ' + str(N) + '\n')
for i in range(N):
    DATA.write('Length_of_' + str(i+1) + '_stage= ' + str(L[i]) + '\n')
for i in range(N):
    DATA.write('Mass_of_' + str(i+1) + '_mass= ' + str(M[i]) + '\n')
DATA.close()