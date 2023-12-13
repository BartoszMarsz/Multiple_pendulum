from numpy import vectorize
import mpmath as mp
import matplotlib.pyplot as plt

mp.mp.dps = 50

g = mp.mpf('9.81')

time = 1 #Time in seconds
h = mp.mpf('0.01') #Step
precision = 0.01 #[%]
N = 5 #Number of stages
L = mp.matrix(N, 1) #Vector of pendulum lengths
M = mp.matrix(N, 1) #Vector of masses
Th = mp.matrix(N, 1) #Vector of angles
Om = mp.matrix(N, 1) #Vector of angular velocities

for i in range(N):
    # input
    L[i, 0] = mp.mpf(1.0)/mp.mpf(N)
    M[i, 0] = mp.mpf(1.0)/mp.mpf(N+4)
    Th[i, 0] = mp.mpf(mp.radians(90+(90/(N-1))*i))
    Om[i, 0] = mp.mpf(0)
M[N-1,0] = mp.mpf(5)/mp.mpf(N+4)

def sum_to(V, start, stop):
    """
    :param V: Vector to sum
    :param start: index start (included)
    :param stop: index of stop (included)
    :return: Sum of elements from start to stop
    """
    sum = 0
    for i in range(start, stop+1):
        sum = sum + V[i]
    return sum


def T_energy(Th, Om):
    T = mp.mpf('0')
    for n in range(1,N+1):
        for i in range(1,n+1):
            for j in range(1,n+1):
                T = T + mp.mpf('0.5') * M[n-1,0] * L[i-1,0] * L[j-1,0] * Om[i-1,0] * Om[j-1,0] * mp.cos(Th[i-1,0]-Th[j-1,0])
    return T


def V_energy(Th, Om):
    V = mp.mpf('0')
    for n in range(1,N+1):
        for i in range(1,n+1):
            V = V - g * M[n-1,0] * L[i-1,0] * mp.cos(Th[i-1,0])
    return V

T0 = T_energy(Th, Om)
V0 = V_energy(Th, Om)
E0 = T0+V0

def diff_energy(Th, Om):
    return abs(100*((V_energy(Th, Om)+T_energy(Th, Om)-E0)/E0))

def F(Th, Om):
    """
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Vector of functions F(Th, Om)
    """
    A = mp.zeros(N)
    B = mp.zeros(N)
    Om_2 = mp.zeros(N, 1)
    C = mp.zeros(N, 1)
    #eq (2.25)
    for i in range(N):
        Om_2[i,0] = Om[i,0] ** 2

    for k in range(N):
        mass_sum = mp.mpf('0')
        for n in range(k,N):
            mass_sum = mass_sum + M[n,0]
        #eq (2.25)
        C[k, 0] = -g * mp.sin(Th[k, 0]) * mass_sum
        for i in range(N):
            if i > k:
                mass_sum = mp.mpf('0')
                for n in range(i, N):
                    mass_sum = mass_sum + M[n,0]
            #eq (2.26), (2.27)
            A[k, i] = L[i, 0] * mp.cos(Th[i, 0] - Th[k, 0]) * mass_sum
            B[k, i] = L[i, 0] * mp.sin(Th[i, 0] - Th[k, 0]) * mass_sum

    return (A ** -1) * ((B * Om_2) + C) #eq (2.24)

def Runge_Kutta(Th, Om, h):
    '''
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :param h: Step
    :return: Vectors of angles and velocities after step
    '''
    #eq (3.17)-(3.24)
    K1_Th = Om
    K1_Om = F(Th, Om)
    K2_Th = Om + mp.mpf('0.5') * h * K1_Om
    K2_Om = F(Th + mp.mpf('0.5') * h * K1_Th, Om + mp.mpf('0.5') * h * K1_Om)
    K3_Th = Om + mp.mpf('0.5') * h * K2_Om
    K3_Om = F(Th + mp.mpf('0.5') * h * K2_Th, Om + mp.mpf('0.5') * h * K2_Om)
    K4_Th = Om + h * K3_Om
    K4_Om = F(Th + h * K3_Th, Om + h * K3_Om)

    #eq (3.25), (3.26)
    return Th + mp.mpf('1/6') * h * (K1_Th + mp.mpf('2') * K2_Th +
                                     mp.mpf('2') * K3_Th + K4_Th), \
           Om + mp.mpf('1/6') * h * (K1_Om + mp.mpf('2') * K2_Om +
                                     mp.mpf('2') * K3_Om + K4_Om)


DATA = open('raw_data.pdb', 'w')
#Initial parameters
DATA.write(str(0))
for i in range(N):
    DATA.write(' ' + str(round(Th[i,0],5)))
for i in range(N):
    DATA.write(' ' + str(round(Om[i,0],5)))
DATA.write(' ' + str(round(T_energy(Th, Om),5)))
DATA.write(' ' + str(round(V_energy(Th, Om),5)))
DATA.write(' ' + str(diff_energy(Th,Om)))
DATA.write(' ' + str(h))
DATA.write('\n')

t = 0
h_1 = mp.mpf('0')
while(t<time+h):
    print('t = '+str(round(t,5)), end='    ')
    Th_1, Om_1 = Runge_Kutta(Th, Om, h)
    h_1=h
    while(diff_energy(Th_1,Om_1)>precision):
        print('dE = ' + str(round(diff_energy(Th_1, Om_1), 5)), end='    ')
        h_1=mp.mpf('0.1')*h_1
        Th_1, Om_1 = Runge_Kutta(Th, Om, h_1)
    print('dE = ' + str(round(diff_energy(Th_1, Om_1), 5)), end='    ')
    print('h = '+str(h_1))
    Th, Om = Th_1, Om_1
    DATA.write(str(round(t,8)))
    for j in range(N):
        DATA.write(' ' + str(round(Th[j,0],5)))
    for j in range(N):
        DATA.write(' ' + str(round(Om[j,0],5)))
    DATA.write(' ' + str(round(T_energy(Th, Om),5)))
    DATA.write(' ' + str(round(V_energy(Th, Om),5)))
    DATA.write(' ' + str(abs(diff_energy(Th,Om))))
    DATA.write(' ' + str(h_1))
    DATA.write('\n')

    t = t+h_1

DATA.close()

DATA = open('params.pdb', 'w')
DATA.write('Time= ' + str(time) + '\n')
DATA.write('Step= ' + str(h) + '\n')
DATA.write('Number_of_stages= ' + str(N) + '\n')
for i in range(N):
    DATA.write('Length_of_' + str(i+1) + '_stage= ' + str(L[i,0]) + '\n')
for i in range(N):
    DATA.write('Mass_of_' + str(i+1) + '_mass= ' + str(M[i,0]) + '\n')
DATA.close()