import mpmath as mp

mp.mp.dps = 50  # number of digits

# Simulation parameters
g = mp.mpf('9.81')  # standard gravity !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
time = 1  # time of simulation in seconds
h0 = mp.mpf('0.01')  # default time step
eta = 0.01  # precision of energy conservation[%]
N = 5  # number of pendulum stages
t = mp.mpf('0')

# Declaration of pendulum stages parameters
L = mp.matrix(N, 1)  # vector of pendulum lengths (2-dim, vertical)
M = mp.matrix(N, 1)  # vector of masses (2-dim, vertical)
Th = mp.matrix(N, 1)  # vector of angles (2-dim, vertical)
Om = mp.matrix(N, 1)  # vector of angular velocities (2-dim, vertical)

# Definition of pendulum stages parameters
for i in range(N):
    L[i, 0] = mp.mpf(1.0) / mp.mpf(N)
    M[i, 0] = mp.mpf(1.0) / mp.mpf(N)
    Th[i, 0] = mp.mpf(mp.radians(90 + (90 / (N - 1)) * i))
    Om[i, 0] = mp.mpf(0)


def T_energy(Th, Om):
    """
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Kinetic energy
    """
    T = mp.mpf('0')
    for n in range(1, N + 1):
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                T = T + mp.mpf('0.5') * M[n - 1, 0] * L[i - 1, 0] * L[j - 1, 0] * Om[i - 1, 0] * Om[j - 1, 0] * mp.cos(
                    Th[i - 1, 0] - Th[j - 1, 0])
    return T


def V_energy(Th, Om):
    """
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Potential energy
    """
    V = mp.mpf('0')
    for n in range(1, N + 1):
        for i in range(1, n + 1):
            V = V - g * M[n - 1, 0] * L[i - 1, 0] * mp.cos(Th[i - 1, 0])
    return V


# Initial energy
T0 = T_energy(Th, Om)  # kinetic energy
V0 = V_energy(Th, Om)  # potential energy
E0 = T0 + V0  # total energy


def diff_energy(Th, Om):
    """
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :return: Difference in energy/E0
    """
    return abs(100 * ((V_energy(Th, Om) + T_energy(Th, Om) - E0) / E0))


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
    # eq (2.25)
    for i in range(N):
        Om_2[i, 0] = Om[i, 0] ** 2

    for k in range(N):
        mass_sum = mp.mpf('0')
        for n in range(k, N):
            mass_sum = mass_sum + M[n, 0]
        # eq (2.25)
        C[k, 0] = -g * mp.sin(Th[k, 0]) * mass_sum
        for i in range(N):
            if i > k:
                mass_sum = mp.mpf('0')
                for n in range(i, N):
                    mass_sum = mass_sum + M[n, 0]
            # eq (2.26), (2.27)
            A[k, i] = L[i, 0] * mp.cos(Th[i, 0] - Th[k, 0]) * mass_sum
            B[k, i] = L[i, 0] * mp.sin(Th[i, 0] - Th[k, 0]) * mass_sum

    return (A ** -1) * ((B * Om_2) + C)  # eq (2.24)


def Runge_Kutta(Th, Om, h):
    """
    :param Th: Vector of angles
    :param Om: Vector of angular velocities
    :param h: Step
    :return: Vectors of angles and velocities after step
    """
    # eq (3.17)-(3.24)
    K1_Th = Om
    K1_Om = F(Th, Om)
    K2_Th = Om + mp.mpf('0.5') * h * K1_Om
    K2_Om = F(Th + mp.mpf('0.5') * h * K1_Th, Om + mp.mpf('0.5') * h * K1_Om)
    K3_Th = Om + mp.mpf('0.5') * h * K2_Om
    K3_Om = F(Th + mp.mpf('0.5') * h * K2_Th, Om + mp.mpf('0.5') * h * K2_Om)
    K4_Th = Om + h * K3_Om
    K4_Om = F(Th + h * K3_Th, Om + h * K3_Om)

    # eq (3.25), (3.26)
    return Th + mp.mpf('1/6') * h * (K1_Th + mp.mpf('2') * K2_Th +
                                     mp.mpf('2') * K3_Th + K4_Th), \
           Om + mp.mpf('1/6') * h * (K1_Om + mp.mpf('2') * K2_Om +
                                     mp.mpf('2') * K3_Om + K4_Om)


# Save parameters to file
DATA = open('params.pdb', 'w')
DATA.write('Time= ' + str(time) + '\n')
DATA.write('Step= ' + str(h0) + '\n')
DATA.write('Number_of_stages= ' + str(N) + '\n')
for i in range(N):
    DATA.write('Length_of_' + str(i + 1) + '_stage= ' + str(L[i, 0]) + '\n')
for i in range(N):
    DATA.write('Mass_of_' + str(i + 1) + '_mass= ' + str(M[i, 0]) + '\n')
DATA.close()

# Save initial data
DATA = open('raw_data.pdb', 'w')
DATA.write(str(0))
for i in range(N):
    DATA.write(' ' + str(round(Th[i, 0], 5)))
for i in range(N):
    DATA.write(' ' + str(round(Om[i, 0], 5)))
DATA.write(' ' + str(round(T_energy(Th, Om), 5)))
DATA.write(' ' + str(round(V_energy(Th, Om), 5)))
DATA.write(' ' + str(diff_energy(Th, Om)))
DATA.write(' ' + str(h0))
DATA.write('\n')

# Simulation loop
while t < time + h0:
    print('t = ' + str(round(t, 5)), end='    ')

    # Calculation of angles and angular velocities after time step
    h = h0  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Th_1, Om_1 = Runge_Kutta(Th, Om, h)
    dE = diff_energy(Th_1, Om_1)

    # Condition of energy conservation
    while dE > eta:
        h = mp.mpf('0.1') * h
        Th_1, Om_1 = Runge_Kutta(Th, Om, h)
        dE = diff_energy(Th_1, Om_1)  # Obliczenie tylko dE!!!!!!!!!!!!!!!!!!
    print('dE = ' + str(round(dE, 5)), end='    ')
    print('h = ' + str(h))
    Th, Om = Th_1, Om_1
    T = T_energy(Th, Om)  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    V = V_energy(Th, Om)  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    # Write data to file
    DATA.write(str(round(t, 8)))
    for j in range(N):
        DATA.write(' ' + str(round(Th[j, 0], 5)))
    for j in range(N):
        DATA.write(' ' + str(round(Om[j, 0], 5)))
    DATA.write(' ' + str(round(T, 5)))
    DATA.write(' ' + str(round(V, 5)))
    DATA.write(' ' + str(dE))
    DATA.write(' ' + str(h))
    DATA.write('\n')

    # Time iteration
    t = t + h
DATA.close()
