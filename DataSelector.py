import numpy as np

raw_data = np.loadtxt('raw_data.pdb', delimiter=' ')  # load data to RAM
fps = 60  # frames per second
n = raw_data.shape[0]  # number of rows
m = raw_data.shape[1]  # number of columns
tau = 1 / fps  # time step
t = 0  # time iterator

DATA = open('data.pdb', 'w')

# loop on every row
for i in range(n):
    # find iteration of tau or nearest approximation
    if raw_data[i, 0] >= round(t, 3):
        print(round(t, 3), end='    ')
        print(raw_data[i, 0])

        # write found line to file
        for j in range(0, m-1):
            DATA.write(str(raw_data[i, j]) + ' ')
        DATA.write(str(raw_data[i, m-1]))
        DATA.write('\n')

        t = t + tau
DATA.close()

"""
        # write angles and angular velocities (rounded to 2 digits)
        for j in range(1, m-4):
            DATA.write(str(round(raw_data[i, j], 2)) + ' ')

        # write T, V (3 digits) and dE (6 digits)
        DATA.write(str(round(raw_data[i, m - 4], 3)) + ' ')
        DATA.write(str(round(raw_data[i, m - 3], 3)) + ' ')
        DATA.write(str(round(raw_data[i, m - 2], 6)) + ' ')

        # write h
        DATA.write(str(raw_data[i, m - 1]))
"""