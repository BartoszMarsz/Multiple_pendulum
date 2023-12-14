import numpy as np
import math

# load data to RAM
raw_data = np.loadtxt('raw_data.pdb', delimiter=' ')

# number of lines
n = raw_data.shape[0]

# frames per second
fps = 60
h = 1 / fps

t = 0
DATA = open('data.pdb', 'w')

# loop on every row
for i in range(n):
    # find iteration of h or nearest approximation
    if raw_data[i, 0] >= round(t, 3):
        print(round(t,3), end='    ')
        print(raw_data[i, 0])

        # write found line to file

        # write time (rounded down to 3 digits)
        DATA.write(str(math.floor(raw_data[i, 0] * 1000) / 1000) + ' ')

        # write angles and angular velocities (rounded to 2 digits)
        for j in range(1,raw_data.shape[1]-4):
            DATA.write(str(round(raw_data[i, j], 2)) + ' ')

        # write T, V (3 digits) and dE (6 digits)
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 4],3)) + ' ')
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 3],3)) + ' ')
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 2],6)) + ' ')

        # write h
        DATA.write(str(raw_data[i, raw_data.shape[1] - 1]))
        DATA.write('\n')

        t = t + h
