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

        # write found line to file
        for j in range(0, m-1):
            DATA.write(str(raw_data[i, j]) + ' ')
        DATA.write(str(raw_data[i, m-1]))
        DATA.write('\n')

        t = t + tau
DATA.close()
