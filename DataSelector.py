import numpy as np
import math
raw_data = np.loadtxt('raw_data.pdb', delimiter=' ')
t_max = raw_data.shape[0]
t = 0
DATA = open('data.pdb', 'w')
for i in range(t_max):
#    print('t % 0.01 = '+str(round(raw_data[i, 0] % 0.01,2)))
    if(raw_data[i, 0]>=t):
        for j in range(raw_data.shape[1]-1):
            DATA.write(str(math.floor(raw_data[i, j]*100)/100)+' ')
        DATA.write(str(math.floor(raw_data[i, raw_data.shape[1]-1] * 100) / 100))
        DATA.write('\n')
        t = t + 0.01