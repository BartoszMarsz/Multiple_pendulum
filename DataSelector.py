import numpy as np
import math
raw_data = np.loadtxt('raw_data.pdb', delimiter=' ')
t_max = raw_data.shape[0]
t = 0
DATA = open('data.pdb', 'w')
for i in range(t_max):
    if(raw_data[i, 0]>=round(t,3)):
        print(round(t,3), end='    ')
        print(raw_data[i, 0])
        DATA.write(str(math.floor(raw_data[i, 0] * 1000) / 1000) + ' ')
        for j in range(1,raw_data.shape[1]-4):
            DATA.write(str(round(raw_data[i, j], 2)) + ' ')
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 4],3)) + ' ')
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 3],3)) + ' ')
        DATA.write(str(round(raw_data[i, raw_data.shape[1] - 2],6)) + ' ')
        DATA.write(str(raw_data[i, raw_data.shape[1] - 1]))
        DATA.write('\n')
        t = t + 1/60