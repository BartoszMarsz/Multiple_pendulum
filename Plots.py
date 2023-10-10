import numpy as np
import matplotlib.pyplot as plt

params = np.loadtxt('params.pdb', delimiter=' ', usecols=1)
time = float(params[0])
h = float(params[1])
N = int(params[2])
L = params[3:3+N]
M = params[3+N:3+2*N]


data = np.loadtxt('trajectory.pdb', delimiter=' ')
iterator = data[:, 0]
Th = data[:, 1:1 + N]
Om = data[:, 1+N:1+2*N]

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

def trajectory(k):
    X_k = np.sum(L[0:k]*np.sin(Th[:, 0:k]), axis=1)
    Y_k = - np.sum(L[0:k]*np.cos(Th[:, 0:k]), axis=1)
    return X_k, Y_k


X_k, Y_k = trajectory(3)

fig, ax = plt.subplots()
plt.title('Trajectory of ' + str(10) + '-mass', fontsize=15)
plt.plot(X_k, Y_k, color='orange', label='mass k')
plt.xticks(np.linspace(-1, 1, 3), fontsize=15)
plt.yticks(np.linspace(-1, 1, 3), fontsize=15)

ratio = 1.0
x_left, x_right = ax.get_xlim()
y_low, y_high = ax.get_ylim()
ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)

plt.legend(loc='upper right', fontsize=15)
plt.xlabel('x', fontsize=15)
plt.ylabel('y', fontsize=15)

ax.grid(color='dimgrey')
ax.set_facecolor(color='black')

figManager = plt.get_current_fig_manager()
figManager.window.showMaximized()
plt.show()
