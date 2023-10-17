import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
params = np.loadtxt('params.pdb', delimiter=' ', usecols=1)
time = float(params[0])
h = float(params[1])
N = int(params[2])
L = params[3:3+N]
M = params[3+N:3+2*N]
tau = int(time/h)+1


data = np.loadtxt('data.pdb', delimiter=' ')
iterator = data[:, 0]
Th = data[:, 1:1 + N]
Om = data[:, 1+N:1+2*N]
T_data = data[:, 1+2*N]
V_data = data[:, 2+2*N]
dE_data = data[:, 3+2*N]
X = np.zeros((tau,1))
Y = np.zeros((tau,1))

E_max = max(10,np.amax(np.concatenate((T_data,V_data))))
E_min = min(-10,np.amin(np.concatenate((T_data,V_data))))
dE_max = max(10,np.amax(dE_data))

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

def XY_pos(k):
    X_k = np.sum(L[0:k]*np.sin(Th[:, 0:k]), axis=1)
    Y_k = - np.sum(L[0:k]*np.cos(Th[:, 0:k]), axis=1)
    return X_k, Y_k

for i in range(1,N+1):
    X = np.concatenate((X, np.sum(L[0:i]*np.sin(Th[:, 0:i]), axis=1, keepdims=True)), axis=1)
    Y = np.concatenate((Y, - np.sum(L[0:i]*np.cos(Th[:, 0:i]), axis=1, keepdims=True)), axis=1)


X_k, Y_k = XY_pos(N)

fig, ax = plt.subplots()
plt.title('Trajectory of ' + str(N) + '-mass', fontsize=15)
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

def animation():
    # defining figure and plots
    fig = plt.figure(figsize=(15, 8))
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 3)
    ax3 = fig.add_subplot(2, 2, (2, 4))
    ax3.set_title('Motion of multiple pendulum', fontsize=20)
    ax1.set_title('Graph of kinetic and potential energy', fontsize=20)
    ax2.set_title(r'Graph of difference in energy', fontsize=20)
    # defining axises
    ax1.set_xlim(0, time)
    ax1.set_ylim(1.1 * E_min, 1.1 * E_max)
    ax2.set_xlim(0, time)
    ax2.set_ylim(0, 1.1 * dE_max)
    ax3.set_xlim(-1.1 * np.sum(L), 1.1 * np.sum(L))
    ax3.set_ylim(-1.1 * np.sum(L), 1.1 * np.sum(L))

    ax1.set_ylabel('[J]', fontsize=15)
    ax2.set_ylabel('[%]', fontsize=15)
    ax2.set_xlabel('t[s]', fontsize=15)
    ax3.set_xlabel('[m]', fontsize=15)
    ax3.set_ylabel('[m]', fontsize=15)

    # specifying appearance
    ax1.grid(color='dimgrey')
    ax1.set_facecolor(color='black')
    ax2.grid(color='dimgrey')
    ax2.set_facecolor(color='black')
    ax3.grid(color='dimgrey')
    ax3.set_facecolor(color='black')

    ratio = 1.0
    x_left, x_right = ax3.get_xlim()
    y_low, y_high = ax3.get_ylim()
    ax3.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)

    # defining objects

    # Energy (graphs)
    T, = ax1.plot(0, 0)
    V, = ax1.plot(0, 0)
    dE, = ax2.plot(0, 0)

    # pendulum (motion)
    trajectory, = ax3.plot(0, 0)
#    trajectory2, = ax3.plot(0, 0)
#    line1, = ax3.plot(0, 0)
#    line2, = ax3.plot(0, 0)
    line, = ax3.plot(0, 0)
    # timer
#    timer, = ax3.text(0.05, 0.95, '', transform=ax3.transAxes, verticalalignment='top', color='black',fontsize=20, bbox=(dict(facecolor='wheat', boxstyle='round')))

    # specifying objects appearance
    T.set_color('darkorange')
    V.set_color('mediumblue')
    dE.set_color('darkorange')


    trajectory.set_color('lightgray')
    trajectory.set_alpha(0.3)
#    trajectory2.set_color('lightgray')
#    trajectory2.set_alpha(0.3)
#    line1.set_color('darkorange')
#    line2.set_color('mediumblue')
#    line1.set_linewidth(2)
#    line2.set_linewidth(2)

    # defining legends
    ax1.legend((T, V), ('Kinetic energy', 'Potential energy'), loc='upper right', shadow=True, labelcolor='white',
               facecolor='black', fontsize=20)
#    ax2.legend(dE, r'$frac{\DeltaE}{E}$', loc='upper right', shadow=True, labelcolor='white',
#               facecolor='black', fontsize=20)

    def animation_frame(i):
        T.set_xdata(iterator[:i]*time/(tau-1))
        T.set_ydata(T_data[:i])
        V.set_xdata(iterator[:i]*time/(tau-1))
        V.set_ydata(V_data[:i])
        dE.set_xdata(iterator[:i]*time/(tau-1))
        dE.set_ydata(dE_data[:i])

#        line1.set_xdata([0.0, X1[i]])
#        line1.set_ydata([0.0, Y1[i]])
#        line2.set_xdata([X1[i], X2[i]])
#        line2.set_ydata([Y1[i], Y2[i]])
        trajectory.set_xdata(X[:i,N])
        trajectory.set_ydata(Y[:i,N])
#        trajectory2.set_xdata(X2[:i])
#        trajectory2.set_ydata(Y2[:i])
        line.set_xdata(X[i,:])
        line.set_ydata(Y[i,:])
#        timer.set_text(str(iterator[:i]/(tau-1)))

        return line, T, V, dE, trajectory,

    anim = FuncAnimation(fig, func=animation_frame, frames=range(0,tau), interval=h*1000, repeat=False, blit=True)
    plt.show()

animation()