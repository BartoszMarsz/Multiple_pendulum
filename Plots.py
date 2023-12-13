import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.animation import FFMpegWriter

params = np.loadtxt('5_stages_bigmass/params.pdb', delimiter=' ', usecols=1)
time = float(params[0])
h = float(params[1])
N = int(params[2])
L = params[3:3+N]
M = params[3+N:3+2*N]



data = np.loadtxt('5_stages_bigmass/data.pdb', delimiter=' ')
t = data[:, 0]
Th = data[:, 1:1 + N]
Om = data[:, 1+N:1+2*N]
T_data = data[:, 1+2*N]
V_data = data[:, 2+2*N]
dE_data = data[:, 3+2*N]
tau = t.shape[0]
X = np.zeros((tau,1))
Y = np.zeros((tau,1))

E_max = max(10,np.amax(np.concatenate((T_data,V_data))))
E_min = min(-10,np.amin(np.concatenate((T_data,V_data))))
dE_max = np.amax(dE_data)

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

def XY_pos(k):
    X_k = np.sum(L[0:k]*np.sin(Th[:, 0:k]), axis=1)
    Y_k = - np.sum(L[0:k]*np.cos(Th[:, 0:k]), axis=1)
    return X_k, Y_k

for i in range(1,N+1):
    X = np.concatenate((X, np.sum(L[0:i]*np.sin(Th[:, 0:i]), axis=1, keepdims=True)), axis=1)
    Y = np.concatenate((Y, - np.sum(L[0:i]*np.cos(Th[:, 0:i]), axis=1, keepdims=True)), axis=1)


X_k, Y_k = XY_pos(N)
def static():
    fig = plt.figure(figsize=(12.8, 6.8))
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 3)
    ax3 = fig.add_subplot(2, 2, (2, 4))
    ax3.set_title('Motion of multiple pendulum', fontsize=16)
    ax1.set_title('Graph of kinetic and potential energy', fontsize=16)
    ax2.set_title('Graph of difference in energy', fontsize=16)
    # defining axises
    ax1.set_xlim(0, time)
    ax1.set_ylim(1.1 * E_min, 1.1 * E_max)
    ax2.set_xlim(0, time)
    ax2.set_ylim(0, 0.011)
    ax3.set_xlim(-1.1 * np.sum(L), 1.1 * np.sum(L))
    ax3.set_ylim(-1.1 * np.sum(L), 1.1 * np.sum(L))

    ax1.set_ylabel('Energy [J]', fontsize=16)
    ax2.set_ylabel(r'$\delta E/E_0 [\%]$', rotation='horizontal', fontsize=16, ha='right')

    ax2.set_xlabel('t[s]', fontsize=16)
    ax3.set_xlabel('[m]', fontsize=16)
    ax3.set_ylabel('[m]', fontsize=16)

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
    line, = ax3.plot(0, 0)
    # specifying objects appearance
    T.set_color('darkorange')
    V.set_color('mediumblue')
    dE.set_color('darkorange')

    line.set_linewidth(2.5)
    trajectory.set_color('green')
    trajectory.set_alpha(0.9)

    # defining legends
    ax1.legend((T, V), ('Kinetic energy', 'Potential energy'), loc='upper right', shadow=True, labelcolor='white',
               facecolor='black', fontsize=10)

    # adding data
    T.set_data(t, T_data)
    V.set_data(t, V_data)
    dE.set_data(t, dE_data)
    trajectory.set_data(X[:, N], Y[:, N])
    line.set_data(X[0, :], Y[0, :])

    plt.savefig('5_stages_bigmass/5_stages_bigmass.png', bbox_inches='tight')
    plt.show()

def animation():
    # defining figure and plots
    fig = plt.figure(figsize=(12.8, 6.8))
    figManager = plt.get_current_fig_manager()
    #figManager.window.showMaximized()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 3)
    ax3 = fig.add_subplot(2, 2, (2, 4))
    ax3.set_title('Motion of multiple pendulum', fontsize=12)
    ax1.set_title('Graph of kinetic and potential energy', fontsize=12)
    ax2.set_title(r'Graph of difference in energy', fontsize=12)
    # defining axises
    ax1.set_xlim(0, time)
    ax1.set_ylim(1.1 * E_min, 1.1 * E_max)
    ax2.set_xlim(0, time)
    ax2.set_ylim(0, 1.1 * dE_max)
    ax3.set_xlim(-1.1 * np.sum(L), 1.1 * np.sum(L))
    ax3.set_ylim(-1.1 * np.sum(L), 1.1 * np.sum(L))

    ax1.set_ylabel('Energy [J]', fontsize=12)
    ax2.set_ylabel(r'$\delta E/E_0 [\%]$', rotation='horizontal', fontsize=12, ha='right')

    ax2.set_xlabel('t[s]', fontsize=12)
    ax3.set_xlabel('[m]', fontsize=12)
    ax3.set_ylabel('[m]', fontsize=12)

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
    line, = ax3.plot(0, 0)
    # timer
    timer = ax3.text(0.05, 0.95, '', transform=ax3.transAxes, va='top', ha='left', color='black',fontsize=14, bbox=(dict(facecolor='wheat', boxstyle='round')))

    # specifying objects appearance
    T.set_color('darkorange')
    V.set_color('mediumblue')
    dE.set_color('darkorange')


    trajectory.set_color('lightgray')
    trajectory.set_alpha(0.3)

    # defining legends
    ax1.legend((T, V), ('Kinetic energy', 'Potential energy'), loc='upper right', shadow=True, labelcolor='white',
               facecolor='black', fontsize=10)

    def init():
        T.set_data([t[0]], [T_data[0]])
        V.set_data([t[0]], [V_data[0]])
        dE.set_data([t[0]], [dE_data[0]])
        trajectory.set_data([X[0, N]],[Y[0, N]])
        line.set_xdata(X[0, :])
        line.set_ydata(Y[0, :])
        timer.set_text(str(0)+'s')
        return line, T, V, dE, trajectory, timer,
    def animation_frame(i):
        print(round(i/60,2))
        T.set_data(t[:i], T_data[:i])
        V.set_data(t[:i], V_data[:i])
        dE.set_data(t[:i], dE_data[:i])
        trajectory.set_data(X[:i,N], Y[:i,N])
        line.set_data(X[i,:],Y[i,:])
        timer.set_text(format(t[i], '.2f'))

        return line, T, V, dE, trajectory, timer,

    anim = FuncAnimation(fig, func=animation_frame, frames=tau, init_func=init, interval=int(1000/60), repeat=False, blit=True)
    writervideo = FFMpegWriter(fps=60)
    anim.save('5_stages_bigmass/Pendulum5.mp4', writer=writervideo)
    #plt.show()
    plt.close()
#static()
#animation()

