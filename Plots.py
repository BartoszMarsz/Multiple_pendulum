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
X = np.zeros((tau,1))
Y = np.zeros((tau,1))

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

X_k, Y_k = XY_pos(3)

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

def animation():
    # defining figure and plots
    fig = plt.figure(figsize=(15, 8))
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 3)
    ax3 = fig.add_subplot(2, 2, (2, 4))
    ax3.set_title('Motion of double pendulum', fontsize=20)
    ax1.set_title(r'Graph of $\theta_1(t)$ and $\theta_2(t)$', fontsize=20)
    ax2.set_title(r'Graph of $\omega_1(t)$ and $\omega_2(t)$', fontsize=20)
    # defining axises
#    ax1.set_xlim(0, time)
#    ax1.set_ylim(th_min-0.5, th_max+0.5)
#    ax2.set_xlim(0, t)
#    ax2.set_ylim(om_min-0.5, om_max+0.5)
    ax3.set_xlim(-1.1 * np.sum(L), 1.1 * np.sum(L))
    ax3.set_ylim(-1.1 * np.sum(L), 1.1 * np.sum(L))

    ax1.set_ylabel('[rad]', fontsize=15)
    ax2.set_ylabel(r'[$\frac{rad}{s}$]', fontsize=15)
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

    # theta and omega (graphs)
#    TH1, = ax1.plot(0, 0)
#    TH2, = ax1.plot(0, 0)
#    OM1, = ax2.plot(0, 0)
#    OM2, = ax2.plot(0, 0)
    # pendulum (motion)
#    trajectory1, = ax3.plot(0, 0)
#    trajectory2, = ax3.plot(0, 0)
#    line1, = ax3.plot(0, 0)
#    line2, = ax3.plot(0, 0)
    line, = ax3.plot(0, 0)
    # timer
    timer = ax3.text(0.05, 0.95, '', transform=ax3.transAxes, verticalalignment='top', color='black',fontsize=20, bbox=(dict(facecolor='wheat', boxstyle='round')))

    # specifying objects appearance
#    TH1.set_color('darkorange')
#    TH2.set_color('mediumblue')
#    OM1.set_color('darkorange')
#    OM2.set_color('mediumblue')

#    trajectory1.set_color('lightgray')
#    trajectory1.set_alpha(0.3)
#    trajectory2.set_color('lightgray')
#    trajectory2.set_alpha(0.3)
#    line1.set_color('darkorange')
#    line2.set_color('mediumblue')
#    line1.set_linewidth(2)
#    line2.set_linewidth(2)

    # defining legends
#    ax1.legend((TH1, TH2), (r'$\theta_1$', r'$\theta_2$'), loc='upper right', shadow=True, labelcolor='white',
#               facecolor='black', fontsize=20)
#    ax2.legend((OM1, OM2), (r'$\omega_1$', r'$\omega_2$'), loc='upper right', shadow=True, labelcolor='white',
#               facecolor='black', fontsize=20)

    def animation_frame(i):
#        TH1.set_xdata(time[:i])
#        TH1.set_ydata(th1[:i])
#        TH2.set_xdata(time[:i])
#        TH2.set_ydata(th2[:i])
#        OM1.set_xdata(time[:i])
#        OM1.set_ydata(om1[:i])
#        OM2.set_xdata(time[:i])
#        OM2.set_ydata(om2[:i])

#        line1.set_xdata([0.0, X1[i]])
#        line1.set_ydata([0.0, Y1[i]])
#        line2.set_xdata([X1[i], X2[i]])
#        line2.set_ydata([Y1[i], Y2[i]])
#        trajectory1.set_xdata(X1[:i])
#        trajectory1.set_ydata(Y1[:i])
#        trajectory2.set_xdata(X2[:i])
#        trajectory2.set_ydata(Y2[:i])
        line.set_xdata(X[i,:])
        line.set_ydata(Y[i,:])
#        timer.set_text(str(time[i]))

        return line,

    anim = FuncAnimation(fig, func=animation_frame, frames=range(0,tau), interval=h*1000, repeat=False, blit=True)
    plt.show()

animation()