import numpy as np
import matplotlib.pyplot as plt

n = 0
g = 9.81
M = 1
coef = -n / M


def vitesse_init(angle, vo):
    xp = vo * np.cos(np.deg2rad(angle))
    yp = vo * np.sin(np.deg2rad(angle))
    return xp, yp


def f1(xp):
    return coef * xp


def f2(yp):
    return (coef * yp) - g


def Euleur(dt, x_init, y_init, xp_init, yp_init):
    x = x_init
    y = y_init
    xp = xp_init
    yp = yp_init
    t = 0

    tab_t = np.array(0)
    tab_x = np.array(x)
    tab_y = np.array(y)
    tab_xp = np.array(xp)
    tab_yp = np.array(yp)

    while (y >= 0):
        t += dt
        x += xp * dt
        xp += f1(xp) * dt
        y += yp * dt
        yp += f2(yp) * dt
        tab_x = np.append(tab_x, x)
        tab_y = np.append(tab_y, y)
        tab_xp = np.append(tab_xp, xp)
        tab_yp = np.append(tab_yp, yp)
        tab_t = np.append(tab_t, t)
    return tab_x, tab_y, tab_xp, tab_yp, tab_t


def Euleur_Cromer(dt,x_init,y_init,xp_init,yp_init):
    x = x_init
    y = y_init
    xp = xp_init
    yp = yp_init
    t = 0

    tab_t = np.array(0)
    tab_x = np.array(x)
    tab_y = np.array(y)
    tab_xp = np.array(xp)
    tab_yp = np.array(yp)

    while (y >= 0):
        t += dt
        xp += f1(xp) * dt
        x += xp * dt
        yp += f2(yp) * dt
        y += yp * dt
        tab_x = np.append(tab_x, x)
        tab_y = np.append(tab_y, y)
        tab_xp = np.append(tab_xp, xp)
        tab_yp = np.append(tab_yp, yp)
        tab_t = np.append(tab_t, t)
    return tab_x, tab_y, tab_xp, tab_yp, tab_t




xp, yp = vitesse_init(45, 0.5)
x, y, xp1, yp1, tab_t = Euleur_Cromer(0.001, 0, 0, xp, yp)
x2,y2,xp2,yp2,tab_t2  = Euleur(0.001, 0, 0, xp, yp)
plt.plot(tab_t, y,label='cromer')
plt.plot(tab_t2,y2,label='euleur')
plt.legend()
plt.show()
