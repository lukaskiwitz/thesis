import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors, DistanceMetric
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation



T = 100
n = 500
s = 500
dt = 0.01


class Sim:

    def __init__(self):
        np.random.seed(0)
        self.r = np.random.uniform(-2*s, 2*s, [2, n])
        self.c = np.random.uniform(-2*s, 2*s, [2, 10])


    def update(self,i):

        F = self.F(self.r, dt)
        self.r = self.r +dt*F
        ax = plt.gca()



        ax.clear()
        ax.set_xlim((-2*s, 2*s))
        ax.set_ylim((-2*s, 2*s))


        ax.scatter(self.r[0],self.r[1],s = 100)

        ax.scatter(self.c[0], self.c[1], s=1000,color = "red")

        # ax.quiver(self.r[0], self.r[1], dt*F[0], dt*F[1],units = "xy")
        return ax

    def F(self, r, dt):

        def get_b(D, r, range):

            D_r = np.where(
                (D >= range[0]) & (D < range[1])
                , 1, 0)

            def f(d, r):
                result = np.mean(r.T[np.ravel(np.argwhere(d))], axis=0)
                return result



            return np.nan_to_num(np.apply_along_axis(f, 1, D_r, r).T)

        def get_F(X, Y, a, b, c,attracted, repulse):

            dist = DistanceMetric.get_metric('euclidean')
            D = dist.pairwise(X, Y)

            r_rep = (a, b)
            r_attracted = (b, c)

            b_rep = get_b(D, Y.T, r_rep)
            b_attracted = get_b(D, Y.T, r_attracted)


            F = -attracted*np.where(b_attracted != 0,1,0)*(X.T - b_attracted) +  repulse*np.where(b_rep != 0,1,0)*(X.T - b_rep)

            return F, b_attracted


        F_self,a_self  = get_F(r.T,r.T,0,30,50,0,20)

        F_c, a_c = get_F(r.T,self.c.T,0,50,200,5,5)

        poison = np.random.poisson(1,a_c.shape)

        jump = np.where(np.linalg.norm(a_c) > 200,0,1000) * np.random.normal(0, 2, a_c.shape)

        diffuse = 10 * np.random.normal(0, 2, a_c.shape)

        F =  F_self + diffuse + poison*jump

        return F

sim = Sim()

fig = plt.figure(figsize=(10,10))
ani = FuncAnimation(fig, sim.update, frames=T)

ani.save("movie.mp4")
plt.show()


