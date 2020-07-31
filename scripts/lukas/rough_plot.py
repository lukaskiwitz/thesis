import matplotlib.pyplot as plt
import numpy as np

def f(x,n,K):

    return x**n/(K**n+x**n)

x_log = np.logspace(0,3,1000)

x = x_log#np.linspace(0,max(x_log),1000)
e = np.e


fig, ax = plt.subplots(2,2)

ax_l = ax[0][0]
ax_l.plot(x,f(x_log,1,4))
ax_l.plot(x,f(x_log,10,100))
ax_l.legend(["no Treg","1/3 Treg"])
ax_l.set_title("uniform Treg")
ax_l.set_xlabel("antigen dose")
ax_l.set_ylabel("tissue response")
ax_l.semilogx()

ax_l = ax[0][1]
ax_l.plot(x,f(x_log,10,4))
ax_l.plot(x,f(x_log,10,100))
ax_l.legend(["foreign antigen","self-antigen"])
ax_l.set_title("Treg Colocalization")
ax_l.set_xlabel("antigen dose")
ax_l.set_ylabel("tissue response")

ax_l.semilogx()


# plt.semilogx()
# plt.ylim([0,1])
plt.savefig("./rough_plot.pdf")
# plt.show()
