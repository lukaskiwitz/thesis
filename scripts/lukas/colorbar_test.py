import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

x = np.linspace(-10,10,100)

X,Y = np.meshgrid(x,x)

Z = np.sin(X*Y)

norm_cytokine = Normalize(vmin=0, vmax=2)
mappable = ScalarMappable(norm=norm_cytokine, cmap="viridis")

fig = plt.figure()
cb_cytokine = fig.colorbar(mappable)

plt.contourf(Z,norm = norm_cytokine)
plt.show()

