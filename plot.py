#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter


#plot implied vol
x,y,z = np.loadtxt("impliedvol.txt", usecols=(0, 1, 2), unpack=True)

x = np.reshape(x, (-1, 50))
y = np.reshape(y, (-1, 50))
z = np.reshape(z, (-1, 50))


fig = plt.figure()
ax = plt.axes(projection='3d')

surf = ax.plot_surface(x, y, z, linewidth=1, rstride=2, cstride=2, cmap=plt.cm.RdPu, antialiased=True )

ax.zaxis.set_major_locator(LinearLocator(6))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_title('Implied Volatility Surface')        # title

vals = ax.get_zticks()
ax.set_zticklabels(['   {:3.2f}%'.format(x*100) for x in vals])

ax.view_init(elev=30,azim=150) 
fig.colorbar(surf, shrink=0.5, aspect=5)     # colour bar


#ax.plot_wireframe(x, y, z) #, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
plt.show()


#plot local vol
x,y,z = np.loadtxt("localvol.txt", usecols=(0, 1, 2), unpack=True)

x = np.reshape(x, (-1, 50))
y = np.reshape(y, (-1, 50))
z = np.reshape(z, (-1, 50))


fig = plt.figure()
ax = plt.axes(projection='3d')

surf = ax.plot_surface(x, y, z, linewidth=1, rstride=2, cstride=2, cmap=plt.cm.RdPu, antialiased=True )

ax.zaxis.set_major_locator(LinearLocator(6))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_title('Local Volatility Surface')        # title

vals = ax.get_zticks()
ax.set_zticklabels(['   {:3.2f}%'.format(x*100) for x in vals])

ax.view_init(elev=30,azim=150) 
fig.colorbar(surf, shrink=0.5, aspect=5)     # colour bar


#ax.plot_wireframe(x, y, z) #, cmap=plt.cm.jet, rstride=1, cstride=1, linewidth=0)
plt.show()
