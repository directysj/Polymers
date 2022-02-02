from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import warnings
from numpy import *
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 4, 3

################ Plot multichain linear polymers ##############

npc = 20
nc = 5
#nc = 250

fig = plt.figure()
ax = plt.axes(projection="3d")

#a = np.loadtxt('chains.out')
a = np.loadtxt('chkCrd5k.dat')

npt = nc*npc

for i in range(0, nc, 1):

      end1 = i*npc
      end2 = (i+1)*npc
      """
      k1 = end1
      k2 = end2 - 1

      for j in range(k1, k2, 1):

            dx = a[j, 0] - a[j+1, 0]
            dy = a[j, 1] - a[j+1, 1]
            dz = a[j, 2] - a[j+1, 2]

            ds2 = dx*dx + dy*dy + dz*dz

            ds = np.sqrt(ds2)

            if(ds > 1.30):

                      end2 = j

                      ax.plot3D(a[end1:end2, 0],a[end1:end2, 1],a[end1:end2, 2], marker = 'o', markersize = 8, markerfacecolor = "w")
                      end1 = end2 + 1

      ax.plot3D(a[end1:end2-1, 0],a[end1:end2-1, 1],a[end1:end2-1, 2], marker = 'o', markersize = 8, markerfacecolor = "w")
    """
      ax.plot3D(a[end1:end2-1, 0],a[end1:end2-1, 1],a[end1:end2-1, 2], marker = 'o', markersize = 8, markerfacecolor = "w")


#plt.xlabel('X')
#plt.ylabel('Y')
#plt.zlabel('Z')

#plt.xticks([0, 1, 4, 8, 10, 12, 20])

#ax.set_xlim(0.0, bl)
#ax.set_ylim(0.0, bl)
#ax.set_zlim(0.0, bl)

#plt.savefig('deform.png', dpi=500, bbox_inches='tight')
#plt.savefig('deform.eps', dpi=1200, bbox_inches='tight', format='eps')
plt.show()
