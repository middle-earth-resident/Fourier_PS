from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import numpy as np

L = 80
N = 256
dx = L/N

x = np.arange(-L/2,L/2+dx,dx)
t = np.arange(0,10,1)
u = np.loadtxt("den.txt",usecols=(0),unpack=True)
X,T = np.meshgrid(x,t)
U = np.reshape(u,(10,257))
print(X.shape)
fig = plt.figure()
ax = plt.axes(projection="3d")
ax.plot_surface(T,X,U)
plt.savefig("den.pdf")
