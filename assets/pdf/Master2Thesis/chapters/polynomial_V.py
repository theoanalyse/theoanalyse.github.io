import numpy as np
import matplotlib.pyplot as plt

mf = 0.87
mb = 0.68
ml = 0.05
me = 0.60
m1 = 5.36
m2 = 9.68
m3 = 17.27

A = -ml*mb - (m3*mb)/me + (m3*mb**2)/(m1*me)
B = -ml*mf + m2*mb - (m2*mb**2)/(m1) + mb**2 - m1*mb - (m3*mf)/(me) + (2*mf*mb*m3)/(me*m1)  
C = m2*mf - (2*mf*mb*m2)/(m1) + mf*mb + (m3*mf**2)/(me*m1)
D = -(m2*mf**2)/(m1)

P = [A, B, C, D]
R = lambda V: m1/(mf + mb*V) - 1/V
S = lambda V: m3/me*(1 - (mf+mb*V)/(m1*V))


roots = np.roots(P)
# roots = list(filter(lambda x : x > 0, roots))

xdata = [R(root) for root in roots]
zdata = [S(root) for root in roots]

for e in zip(xdata, roots, zdata):
	print(np.round(e, 5))


# fig = plt.figure()
# ax = fig.add_subplot(projection="3d")
# ax.scatter3D(xdata, roots, zdata, color="black", marker="*")
# ax.scatter3D([0], [0], [0], color="black", marker="*")

# X = np.linspace(0, np.max(xdata), 5)
# Y = np.linspace(0, np.max(roots), 5)
# Z = np.linspace(0, np.max(zdata), 5)

# xx_z, yy_z = np.meshgrid(X,Y)
# xx_y, zz_y = np.meshgrid(X,Z)
# yy_x, zz_x = np.meshgrid(Y,Z)
# zero = 0*xx_z - 0*yy_z
# ax.plot_surface(xx_z, yy_z, zero, alpha=0.1, color="gray")
# ax.plot_surface(xx_y, zero, zz_y, alpha=0.1, color="gray")
# ax.plot_surface(zero, yy_x, zz_x, alpha=0.1, color="gray")

# plt.show()