#!/usr/bin/env python3
#
# Newton-Rhapson minimizer of the Rosenbrock function

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation


from Rosenbrock import rosenbrock

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

parameter_file = open("Rosenbrock_NR.log")
x = []
y = []
z = []
for line in parameter_file.readlines():
    parameter = line.split(', ')
    x.append(float(parameter[0]))
    y.append(float(parameter[1]))
    z.append(float(parameter[2]))
data = [np.array((x,y,z))]
print(len(x))
#print(data[0])
#fig = plt.figure()
#ax = fig.gca(projection='3d')
# Make data.
s = 0.05  # Try s=1, 0.25, 0.1, or 0.05
xx = np.arange(-2, 2. + s, s)  # Could use linspace instead if dividing
yy = np.arange(-2, 3. + s, s)  # evenly instead of stepping...

# Create the mesh grid(s) for all X/Y combos.
xx, yy = np.meshgrid(xx, yy)

zz = rosenbrock(xx, yy)

# Try coolwarm vs jet


# surf = ax.scatter(x, y, z, s=0.001, c=np.arange(len(x)))
# line, = ax.scatter(x, y, z, s=1, c=np.arange(len(x)))

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)
ax.elev = 24
ax.azim = 26
surf = ax.plot_surface(xx, yy, zz, norm=LogNorm(), rstride = 1, cstride = 1, cmap = cm.coolwarm,
                       linewidth=0, antialiased=False)
# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
#lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
#ax.set_xlim3d([min(x), max(x)])
#ax.set_xlabel('X')

#ax.set_ylim3d([min(y), max(y)])
#ax.set_ylabel('Y')

#ax.set_zlim3d([min(z), max(z)])
#ax.set_zlabel('Z')

#ax.set_title('3D Test')



# Creating the Animation object
#ani = animation.FuncAnimation(fig, update_lines, 70688, fargs=(data, lines),
#                                   interval=5, blit=False, save_count=70688)

#ani.save('./line.mp4', writer=animation.FFMpegWriter(fps=20, metadata=dict(artist='Me'), bitrate=1800))
plt.show()
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.01f'))

# for ii in range(0,360):
#         ax.elev += 1
#         ax.azim += 1
#         plt.savefig("./Movie/movie%d.png" % ii)

fig.colorbar(surf, shrink=0.5, aspect=5)

# Displays the figure, handles user interface, returns when user closes window
plt.show()
