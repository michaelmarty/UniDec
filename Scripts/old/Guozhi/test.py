import numpy as np
import matplotlib.pyplot as plt


fig, ax = plt.subplots()

size = 0.5
vals = np.array([[60., 32.], [37., 40.], [29., 10.]])

cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(3)*4)
inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,
       wedgeprops=dict(width=size, edgecolor='w'))

ax.pie(vals.flatten(), radius=1-size, colors=inner_colors,
       wedgeprops=dict(width=size, edgecolor='w'))

ax.set(aspect="equal", title='Pie plot with `ax.pie`')
plt.show()

"""
fig, ax = plt.subplots(subplot_kw=dict(polar=True))
size = 0.3
vals = np.array([[60., 32.], [37., 40.], [29., 10.]])
#normalize vals to 2 pi
valsnorm = vals/np.sum(vals)*2*np.pi
#obtain the ordinates of the bar edges
valsleft = np.cumsum(np.append(0, valsnorm.flatten()[:-1])).reshape(vals.shape)

cmap = plt.get_cmap("tab20c")
outer_colors = cmap(np.arange(3)*4)
inner_colors = cmap(np.array([1, 2, 5, 6, 9, 10]))

ax.bar(x=valsleft[:, 0],
       width=valsnorm.sum(axis=1), bottom=1-size, height=size, xerr=0.1,
       color=outer_colors, edgecolor='w', linewidth=1, align="edge")

ax.bar(x=valsleft.flatten(),
       width=valsnorm.flatten(), bottom=1-2*size, height=size, xerr=0.1,
       color=inner_colors, edgecolor='w', linewidth=1, align="edge")

ax.set(title="Pie plot with `ax.bar` and polar coordinates")
ax.set_axis_off()
plt.show()"""

#xdata = np.array([0.273571103, 0.427924, 0.521518472, 0.606904818, 0.721550878, 0.818676778, 0.169139,  # batch a, 20
                  0.264946629, 0.377567405, 0.499840175, 0.611320635, 0.691793579, 0.791046772, 0.135019464,  # batch b
                  0.257864562, 0.402562646, 0.545817, 0.620564076, 0.662262, 0.789596279, 0.18042874])  # batch c, 40, 80
#ydata = np.array([0.014557637, 0.039132, 0.076345411, 0.093313905, 0.106628672, 0.076113749, -0.02814,  # batch a, 20
                  0.009392632, 0.075266738, 0.115078044, 0.134266394, 0.104007865, 0.128838869, -0.058457228,  # batch b
                  0.034923576, 0.048626812, 0.114883, 0.138451348, 0.062596, 0.09148421, -0.024295153])  # batch c, 40, 80