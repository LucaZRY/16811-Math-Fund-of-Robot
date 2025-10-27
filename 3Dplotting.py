from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.patches import Circle, Ellipse
import matplotlib.pyplot as plt

def add_point(ax, x, y, z, fc = None, ec = None, radius = 0.005):
    # https://stackoverflow.com/a/65115447/5487412
       xy_len, z_len = ax.get_figure().get_size_inches()
       axis_length = [x[1] - x[0] for x in [ax.get_xbound(), ax.get_ybound(), ax.get_zbound()]]
       axis_rotation =  {'z': ((x, y, z), axis_length[1]/axis_length[0]),
                         'y': ((x, z, y), axis_length[2]/axis_length[0]*xy_len/z_len),
                         'x': ((y, z, x), axis_length[2]/axis_length[1]*xy_len/z_len)}
       for a, ((x0, y0, z0), ratio) in axis_rotation.items():
           p = Ellipse((x0, y0), width = radius, height = radius*ratio, fc=fc, ec=ec)
           ax.add_patch(p)
           art3d.pathpatch_2d_to_3d(p, z=z0, zdir=a)

def plot_path(path, figsize=(7,7)):
    tt = path.shape[0]
    path_values = np.zeros((tt, 1))
    for i in range(tt):
        path_values[i] = obs_cost[int(np.floor(path[i, 0])), int(np.floor(path[i, 1]))]

    # Plot 2D
    plt.figure(figsize=figsize)
    plt.imshow(obs_cost.T)
    plt.plot(path[:, 0], path[:, 1], 'ro')

    # Plot 3D
    fig3d = plt.figure(figsize=figsize)
    ax3d = fig3d.add_subplot(111, projection='3d')
    xx, yy = np.meshgrid(range(N), range(N))
    ax3d.plot_surface(xx, yy, obs_cost.T, cmap=plt.get_cmap('coolwarm'))
    ax3d.scatter(path[:, 0], path[:, 1], path_values, s=20, c='r', alpha=1)
    for i,(x,y) in enumerate(path):
        z = path_values[i][0]
        add_point(ax3d, x, y, z, fc="r", radius=1)
    ax3d.view_init(elev=47, azim=27)
    plt.show()
