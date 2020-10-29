import numpy as np
import matplotlib.pyplot as plt
from numpy import *
import matplotlib.animation as manimation
import time
import sys
from mpl_toolkits.mplot3d import Axes3D

start_time = time.time()

plt.rcParams['animation.ffmpeg_path']='C:/Users/d-w-h/Downloads/ffmpeg-20200818-1c7e55d-win64-static/ffmpeg-20200818-1c7e55d-win64-static/bin/ffmpeg.exe'
writer=manimation.FFMpegWriter(bitrate=20000, fps=15)

fig = plt.figure(figsize=(8,8))

grid_data_file = 'grid_size_uy.txt'
nx, ny, nt = np.genfromtxt(grid_data_file, unpack=True)

nx = int(nx)
ny = int(ny)
nt = int(nt)

def animate(i):
    ref_file = 'data_uy_'+ str(i) +'.txt'
    UX_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_x_uy_' + str(i) +'.txt'
    X_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_y_uy_' + str(i) + '.txt'
    Y_ref = np.genfromtxt(ref_file, unpack=True)
    print(i)
    fig.clear()

    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X_ref, Y_ref, UX_ref)
    #cont = plt.contourf(X_ref, Y_ref, UX_ref)
    ax.set_xlabel('X')
    ax.set_xlim(0, 1)
    #plt.xlim((0, 1))
    ax.set_ylabel('Y')
    ax.set_ylim(0, 0.1)
    #plt.ylim((0, 0.1))
    ax.set_zlabel('T')
    ax.set_zlim(-0.05, 0.05)
    #plt.clim((0.0, 0.00012))
    return ax

size_t = nt-1
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('2DFlowSimUy.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
