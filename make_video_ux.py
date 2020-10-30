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

max_and_min_data_file = 'max_and_min_ux.txt'
max_ux, min_ux = np.genfromtxt(max_and_min_data_file, unpack=True)

grid_data_file = 'grid_size.txt'
nx, ny, nt = np.genfromtxt(grid_data_file, unpack=True)

nx = int(nx)
ny = int(ny)
nt = int(nt)

def animate(i):
    ref_file = 'data_ux_'+ str(i) +'.txt'
    UX_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_x_' + str(i) +'.txt'
    X_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_y_' + str(i) + '.txt'
    Y_ref = np.genfromtxt(ref_file, unpack=True)
    print(i)
    fig.clear()

    #ax = fig.gca(projection='3d')
    #surf = ax.plot_surface(X_ref, Y_ref, UX_ref)
    cont = plt.contourf(X_ref, Y_ref, UX_ref, N=256)
    #ax.set_xlabel('X')
    #ax.set_xlim(0, 1)
    plt.xlim((0, 1))
    #ax.set_ylabel('Y')
    #ax.set_ylim(0, 0.1)
    plt.ylim((0, 0.1))
    #ax.set_zlabel('T')
    #ax.set_zlim(min_ux*1.01, max_ux*1.01)
    plt.clim((min_ux*1.01, max_ux*1.01))
    return cont

size_t = nt-1
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('2DFlowSimUx.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
