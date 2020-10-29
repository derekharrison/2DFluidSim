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

grid_data_file = 'grid_size_p.txt'
nx, ny, nt = np.genfromtxt(grid_data_file, unpack=True)

nx = int(nx)
ny = int(ny)
nt = int(nt)

print(nx)
print(ny)

ux_p_norm = np.zeros((nx, ny))
uy_p_norm = np.zeros((nx, ny))

def animate(i):
    ref_file = 'data_ux_p_'+ str(i) +'.txt'
    ux_p_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_uy_p_'+ str(i) +'.txt'
    uy_p_ref = np.genfromtxt(ref_file, unpack=True)

    ux_p_norm = ux_p_ref/sqrt(ux_p_ref**2 + uy_p_ref**2 + 1e-20)
    uy_p_norm = uy_p_ref/sqrt(ux_p_ref**2 + uy_p_ref**2 + 1e-20)
    ref_file = 'data_x_p_' + str(i) + '.txt'
    x_p_ref = np.genfromtxt(ref_file, unpack=True)
    ref_file = 'data_y_p_' + str(i) + '.txt'
    y_p_ref = np.genfromtxt(ref_file, unpack=True)
    print(i)
    fig.clear()

    #ax = fig.gca(projection='3d')
    quiv = plt.quiver(x_p_ref, y_p_ref, ux_p_norm, uy_p_norm)
    #surf = ax.plot_surface(X_ref, Y_ref, UX_ref)
    #cont = plt.contourf(X_ref, Y_ref, UX_ref)
    #ax.set_xlabel('X')
    #ax.set_xlim(0, 1)
    #plt.xlim((0, 1))
    #ax.set_ylabel('Y')
    #ax.set_ylim(0, 0.1)
    #plt.ylim((0, 0.1))
    #ax.set_zlabel('T')
    #ax.set_zlim(-0.00012, 0.00012)
    #plt.clim((0.0, 0.00012))
    return quiv

size_t = nt-1
anim = manimation.FuncAnimation(fig, animate, frames=size_t, repeat=False)

print("Done Animation, start saving")

anim.save('2DFlowSimQuiv.mp4', writer=writer, dpi=200)
    
print("--- %s seconds ---" % (time.time() - start_time))
