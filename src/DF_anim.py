import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re
import sys

column = 1

def extract_index(filename):
    match = re.search(r'(\d+)(?=\.\w+$)', filename)
    return int(match.group(1)) if match else -1

def update(frame):
    data = np.loadtxt(file_list[frame])
    Z = data[:, 2].reshape(len(xmesh), len(ymesh))
    mesh.set_array(Z.ravel())  # pcolormesh ÌFXV
    ax.set_title(file_list[frame])
    return mesh,

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        print("Usage: python "+sys.argv[0]+" \"model_??.dat\" column")
        print("file list is given using the wildcard and must be double quated")
        sys.exit(1)

    model_name = sys.argv[1]
    
    file_list = sorted(glob.glob(sys.argv[1]), key=extract_index)

    fig = plt.figure(dpi=200)

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
    ax.grid(axis='x',which='major', color='#e9e9e9')
    ax.grid(axis='y',which='major', color='#e9e9e9')

    # read the first data
    data0 = np.loadtxt(file_list[0])
    x = data0[:,0]
    y = data0[:,1]
    z = data0[:,2]

    xmesh = np.unique(x)
    ymesh = np.unique(y)

    X, Y = np.meshgrid(ymesh, xmesh)
    Z = z.reshape(len(xmesh), len(ymesh))

    all_min, all_max = np.inf, -np.inf
    for f in file_list:
        d = np.loadtxt(f)[:, 2]
        all_min = min(all_min, np.min(d))
        all_max = max(all_max, np.max(d))
    
    
    mesh = ax.pcolormesh(Y, X, Z, cmap='viridis', shading='auto',
                         vmin=all_min, vmax=all_max)
    
    fig.colorbar(mesh,ax=ax)
    ax.set_title(file_list[0])


    ani = animation.FuncAnimation(
        fig, update, frames=len(file_list), blit=False, interval=500
    )

    plt.show()
    plt.close(fig)

    for i, filename in enumerate(file_list):
        data = np.loadtxt(file_list[i])
        x = data[:,0]
        y = data[:,1]
        z = data[:,2]

        xmesh = np.unique(x)
        ymesh = np.unique(y)

        X, Y = np.meshgrid(ymesh, xmesh)
        Z = z.reshape(len(xmesh), len(ymesh))

        fig = plt.figure(dpi=200)

        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'

        ax = fig.add_axes((0.1, 0.1, 0.8, 0.8))
        ax.grid(axis='x',which='major', color='#e9e9e9')
        ax.grid(axis='y',which='major', color='#e9e9e9')

        mesh = ax.pcolormesh(Y, X, Z, cmap='viridis', shading='auto',
                             vmin=all_min, vmax=all_max)
    
        fig.colorbar(mesh,ax=ax)

        output_base, _= os.path.splitext(filename)
        output_filename = output_base+".png"

        plt.savefig(output_filename, bbox_inches='tight')
        plt.close(fig)
