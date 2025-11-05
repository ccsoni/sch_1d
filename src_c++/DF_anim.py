import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import glob
import os
import re
import sys

column = 1

def read_header(filename):
    with open(filename, "r") as f:
        nmesh_x, nmesh_v = map(int,f.readline().split())
        sigma_x, sigma_v = map(float, f.readline().split())
        hbar = float(f.readline().strip())
        tnow = float(f.readline().strip())

    return nmesh_x, nmesh_v, sigma_x, sigma_v, hbar, tnow

def extract_index(filename):
    match = re.search(r'(\d+)(?=\.\w+$)', filename)
    return int(match.group(1)) if match else -1

def update(frame):
    nmesh_x, nmesh_v, hbar, tnow = read_header(file_list[frame])
    data = np.loadtxt(file_list[frame],skiprows=4)
    Z = data[:, 2].reshape(len(xmesh), len(ymesh))
    mesh.set_array(Z.ravel())  # pcolormesh ÌFXV

    ax.set_title(f"$(N_x, N_v)$ = {nmesh_x, nmesh_v} / $\hbar$ = {hbar:.2e} / $t$ = {tnow:.2f}", fontsize=10)

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
    nmesh_x, nmesh_v, hbar, tnow = read_header(file_list[0])
    data0 = np.loadtxt(file_list[0], skiprows=4)
    x = data0[:,0]
    y = data0[:,1]
    z = data0[:,2]

    xmesh = np.unique(x)
    ymesh = np.unique(y)

    X, Y = np.meshgrid(ymesh, xmesh)
    Z = z.reshape(len(xmesh), len(ymesh))

    all_min, all_max = np.inf, -np.inf
    for f in file_list:
        d = np.loadtxt(f,skiprows=4)[:, 2]
        all_min = min(all_min, np.min(d))
        all_max = max(all_max, np.max(d))


    ax.set_xlabel("x",fontsize=12);
    ax.set_ylabel("v",fontsize=12);
    mesh = ax.pcolormesh(Y, X, Z, cmap='viridis', shading='auto',
                         vmin=all_min, vmax=all_max)

    fig.colorbar(mesh,ax=ax)

    ax.set_title(f"$(N_x, N_v)$ = {nmesh_x, nmesh_v}/ $\hbar$ = {hbar:.2e}/ $t$ = {tnow:.2f}", fontsize=10)

    ani = animation.FuncAnimation(
        fig, update, frames=len(file_list), blit=False, interval=500
    )

    #plt.show()
    plt.close(fig)

    for i, filename in enumerate(file_list):
        nmesh_x, nmesh_v, hbar, tnow = read_header(filename)
        data = np.loadtxt(file_list[i],skiprows=4)
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

        ax.set_xlabel("x")
        ax.set_ylabel("v")

        mesh = ax.pcolormesh(Y, X, Z, cmap='viridis', shading='auto',
                             vmin=all_min, vmax=all_max)

        # coordinate of center of the cross bar (sigma_x, sigma_v)
        if sigma_x < 0.05*(xmesh[int(len(xmesh))-1]-xmesh[0]):
            fc = 0.9
        else:
            fc = 1.0 - 2.0*sigma_x/(xmesh[int(len(xmesh))]-xmesh[0])

        xc = xmesh[int(len(xmesh)*0.9)]
        vc = ymesh[int(len(ymesh)*0.9)]
        xc_lo = xc - 0.5*sigma_x
        xc_hi = xc + 0.5*sigma_x
        vc_lo = vc - 0.5*sigma_v
        vc_hi = vc + 0.5*sigma_v
        # draw the cross bar presenting (sigma_x, sigma_v)
        ax.plot([xc_lo, xc_hi], [vc, vc], color='white', lw=1.5)
        ax.plot([xc, xc], [vc_lo, vc_hi], color='white', lw=1.5)

        fig.colorbar(mesh,ax=ax)

        ax.set_title(f"$(N_x, N_v)$ = {nmesh_x, nmesh_v} / $\hbar$ = {hbar:.2e} / $t$ = {tnow:.2f}", fontsize=10)

        output_base, _= os.path.splitext(filename)
        output_filename = output_base+".png"

        plt.savefig(output_filename, bbox_inches='tight')
        plt.close(fig)
