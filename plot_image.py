from mpl_toolkits.mplot3d.axes3d import Axes3D
import matplotlib.pyplot as plt
import struct
import numpy as np
import seaborn as sns


# axes3d = Axes3D(fig)

def plot(filename,ax,title=""):
    with open(filename,'rb') as fp:
        (M,),(N,) = struct.iter_unpack('i',fp.read(4*2))
        print(M,N)
        data = np.ndarray((M,N))
        for i in range(M):
            for j in range(N):
                (data[i][j],) = struct.unpack('d',fp.read(8))
        # data[1][1]=0
    # axes3d = Axes3D(fig)
    # axes3d.plot_surface(range(M),range(N),data)

    # ax.imshow(data)
    sns.heatmap(data)
    plt.title(title)
    return data

if __name__ == '__main__':

    filename1 = f"Vr.dat"
    fig = plt.figure()
    # ax = fig.add_subplot(121)
    # plot(filename1,ax,"out1")
    ax = fig.add_subplot(111)
    plot(filename1,ax,"")

    plt.xlabel("a1/RCount")
    plt.ylabel("a2/RCount")
    plt.show()


    # plt.show()
