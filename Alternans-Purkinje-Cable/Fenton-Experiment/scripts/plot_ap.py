import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def get_state_vector_files (output_dir, dir_files, cell_id):
    t = []
    v = []
    m = []
    h = []
    n = []
    for filename in dir_files:
        data = np.genfromtxt(output_dir+filename)
        t.append(data[cell_id][0])
        v.append(data[cell_id][1])
        m.append(data[cell_id][2])
        h.append(data[cell_id][3])
        n.append(data[cell_id][4])
    T = np.asarray(t)
    V = np.asarray(v)
    M = np.asarray(m)
    H = np.asarray(h)
    N = np.asarray(n)

    return T, V, M, H, N

def plot_state_vector (out_filename,t,sv,sv_name):
    plt.grid()
    plt.plot(t,sv,label="aprox",c="black")
    plt.xlabel(u"t",fontsize=15)
    plt.ylabel(sv_name,fontsize=15)
    plt.title(out_filename,fontsize=14)
    plt.legend(loc=0,fontsize=14)
    plt.show()

def main():

    if (len(sys.argv) != 3):
        print("Usage:> %s <output_dir> <cell_id>" % (sys.argv[0]))
        sys.exit(1)
    else:
        output_dir = sys.argv[1]
        cell_id = int(sys.argv[2])
        out_filename = "Cell " + str(cell_id)

        dir_files = os.listdir(output_dir)
        dir_files.sort()

        T, V, M, H, N = get_state_vector_files(output_dir,dir_files,cell_id)

        plot_state_vector(out_filename,T,V,"V")        

        #plt.clf()
        #plot_state_vector(out_filename,T,M,"m")

        #plt.clf()
        #plot_state_vector(out_filename,T,H,"h")

        #plt.clf()
        #plot_state_vector(out_filename,T,N,"n")        

if __name__ == "__main__":
    main()