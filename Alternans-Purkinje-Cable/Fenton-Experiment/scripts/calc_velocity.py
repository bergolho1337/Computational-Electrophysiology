import os
import sys
import numpy as np
import matplotlib.pyplot as plt

def calc_maximum_derivative (t, v, start_time, end_time):
    max_t = 0
    max_dvdt = 0
    for i in range(len(t)-1):
        if (t[i] >= start_time and t[i] < end_time):
            diff = v[i+1] - v[i]
            if (diff > max_dvdt):
                max_dvdt = diff
                max_t = t[i]
    return max_t, max_dvdt

def get_start_discretization (output_mesh):
    file = open(output_mesh,"r")
    line = file.readline()
    tokens = line.split()
    
    return float(tokens[2])

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

    if (len(sys.argv) != 5):
        print("Usage:> %s <output_dir> <output_mesh> <cell_id_1> <cell_id_2>" % (sys.argv[0]))
        sys.exit(1)
    else:
        output_dir = sys.argv[1]
        output_mesh = sys.argv[2]
        cell_id_1 = int(sys.argv[3])
        cell_id_2 = int(sys.argv[4])
        out_filename = "Velocity between " + str(cell_id_1) + " and " + str(cell_id_2)

        dir_files = os.listdir(output_dir)
        dir_files.sort()

        T1, V1, M1, H1, N1 = get_state_vector_files(output_dir,dir_files,cell_id_1)
        T2, V2, M2, H2, N2 = get_state_vector_files(output_dir,dir_files,cell_id_2)
        
        t1, dvdt1 = calc_maximum_derivative(T1,V1,0.0,300.0)
        t2, dvdt2 = calc_maximum_derivative(T2,V2,0.0,300.0)
        t = t2 - t1

        print ("t1 = %.10lf" % t1)
        print ("dvdt1 = %.10lf" % dvdt1)
        print ("t2 = %.10lf" % t2)
        print ("dvdt2 = %.10lf" % dvdt2)

        h = get_start_discretization(output_mesh)  
        d = h*(cell_id_2 - cell_id_1)

        v = d / t * 1.0e+03
        print ("d = %.10lf" % d)
        print ("t = %.10lf" % t)
        print ("v = %.10lf" % v)

        #plot_state_vector(out_filename,T,V,"V")        

        #plt.clf()
        #plot_state_vector(out_filename,T,M,"m")

        #plt.clf()
        #plot_state_vector(out_filename,T,H,"h")

        #plt.clf()
        #plot_state_vector(out_filename,T,N,"n")        

if __name__ == "__main__":
    main()