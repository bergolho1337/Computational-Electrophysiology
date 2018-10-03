import os
import sys
import numpy as np
import matplotlib.pyplot as plt

TOLERANCE = 1.0e-01

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

def calc_apd (t,v,start_time,end_time):
    t1 = start_time
    t2 = start_time
    min_diff = 100
    v_base = v[0]
    for i in range(1,len(t)):
        if (t[i] >= start_time and t[i] < end_time):
            diff = abs(v[i]-v_base)
            diff = np.sqrt(diff) 
            if (diff < min_diff):
                min_diff = diff
                t2 = t[i]
    return t2 - t1

def calc_apds (output_dir,dir_files,cells_ids,start_time,end_time):
    apds = []
    for i in range(len(cells_ids)):
        print "Cell = %d" % cells_ids[i]
        T, V, M, H, N = get_state_vector_files(output_dir,dir_files,cells_ids[i])
        apd = calc_apd(T,V,start_time,end_time)
        apds.append(apd)
    return apds

def main():

    if (len(sys.argv) != 8):
        print("Usage:> %s <output_dir> <cell_id_1> <cell_id_2>\
                          <cell_id_3> <cell_id_4> <cell_id_5> <cell_id_6>" % (sys.argv[0]))
        sys.exit(1)
    else:
        output_dir = sys.argv[1]
        cells_ids = []
        for i in range(6):
            cells_ids.append(int(sys.argv[2+i]))

        dir_files = os.listdir(output_dir)
        dir_files.sort()

        apds = calc_apds(output_dir,dir_files,cells_ids,0,300)
        print apds


        #T, V, M, H, N = get_state_vector_files(output_dir,dir_files,cell_id)
        
        #apd = calc_apd(T,V,0.0,300.0)

        #print ("APD = %.10lf" % apd)

        #plot_state_vector(out_filename,T,V,"V")        

        #plt.clf()
        #plot_state_vector(out_filename,T,M,"m")

        #plt.clf()
        #plot_state_vector(out_filename,T,H,"h")

        #plt.clf()
        #plot_state_vector(out_filename,T,N,"n")        

if __name__ == "__main__":
    main()