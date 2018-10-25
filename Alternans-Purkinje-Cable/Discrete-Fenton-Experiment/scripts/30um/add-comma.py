import os
import sys
from glob import glob

def add_comma_to_files_in_dir (dir_name):
    files = [f for f in glob(dir_name + '*.dat')]
    for f in files:
        aux = f.split('/')
        aux2 = aux[1].split('.')
        filename = aux2[0]
        extension = aux2[1]
        new_filename = dir_name + filename[0:4] + "-" + filename[4:] + "." + extension
        
        print("[!] Changing name of file '%s' to '%s'" % (dir_name+aux[1],new_filename))
        os.system("mv %s %s" % (dir_name+aux[1],new_filename))
        

def main():

    if (len(sys.argv) != 2):
        print("Usage:> %s <input_directory>" % (sys.argv[0]))
        sys.exit(1)
    else:
        input_dir = sys.argv[1]

        add_comma_to_files_in_dir(input_dir)
    

if __name__ == "__main__":
    main()