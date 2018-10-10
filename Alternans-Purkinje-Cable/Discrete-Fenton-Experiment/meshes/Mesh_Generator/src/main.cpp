#include <cstdio>
#include <cstdlib>
#include "../include/mesh.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 < 3)
    {
        printf("============================================================================================\n");
        printf("Generate the mesh where the problem will be solved.\n");
        printf("You need to provide a skeleton of the network.\n");
        printf("And the length of cells that make up the fibers.\n");
        printf("alien = 100 um\n");
        printf("dog = 164 um\n");
        printf("orc = 200 um\n");
        printf("pig = 68 um\n");
        printf("--------------------------------------------------------------------------------------------\n");
        printf("Usage:> %s <in_VTK_file> <out_MSH_file> --<type_cell>\n",argv[0]);
        printf("<in_VTK_file> = Name of the input .vtk file (skeleton)\n");
        printf("<out_MSH_file> = Name of the output .msh file\n");
        printf("<type_cell> = alien/dog/orc/pig\n");
        printf("Try for example:> %s ./Skeleton-Networks/test1.vtk ./test1.msh --alien\n",argv[0]);
        printf("============================================================================================\n");
        exit(EXIT_FAILURE);
    }
    else
    {
        printf("============================================================================================\n");
        Mesh *mesh = newMesh(argc,argv);
        writeMeshToFile(mesh,argv[2]);
        writeMeshToVTK(mesh,"mesh.vtk");
        printf("[+] Done\n");
        printf("============================================================================================\n");
        return 0;
    }
}