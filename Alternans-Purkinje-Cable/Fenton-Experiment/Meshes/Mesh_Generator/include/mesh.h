#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include "../include/graph.h"

using namespace std;

struct Mesh;
struct Point;
struct Element;

// Ratio of reduction of the diameter at each growing iteration
const double DELTA = 0.7;      

// Initial diameter of the Purkinje cells
const double D = 0.0055;                         

// Length of Purkinje cells
const double PIG_DX = 0.0068;	                      // T. Stankovičová, 2003 (cm) -- Pig Purkinje cell
const double DOG_DX = 0.0164;                         // Michael F. Sheets (1983)   -- Dog Purkinje cell
const double ALIEN_DX = 0.01;                         // Test 1 -- Alien Purkinje cell
const double ORC_DX = 0.02;				              // Test 2 -- Orc Purinje cell

const double TEST1_DX = 0.0075;
const double TEST2_DX = 0.0100;
const double TEST3_DX = 0.0125;
const double TEST4_DX = 0.0150;
const double TEST5_DX = 0.0175;
const double TEST6_DX = 0.0200;

struct Point
{
    double x, y, z;                 // Coordinates
    double d;                       // Diameter
}typedef Point;

struct Element
{
    int left;                       // Identifier of the left end Node
    int right;                      // Identifier of the right end Node
}typedef Element;

struct Mesh
{
    int nElem;                      // Number of elements
    int nPoints;                    // Number of points
    double h;                       // Distance between to points
    int *map_graph_elem;            // Point mapping (graph -> element)
    vector<Point> points;           // Vector of points
    vector<Element> elements;       // Vector of elements
}typedef Mesh;

Mesh* newMesh (int argc, char *argv[]);
double setTypeCell (char cName[]);
void GraphToMesh (Mesh *mesh, Graph *g);
void calcDirection (double d_ori[], Node *n1, Node *n2);
void writeMeshToFile (Mesh *mesh, char *filename);
void writeMeshToVTK (Mesh *mesh, const char *filename);
void writeLevelToFile (Mesh *mesh, Graph *g);
void writeMeshInfo (Mesh *mesh);
void printMeshInfo (Mesh *mesh);
void changeExtension (char *filename);
void DFS (Mesh *mesh, Node *u, int lvl);
void growSegment (Mesh *mesh, Node *u, Edge *v, double diam);
