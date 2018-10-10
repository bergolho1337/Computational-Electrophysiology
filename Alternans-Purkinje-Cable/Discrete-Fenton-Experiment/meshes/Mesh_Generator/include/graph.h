#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <map>
#include <queue>

using namespace std;

struct Edge;
struct Node;

// =============================================================================================================
// =============================================================================================================
// Structure of an Edge
struct Edge
{
	int id;				// Identifier of the destination Node
	double w;		    // Size of the Edge, Euclidean distance
	Edge *next;			// Pointer to the next Edge
	Node *dest;			// Pointer to the destination Node
}typedef Edge;

Edge* newEdge (int id, double w, Node *dest);
// =============================================================================================================
// =============================================================================================================
// Structure of a Node
struct Node
{
	int id;					// Identifier of the Node
	double x, y, z;			// Coordinates (x,y,z)
	int num_edges;			// Number of edges of the Node
	Node *next;				// Pointer to the next Node
	Edge *edges;			// Pointer to the list of Edge
}typedef Node;
// =============================================================================================================
Node* newNode (int id, double x, double y, double z);
// =============================================================================================================
// Structure of the Graph
struct Graph
{
	Node *listNodes;			// Pointer to the list of Node
	Node *lastNode;				// Pointer to the last Node of the list
	int total_nodes;			// Total number of Nodes
	int total_edges;			// Total number of Edges
}typedef Graph;

// Graph functions
void initGraph (Graph **g);
Graph* readPurkinjeNetworkFromFile (char *filename);
Node* searchNode (Graph *g, int id);
Node* insertNodeGraph (Graph **g, double p[]);
void insertEdgeGraph (Graph **g, int id_1, int id_2);
map<int,int> BFS (Graph *g, int s);
void printGraph (Graph *g);
// =============================================================================================================
// =============================================================================================================
// Auxiliary functions
double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2);
void error (char *msg);
// =============================================================================================================