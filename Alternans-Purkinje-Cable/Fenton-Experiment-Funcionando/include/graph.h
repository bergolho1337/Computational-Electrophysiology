#ifndef GRAPH_H_
#define GRAPH_H_

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <fstream>

using namespace std;

const double INF = __DBL_MAX__;

class Edge;
class Node;

// =============================================================================================================
// =============================================================================================================
// Structure for an Edge of the graph
class Edge
{
public:
    Edge (int id, double w, Node *dest);

public:
	int id;				// Identifier of the destination node
	double w;		    // Size of the edge, euclidean distance
	Edge *next;			// Pointer to the next Edge
	Node *dest;			// Pointer to the destination Node
};

// =============================================================================================================
// =============================================================================================================
// Strucuture for a Node of the graph
class Node
{
public:
    Node (int id, int type, double x, double y, double z);
	Node (int id, int type, double x, double y, double z, double d);

public:
	int type;				// 0 = Purkinje cell || 1 = PMJ
	int id;					// Identifier of the Node 
	double x, y, z;			// Coordinates (x,y,z)
	double d;				// Diameter of the cell
	int num_edges;			// Number of edges
	Node *next;				// Pointer to the next Node
	Edge *edges;			// Pointer to the list of Edges
};
// =============================================================================================================
// =============================================================================================================
// Structure of the Graph
class Graph
{
public:
    Graph (string filename, double &dx);
    void print ();
	void printterm ();
    void error (const char msg[]);
	void dijkstra (int s);
	// Inline
	int getTotalNodes () { return total_nodes; }
	int getTotalEdges () { return total_edges; }
	Node* getListNodes () { return listNodes; }
	double* getDist () { return dist; }
	int* getTerm () { return term; }
	int getNTerm() { return nterm; }
private:
	Node *listNodes;			// Pointer to the lists of Nodes
	Node *lastNode;				// Pointer to the last Node of the list
	int total_nodes;			// Total number of Nodes
	int total_edges;			// Total number of Edges
	int nterm;					// Number of terminals
	double *dist;				// Distance from the source node to all the others
	int *term;					// Pointer to the terminals 

    void initGraph ();
    void insertNodeGraph (int type, double p[]);
	void insertEdgeGraph (int id_1, int id_2);
	Node* searchNode (int id);
	void insertPMJ ();
	void calcPosition (Node *p1, Node *p2, double p[]);
	void setTerm ();

};

// =============================================================================================================
// =============================================================================================================
// Funcoes auxiliares
double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2);
void calcPosition (Node *p1, Node *p2, double p[]);
// =============================================================================================================

#endif