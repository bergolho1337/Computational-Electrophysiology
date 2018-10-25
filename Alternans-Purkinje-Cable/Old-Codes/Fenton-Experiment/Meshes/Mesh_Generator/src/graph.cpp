#include "../include/graph.h"

// Initialize the parameters from the graph
void initGraph (Graph **g)
{
    (*g)->total_nodes = 0;
    (*g)->total_edges = 0;
    (*g)->listNodes = NULL;
    (*g)->lastNode = NULL;
}

// Read a .vtk file and convert it to a graph structure
Graph* readPurkinjeNetworkFromFile (char *filename)
{
    FILE *inFile = fopen(filename,"r");
    if (inFile == NULL) error("Cannot open VTK file!");
    Graph *g = (Graph*)malloc(sizeof(Graph));
    initGraph(&g);
    // Read nodes
    int N;
    char str[100];
    while (fscanf(inFile,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    if (!fscanf(inFile,"%d",&N)) error("Reading file");
    if (!fscanf(inFile,"%s",str)) error("Reading file");
    for (int i = 0; i < N; i++)
    {
        double p[3];
        if (!fscanf(inFile,"%lf %lf %lf",&p[0],&p[1],&p[2])) error("Reading file");
        insertNodeGraph(&g,p);
    }
    // Read edges
    int A, trash;
    while (fscanf(inFile,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(inFile,"%d %d",&A,&trash)) error("Reading file");
    for (int i = 0; i < A; i++)
    {
        int e[2];
        if (!fscanf(inFile,"%d %d %d",&trash,&e[0],&e[1])) error("Reading file");
        insertEdgeGraph(&g,e[0],e[1]);
        //insertEdgeGraph(&g,e[1],e[0]);
    }
    fclose(inFile);
    return g;
}

// Node constructor
Node* newNode (int id, double x, double y, double z)
{
    Node *node = (Node*)malloc(sizeof(Node));
    node->id = id;
    node->x = x;
    node->y = y;
    node->z = z;
    node->num_edges = 0;
    node->next = NULL;
    node->edges = NULL;
    return node;
}

// Edge constructor
Edge* newEdge (int id, double w, Node *dest)
{
	Edge *edge = (Edge*)malloc(sizeof(Edge));
	edge->id = id;
	edge->w = w;
	edge->dest = dest;
	edge->next = NULL;
	return edge;
}


// Insert a Node to graph given its position (px,py,pz)
Node* insertNodeGraph (Graph **g, double p[])
{
    Node *ptr = (*g)->listNodes;
    Node *ptrNew = newNode((*g)->total_nodes++,p[0],p[1],p[2]);
    // First node of the list
    if (ptr == NULL)
        (*g)->listNodes = ptrNew;
    // Iterate over the list and insert to the last
    else
    {
        while (ptr->next != NULL)
            ptr = ptr->next;
        ptr->next = ptrNew;
    }
    return ptrNew;
}

// Insert an edge to the graph given the identifiers of both points
void insertEdgeGraph (Graph **g, int id_1, int id_2)
{
	Node *ptr1, *ptr2;
	Edge *edge;
	double norm;
	// Check if the edge is invalid
	if (id_1 == id_2) return;

	ptr1 = searchNode(*g,id_1);
	ptr2 = searchNode(*g,id_2);
	
    norm = calcNorm(ptr1->x,ptr1->y,ptr1->z,ptr2->x,ptr2->y,ptr2->z);
    edge = newEdge(id_2,norm,ptr2);
    // First edge
    if (ptr1->edges == NULL)
        ptr1->edges = edge;
    // Iterate over the list and insert to the last 
    else
    {
        Edge *ptrl = ptr1->edges;
        while (ptrl->next != NULL)
            ptrl = ptrl->next;
        ptrl->next = edge;
    }
    // Increment the number of edges of origin Node
    ptr1->num_edges++;
    // Increment the total number of edges from the graph
    (*g)->total_edges++;
}

// Search for a Node given an identifier
Node* searchNode (Graph *g, int id)
{
	Node *ptr = g->listNodes;
	while (ptr != NULL)
	{
		if (ptr->id == id)
			return ptr;
		ptr = ptr->next;
	}
    printf("[-] ERROR! Node %d was not found!\n",id);
    error("Node not found!");
    return NULL;
}


void error (char *msg)
{
    printf("[-] ERROR! %s\n",msg);
    exit(-1);   
}

// Euclidean norm
double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
}

void printGraph (Graph *g)
{
	Node *ptr;
	Edge *ptrl;
	ptr = g->listNodes;
	printf("======================= PRINTING GRAPH ================================\n");
	while (ptr != NULL)
	{
	    printf("|| %d (%.2lf %.2lf %.2lf) %d ||",ptr->id,ptr->x,ptr->y,ptr->z,ptr->num_edges);
		ptrl = ptr->edges;
		while (ptrl != NULL)
		{
			printf(" --> || %d %.2lf (%.2lf %.2lf %.2lf) ||",ptrl->id,ptrl->w,ptrl->dest->x,ptrl->dest->y, \
					ptrl->dest->z);
			ptrl = ptrl->next;
		}
		printf("\n");
		ptr = ptr->next;
	}
	printf("=======================================================================\n");
}

// Bread-First-Search
map<int,int> BFS (Graph *g, int s)
{
    Node *ptr;
    Edge *ptrl;
    map<int,int> dist;
    queue<int> q;

    // Inicialization variables
    dist[0] = 0;
    q.push(0);

    // Run BFS
    while (!q.empty())
    {
        int u = q.front(); q.pop();
        ptr = searchNode(g,u);
        ptrl = ptr->edges;
        while (ptrl != NULL)
        {
            if (!dist.count(ptrl->id))
            {
                dist[ptrl->id] = dist[ptr->id] + 1;
                q.push(ptrl->id);
            }
            ptrl = ptrl->next;
        }
    }

    return dist;
}