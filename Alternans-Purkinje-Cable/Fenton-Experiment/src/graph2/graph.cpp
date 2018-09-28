#include "graph.h"

Graph::Graph (string filename, double &dx)
{
    ifstream in(filename.c_str());
    if (!in) error("Cannot open MSH file");

    int V, E;
    initGraph();
    in >> E >> V >> dx;
    // Read vertices
    for (int i = 0; i < V; i++)
    {
        double p[4];
        #ifdef DIAMETER
        in >> p[0] >> p[1] >> p[2] >> p[3];
        #else
        in >> p[0] >> p[1] >> p[2];
        #endif
        insertNodeGraph(0,p);
    }
    // Read edges
    for (int i = 0; i < E; i++)
    {
        int e[2];
        in >> e[0] >> e[1];
        insertEdgeGraph(e[0],e[1]);
        insertEdgeGraph(e[1],e[0]);
    }
    in.close();
    printterm();
    setTerm();
    
    #ifdef PMJ
    insertPMJ();
    #endif

    // Alocate memory for the distance vector
    dist = (double*)malloc(sizeof(double)*total_nodes);
}

Node::Node (int id, int type, double x, double y, double z)
{
    this->type = type;
    this->id = id;
    this->x = x;
    this->y = y;
    this->z = z;
    this->num_edges = 0;
    this->next = NULL;
    this->edges = NULL;
}

Node::Node (int id, int type, double x, double y, double z, double d)
{
    this->type = type;
    this->id = id;
    this->x = x;
    this->y = y;
    this->z = z;
    this->d = d;
    this->num_edges = 0;
    this->next = NULL;
    this->edges = NULL;
}

Edge::Edge (int id, double w, Node *dest)
{
	this->id = id;
	this->w = w;
	this->dest = dest;
	this->next = NULL;
}

void Graph::initGraph ()
{
    total_nodes = 0;
    total_edges = 0;
    listNodes = NULL;
    lastNode = NULL;
}

void Graph::insertNodeGraph (int type, double p[])
{
    Node *ptr = listNodes;
    #ifdef DIAMETER
    Node *ptrNew = new Node(total_nodes++,type,p[0],p[1],p[2],p[3]);
    #else
    Node *ptrNew = new Node(total_nodes++,type,p[0],p[1],p[2]);
    #endif
    // First node
    if (ptr == NULL)
        listNodes = ptrNew;
    // Transverse the list of nodes
    else
    {
        while (ptr->next != NULL) ptr = ptr->next;
        ptr->next = ptrNew;
    }
}

void Graph::insertEdgeGraph (int id_1, int id_2)
{
	Node *ptr1, *ptr2;
	// Check for invalid edge
	if (id_1 == id_2) return;

	ptr1 = searchNode(id_1);
	ptr2 = searchNode(id_2);
	
    double norm = calcNorm(ptr1->x,ptr1->y,ptr1->z,ptr2->x,ptr2->y,ptr2->z);
    Edge *edge = new Edge(id_2,norm,ptr2);
    // First edge
    if (ptr1->edges == NULL)
        ptr1->edges = edge;
    // Transverse to the final
    else
    {
        Edge *ptrl = ptr1->edges;
        while (ptrl->next != NULL) ptrl = ptrl->next;
        ptrl->next = edge;
    }
    // Increment the number of edges
    ptr1->num_edges++;
    // Increment the total number of edges
    total_edges++;
}

Node* Graph::searchNode (int id)
{
	Node *ptr = listNodes;
	while (ptr != NULL)
	{
		if (ptr->id == id) 
            return ptr;
		ptr = ptr->next;
	}
    error("Node not found!");
    return NULL;
}

void Graph::calcPosition (Node *p1, Node *p2, double p[])
{
    double size = 2.00;
    double norm = calcNorm(p1->x,p1->y,p1->z,p2->x,p2->y,p2->z);
    p[0] = (p1->x - p2->x)/norm; p[1] = (p1->y - p2->y)/norm; p[2] = (p1->z - p2->z)/norm;
    p[0] = p1->x + size*p[0]; p[1] = p1->y + size*p[1]; p[2] = p1->z + size*p[2];

    #ifdef DIAMETER
    p[3] = p2->d;
    #endif
}

void Graph::insertPMJ ()
{
    Node *ptr = listNodes;
    while (ptr != NULL)
    {
        double p[4];
        // Node is a leaf, not the root and of type Purkinje cell
        if (ptr->type == 0 && ptr->num_edges == 1 && ptr->id != 0)
        {
            calcPosition(ptr,ptr->edges->dest,p);
            insertNodeGraph(1,p);
            insertEdgeGraph(ptr->id,total_nodes-1);
            insertEdgeGraph(total_nodes-1,ptr->id);
        }
        ptr = ptr->next;
    }
}

void Graph::dijkstra (int s)
{
    printf("[!] Running Dijkstra ... ");
    fflush(stdout);

    for (int i = 0; i < total_nodes; i++) dist[i] = INF;
    dist[s] = 0;
    priority_queue< pair<double,int>, vector< pair<double,int> >, greater< pair<double,int> > > pq;
    pq.push(make_pair(0,s));

    while (!pq.empty())
    {
        pair<double,int> front = pq.top(); pq.pop();
        double d = front.first;
        int u = front.second;
        if (d > dist[u]) continue;
        Edge *ptrl = searchNode(u)->edges;
        while (ptrl != NULL)
        {
            int id = ptrl->id;
            double w = ptrl->w; 
            if (dist[u] + w < dist[id])
            {
                dist[id] = dist[u] + w;
                pq.push(make_pair(dist[id],id));
            }
            ptrl = ptrl->next;
        }
    }
    printf("ok\n");

}

void Graph::setTerm ()
{
    Node *ptr = listNodes;
    nterm = 0;
    while (ptr != NULL)
    {
        if (ptr->num_edges == 1 && ptr->id != 0) nterm++;
        ptr = ptr->next;
    }
    int j = 0;
    term = (int*)malloc(sizeof(int)*nterm);
    ptr = listNodes;
    while (ptr != NULL)
    {
        if (ptr->num_edges == 1 && ptr->id != 0) term[j] = ptr->id, j++;
        ptr = ptr->next;
    }
}

void Graph::print ()
{
    Edge *ptrl;
    Node *ptr = listNodes;
	printf("======================= PRINTING GRAPH ================================\n");
	while (ptr != NULL)
	{
        #ifdef DIAMETER
	    printf("|| %d (%d) (%.4lf %.4lf %.4lf) [%.4lf] %d ||",ptr->id,ptr->type,ptr->x,ptr->y,ptr->z,ptr->d,ptr->num_edges);
		#else
        printf("|| %d (%d) (%.4lf %.4lf %.4lf) %d ||",ptr->id,ptr->type,ptr->x,ptr->y,ptr->z,ptr->num_edges);
		#endif

        ptrl = ptr->edges;
		while (ptrl != NULL)
		{
			printf(" --> || %d %.4lf (%.4lf %.4lf %.4lf) ||",ptrl->id,ptrl->w,ptrl->dest->x,ptrl->dest->y, \
					ptrl->dest->z);
			ptrl = ptrl->next;
		}
		printf("\n");
		ptr = ptr->next;
	}
	printf("=======================================================================\n");
    printf("Number of nodes = %d\n",total_nodes);
    printf("Number of edges = %d\n",total_edges);
}

void Graph::error (const char msg[])
{
    printf("[-] ERROR on Graph ! %s !\n",msg);
    exit(EXIT_FAILURE);
}

double calcNorm (double x1, double y1, double z1, double x2, double y2, double z2)
{
    return sqrt(pow((x1-x2),2) + pow((y1-y2),2) + pow((z1-z2),2));
}

void Graph::printterm ()
{
    Node *ptr = listNodes;
    while (ptr != NULL)
    {
        if (ptr->num_edges == 1)
            printf("Terminal node %d\n",ptr->id);
        else if (ptr->num_edges == 3)
            printf("Bifurcation node %d\n",ptr->id);
        ptr = ptr->next;
    }
}