#include "../include/graph.h"

Graph::Graph ()
{
    last_node = NULL;
    list_nodes = NULL;
    total_nodes = 0;
    total_edges = 0;
}

void Graph::free_list_edges (Node *node)
{
    Edge *e1 = node->list_edges;
    Edge *e2 = node->list_edges->next;

    while (e1 != NULL)
    {
        e1->next = NULL;
        delete e1;
        e1 = e2;
        if (e2 != NULL)
            e2 = e2->next;
    }

    node->list_edges = NULL;
}

void Graph::free_list_nodes ()
{
    Node *n1 = list_nodes;
    Node *n2 = list_nodes->next;

    while (n1 != NULL)
    {
        if (n1->list_edges)
            free_list_edges(n1);
        
        n1->next = NULL;
        delete n1;
        n1 = n2;
        if (n2 != NULL)
            n2 = n2->next;

    }
}

Graph::~Graph ()
{
    if (list_nodes)
        free_list_nodes();
}

Node::Node (const int i, const double pos[])
{
    type = 0;
    id = i;
    x = pos[0];
    y = pos[1];
    z = pos[2];
    num_edges = 0;
    next = NULL;
    list_edges = NULL;
}

Edge::Edge (const int i, const double weight, Node *destination)
{
    id = i;
    w = weight;
    dest = destination;
    next = NULL;
}

void Graph::insert_node_graph (const double pos[])
{
    Node *tmp = list_nodes;
    Node *node = new Node(total_nodes++,pos);

    // First node of the list
    if (!tmp)
    {
        list_nodes = node;
        last_node = node;
    }
    // Insert after the last node and update the pointer
    else
    {
        last_node->next = node;
        last_node = last_node->next;
    }
}

Node* Graph::search_node (const int id)
{
    Node *tmp = list_nodes;
    while (tmp != NULL)
    {
        if (tmp->id == id)
            return tmp;
        tmp = tmp->next;
    }
    cerr << "[-] ERROR! Node " << id << " not found!" << endl;
    return NULL;
}

void Graph::insert_edge_graph (const int id_1, const int id_2)
{
    Node *n1, *n2;
    Edge *edge;
    double norm;

    // Check if the edge is invalid
    if (id_1 == id_2) return;

    n1 = search_node(id_1);
    n2 = search_node(id_2);

    norm = calc_norm(n1->x,n1->y,n1->z,n2->x,n2->y,n2->z);
    edge = new Edge(id_2,norm,n2);

    // First edge
    if (!n1->list_edges)
        n1->list_edges = edge;
    // Iterate over the list and insert to the last edge
    else
    {
        Edge *tmp = n1->list_edges;
        while (tmp->next != NULL)
            tmp = tmp->next;
        tmp->next = edge;
    }

    // Increment the number of edges of origin Node
    n1->num_edges++;
    // Increment the total number of edges of the graph
    total_edges++;
}

double calc_norm (const double x1, const double y1, const double z1,\
                  const double x2, const double y2, const double z2)
{
    return sqrt(pow(x2-x1,2) + pow(y2-y1,2) + pow(z2-z1,2));
}

void Graph::dijkstra (int s)
{
    printf("[!] Running Dijkstra ... \n");

    if (!dist)
        dist = (double*)malloc(sizeof(double)*total_nodes);

    for (int i = 0; i < total_nodes; i++) 
        dist[i] = INF;
    dist[s] = 0;

    priority_queue< pair<double,int>, vector< pair<double,int> >, greater< pair<double,int> > > pq;
    pq.push(make_pair(0,s));

    while (!pq.empty())
    {
        pair<double,int> front = pq.top(); pq.pop();
        double d = front.first;
        int u = front.second;
        if (d > dist[u]) 
            continue;
        Edge *ptrl = search_node(u)->list_edges;
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

}

void Graph::print ()
{
    Edge *ptrl;
    Node *ptr = list_nodes;
	printf("======================= PRINTING GRAPH ================================\n");
	while (ptr != NULL)
	{
        #ifdef DIAMETER
	    printf("|| %d (%.4lf %.4lf %.4lf) [%.4lf] %d ||",ptr->id,ptr->x,ptr->y,ptr->z,ptr->d,ptr->num_edges);
		#else
        printf("|| %d (%.2lf %.2lf %.2lf) %d ||",ptr->id,ptr->x,ptr->y,ptr->z,ptr->num_edges);
		#endif

        ptrl = ptr->list_edges;
		while (ptrl != NULL)
		{
			printf(" --> || %d %.2lf (%.2lf %.2lf %.2lf) * %d * ||",ptrl->id,ptrl->w,ptrl->link_type,\
                                ptrl->dest->x,ptrl->dest->y,ptrl->dest->z);
			ptrl = ptrl->next;
		}
		printf("\n");
		ptr = ptr->next;
	}
	printf("=======================================================================\n");
    printf("Number of nodes = %d\n",total_nodes);
    printf("Number of edges = %d\n",total_edges);
}

void Graph::set_gap_junctions (const int num_div_cell)
{
    int count = 0;
    Node *ptr = list_nodes;

    ptr->list_edges->link_type = 0;
    ptr = ptr->next;
    count++;

    while (ptr != NULL)
    {
        int u = ptr->id;
        Edge *ptrl = ptr->list_edges;
        while (ptrl != NULL)
        {
            int v = ptrl->id;
            // The volumes will be linked by a gap junction
            if (v > u && count % num_div_cell == 0)
            {
                ptrl->link_type = 1;
                count = 0;
                ptr->next->list_edges->link_type = 1;
            }
            // The volumes will be linker by citoplasm
            else
            {
                if (ptrl->link_type != 1)
                    ptrl->link_type = 0;
            }
            ptrl = ptrl->next;
        }
        count++;
        ptr = ptr->next;
    }
}

/*
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
*/