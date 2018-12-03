#include "purkinje.h"

Graph* set_purkinje_mesh_from_file (const string network_filename, const double h)
{
    Graph *the_purkinje_network = initialize_purkinje_mesh(network_filename,h);
    
    return the_purkinje_network;
}

Graph* initialize_purkinje_mesh (const string filename, const double h)
{
    Graph *mesh_network = new Graph();
    Graph *skeleton_network = new Graph();

    build_skeleton_purkinje(skeleton_network,filename);
    //skeleton_network->print();

    build_mesh_purkinje(mesh_network,skeleton_network,h);
    //mesh_network->print();

    delete skeleton_network;

    return mesh_network;
}

void build_mesh_purkinje (Graph *mesh, Graph *sk, const double h)
{
    int n = sk->get_total_nodes();
    int *map_skeleton_to_mesh = new int[n]();

    // Insert the first node
    Node *tmp = sk->get_list_nodes();
    double pos[3];
    pos[0] = tmp->x; pos[1] = tmp->y; pos[2] = tmp->z;
    mesh->insert_node_graph(pos);

    // Make a DFS to build the mesh of the Purkinje network
    depth_first_search(mesh,tmp,0,map_skeleton_to_mesh,h);

    delete [] map_skeleton_to_mesh;

}

void depth_first_search(Graph *mesh, Node *u, int level, int *map_skeleton_to_mesh, const double h)
{
    Edge *v = u->list_edges;
    while (v != NULL)
    {
        grow_segment(mesh,u,v,map_skeleton_to_mesh,h);
        depth_first_search(mesh,v->dest,level+1,map_skeleton_to_mesh,h);
        v = v->next;
    }
}

void grow_segment (Graph *mesh, Node *u, Edge *v, int *map_skeleton_to_mesh, const double h)
{
    double d_ori[3], d[3];
    double segment_length = v->w;
    uint32_t n_points = segment_length / h;

    // Capture the index of the growing node on the mesh
    uint32_t id_source = map_skeleton_to_mesh[u->id];

    // Calculate a unitary direction vector of the segment
    calc_unitary_vector(d_ori,u,v->dest);

    // Copy the position of the source node
    d[0] = u->x;
    d[1] = u->y;
    d[2] = u->z;

    fprintf(stdout,"[Purkinje] Node %d will grow %d points\n",u->id,n_points);

    // Grow the number of points of size 'h' until reaches the size of the segment
    for (int k = 1; k <= n_points; k++)
    {
        double pos[3];
        pos[0] = d[0] + d_ori[0]*h*k;
        pos[1] = d[1] + d_ori[1]*h*k;
        pos[2] = d[2] + d_ori[2]*h*k;

        mesh->insert_node_graph(pos);
        mesh->insert_edge_graph(id_source,mesh->get_total_nodes()-1);
        mesh->insert_edge_graph(mesh->get_total_nodes()-1,id_source);
        
        id_source = mesh->get_total_nodes()-1;
    }

    // Save the last inserted node index, in case this node generates offsprings
    map_skeleton_to_mesh[v->id] = id_source;
}

void calc_unitary_vector (double d_ori[], Node *u, Node *v)
{
    d_ori[0] = v->x - u->x;
    d_ori[1] = v->y - u->y;
    d_ori[2] = v->z - u->z;
    double norm = sqrt(d_ori[0]*d_ori[0] + d_ori[1]*d_ori[1] + d_ori[2]*d_ori[2]);
    for (int i = 0; i < 3; i++)
        d_ori[i] /= norm;
}

void build_skeleton_purkinje (Graph *sk, const string filename)
{
    FILE *file = fopen(filename.c_str(),"r");
    if (!file)
    {
        fprintf(stderr,"Error opening Purkinje mesh described in %s!!\n", filename);
        exit (EXIT_FAILURE);
    }

    int N;
    char str[100];

    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"POINTS") == 0) break;
    
    if (!fscanf(file,"%d",&N))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    if (!fscanf(file,"%s",str))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }

    // Read points
    for (int i = 0; i < N; i++)
    {
        double pos[3];
        if (!fscanf(file,"%lf %lf %lf",&pos[0],&pos[1],&pos[2]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        sk->insert_node_graph(pos);
    }

    // Read edges
    int trash, E;
    while (fscanf(file,"%s",str) != EOF)
        if (strcmp(str,"LINES") == 0) break;
    if (!fscanf(file,"%d %d",&E,&trash))
    {
        fprintf(stderr,"Error reading file.\n");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < E; i++)
    {
        int e[2];
        if (!fscanf(file,"%d %d %d",&trash,&e[0],&e[1]))
        {
            fprintf(stderr,"Error reading file.\n");
            exit(EXIT_FAILURE);
        }
        sk->insert_edge_graph(e[0],e[1]);
    }

    fclose(file);
}