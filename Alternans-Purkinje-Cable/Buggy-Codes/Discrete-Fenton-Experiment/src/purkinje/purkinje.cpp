#include "purkinje.h"

void configure_purkinje_network_from_options (struct graph *the_purkinje_network, struct user_options *config)
{
    assert(the_purkinje_network);
    assert(config);

    char *network_filename = config->pk_network_filename;
    double size_cell = config->start_h;
    int num_div_cell = config->num_div_cell;
    double size_volume = size_cell / num_div_cell;
    double diameter = config->diameter;

    fprintf(stdout,"*********************************************************************************\n");
    fprintf(stdout,"[Purkinje] Loading Purkinje Network: %s\n",network_filename);
    fprintf(stdout,"[Purkinje] Size cell = %.10lf cm\n",size_cell);
    fprintf(stdout,"[Purkinje] Size control volume = %.10lf cm\n",size_volume);

    set_purkinje_network_from_file(the_purkinje_network,network_filename,size_volume,diameter,num_div_cell);

    set_gap_junctions(the_purkinje_network);

    fprintf(stdout,"*********************************************************************************\n");

    //print_graph(the_purkinje_network);
}

void set_purkinje_network_from_file(struct graph *the_purkinje_network, char *network_filename, const double size_volume, const double diameter, const int num_div_cell)
{
    struct graph *skeleton_mesh = new_graph();

    build_skeleton_mesh(network_filename,skeleton_mesh);

    build_mesh_purkinje(the_purkinje_network,skeleton_mesh,size_volume,diameter,num_div_cell);

    free_graph(skeleton_mesh);
}

void build_skeleton_mesh (char *filename, struct graph *sk)
{
    assert(sk);

    FILE *file = fopen(filename,"r");
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
        insert_node_graph(sk,pos);
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
        insert_edge_graph(sk,e[0],e[1]);
    }

    fclose(file);
}

void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *the_skeleton, const double size_volume, const double diameter, const int num_div_cell)
{
    assert(the_purkinje_network);
    assert(the_skeleton);

    uint32_t n = the_skeleton->total_nodes;
    // This map is needed to deal with bifurcations
    uint32_t *map_skeleton_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));

    // Initialize the size of a control volume and the number of divisions of a cell
    the_purkinje_network->dx = size_volume;
    the_purkinje_network->diameter = diameter;
    the_purkinje_network->num_div_cell = num_div_cell;

    // Construct the first node
    struct node *tmp = the_skeleton->list_nodes;
    double pos[3]; pos[0] = tmp->x; pos[1] = tmp->y; pos[2] = tmp->z;
    insert_node_graph(the_purkinje_network,pos);

    // Make a Depth-First-Search to build the mesh of the Purkinje network
    depth_first_search(the_purkinje_network,tmp,0,map_skeleton_to_mesh);

    free(map_skeleton_to_mesh);
}

void depth_first_search (struct graph *the_purkinje_network, struct node *u, int level, uint32_t *map_skeleton_to_mesh)
{

    struct edge *v = u->list_edges;
    while (v != NULL)
    {
        grow_segment(the_purkinje_network,u,v,map_skeleton_to_mesh);
        depth_first_search(the_purkinje_network,v->dest,level+1,map_skeleton_to_mesh);
        v = v->next;
    }

}

void grow_segment (struct graph *the_purkinje_network, struct node *u, struct edge *v, uint32_t *map_skeleton_to_mesh)
{
    double h = the_purkinje_network->dx;
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
    for (uint32_t k = 1; k <= n_points; k++)
    {
        double pos[3];
        pos[0] = d[0] + d_ori[0]*h*k;
        pos[1] = d[1] + d_ori[1]*h*k;
        pos[2] = d[2] + d_ori[2]*h*k;

        insert_node_graph(the_purkinje_network,pos);
        insert_edge_graph(the_purkinje_network,id_source,the_purkinje_network->total_nodes-1);
        insert_edge_graph(the_purkinje_network,the_purkinje_network->total_nodes-1,id_source);
        
        id_source = the_purkinje_network->total_nodes-1;
    }

    // Save the last inserted node index, in case this node generates offsprings
    map_skeleton_to_mesh[v->id] = id_source;
}

void calc_unitary_vector (double d_ori[], struct node *u, struct node *v)
{
    d_ori[0] = v->x - u->x;
    d_ori[1] = v->y - u->y;
    d_ori[2] = v->z - u->z;
    double norm = sqrt(d_ori[0]*d_ori[0] + d_ori[1]*d_ori[1] + d_ori[2]*d_ori[2]);
    for (int i = 0; i < 3; i++)
        d_ori[i] /= norm;
}

// TO DO: Find a way to apply this function over a bifurcation case ...
void set_gap_junctions (struct graph *the_purkinje_network)
{
    uint32_t current_num_volume = 0;
    uint32_t num_div_cell = the_purkinje_network->num_div_cell;

    // The first node will always be linked in a citoplasmatic manner
    the_purkinje_network->list_nodes->list_edges->link_type = 0;
    current_num_volume++;

    // We start from the second node
    struct node *ptr = the_purkinje_network->list_nodes->next;

    while (ptr != NULL)
    {
        struct edge *ptrl = ptr->list_edges;
        uint32_t u = ptr->id;
        while (ptrl != NULL)
        {
            uint32_t v = ptrl->id;
            // Check if the edge is leaving the node and if 
            // the current volume id is a multiple of the number of divisions
            if (v > u && current_num_volume % num_div_cell == 0)
            {
                // Then the edge is a gap junction link
                ptrl->link_type = 1;
                // Reset the counter to work for the next cell ...
                current_num_volume = 0;
                //Remebering to do the same for the back edge ...
                ptr->next->list_edges->link_type = 1;
            }
            // Else the edge is a normal citoplasmatic link
            else
            {
                // Check if the edge isn't already changed
                if (ptrl->link_type != 1)
                    ptrl->link_type = 0;
            }
            ptrl = ptrl->next;
        }
        current_num_volume++;
        ptr = ptr->next;
    }

}