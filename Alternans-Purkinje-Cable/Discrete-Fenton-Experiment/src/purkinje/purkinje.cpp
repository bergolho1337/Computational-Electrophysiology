#include "purkinje.h"

void configure_purkinje_network_from_options (struct graph *the_purkinje_network, struct user_options *config)
{
    assert(the_purkinje_network);
    assert(config);

    char *network_filename = config->pk_network_filename;
    double size_cell = config->start_h;
    int num_div_cell = config->num_div_cell;
    double size_volume = size_cell / num_div_cell;

    fprintf(stdout,"Loading Purkinje Network: %s\n",network_filename);
    fprintf(stdout,"Size cell = %.10lf cm\n",size_cell);
    fprintf(stdout,"Size control volume = %.10lf\n",size_volume);

    set_purkinje_network_from_file(the_purkinje_network,network_filename,size_volume,num_div_cell);
}

void set_purkinje_network_from_file(struct graph *the_purkinje_network, char *network_filename, const double size_volume, const int num_div_cell)
{
    struct graph *skeleton_mesh = new_graph();

    build_skeleton_mesh(network_filename,skeleton_mesh);
    //print_graph(skeleton_mesh);

    build_mesh_purkinje(the_purkinje_network,skeleton_mesh,size_volume,num_div_cell);

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

void build_mesh_purkinje (struct graph *the_purkinje_network, struct graph *the_skeleton, const double size_volume, const int num_div_cell)
{
    assert(the_purkinje_network);
    assert(the_skeleton);

    uint32_t n = the_skeleton->total_nodes;
    // This map is needed to deal with bifurcations
    uint32_t *map_skeleton_to_mesh = (uint32_t*)calloc(n,sizeof(uint32_t));

    // Construct the first node
    struct node *tmp = the_skeleton->list_nodes;
    double pos[3]; pos[0] = tmp->x; pos[1] = tmp->y; pos[2] = tmp->z;
    insert_node_graph(the_purkinje_network,pos);
}