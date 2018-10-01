#include "grid.h"

struct grid* new_grid()
{
    struct grid* result = (struct grid*) malloc(sizeof(struct grid));
    result->first_cell = NULL;
    result->active_cells = NULL;

    // Purkinje
    result->the_purkinje_network = new_graph();

    return result;
}

void initialize_and_construct_grid_purkinje (struct grid *grid)
{
    assert(grid);

    initialize_grid_purkinje(grid);
    construct_grid_purkinje(grid);
}

void initialize_grid_purkinje (struct grid *grid)
{
    assert(grid);

    grid->number_of_cells = 0;
}

void construct_grid_purkinje (struct grid *the_grid)
{
    assert(the_grid);

    double side_length = the_grid->the_purkinje_network->dx;
    double half_side_length = side_length / 2.0f;
    double quarter_side_length = side_length / 4.0f;
    printf("Side length = %lf\n",side_length);
    printf("Half_side length = %lf\n",half_side_length);
    printf("Quarter_side length = %lf\n",quarter_side_length);

    int total_nodes = the_grid->the_purkinje_network->total_nodes;
    
    // Create an array of cell nodes
    struct cell_node **cells = (struct cell_node**)malloc(sizeof(struct cell_node*)*total_nodes);
    for (int i = 0; i < total_nodes; i++)
        cells[i] = new_cell_node();
    
    // Pass through the Purkinje graph and set the cell nodes.
    struct node *n = the_grid->the_purkinje_network->list_nodes;
    for (int i = 0; i < total_nodes; i++)
    {
        
        if (i == 0)
            set_cell_node_data(cells[i],side_length,half_side_length,\
                                NULL,cells[i+1],i,\
                                n->x,n->y,n->z);
        else if (i == total_nodes-1)
            set_cell_node_data(cells[i],side_length,half_side_length,\
                                cells[i-1],NULL,i,\
                                n->x,n->y,n->z);
        else
            set_cell_node_data(cells[i],side_length,half_side_length,\
                                cells[i-1],cells[i+1],i,\
                                n->x,n->y,n->z);

        n = n->next;
    }
    
    // Grid initialization
    the_grid->first_cell = cells[0];
    the_grid->number_of_cells = total_nodes;
    
}

void order_grid_cells (struct grid *the_grid) 
{

    struct cell_node *grid_cell;
    grid_cell = the_grid->first_cell;

    //Here we allocate the maximum number of cells we will need for the whole simulation
    if (the_grid->active_cells == NULL) 
    {
        the_grid->active_cells = (struct cell_node **)malloc (sizeof (struct cell_node *) * the_grid->number_of_cells);
    }

    uint32_t counter = 0;
    while (grid_cell != 0) 
    {
        if (grid_cell->active) 
        {
            grid_cell->grid_position = counter;
            the_grid->active_cells[counter] = grid_cell;
            counter++;
        }

        grid_cell = grid_cell->next;
    }

    the_grid->num_active_cells = counter;
}

void save_old_cell_positions (struct grid *the_grid) 
{

    uint32_t n_active = the_grid->num_active_cells;
    struct cell_node **ac = the_grid->active_cells;

	int i;

	#pragma omp parallel for
    for (i = 0; i < n_active; i++) 
    {
        ac[i]->sv_position = ac[i]->grid_position;
    }
}

bool print_grid_and_check_for_activity(const struct grid *the_grid, FILE *output_file, const int count, const bool binary) 
{

    struct cell_node *grid_cell = the_grid->first_cell;

    double center_x, center_y, center_z, half_face;
    double v;
    bool act = false;

    while (grid_cell != 0) 
    {

        if (grid_cell->active) 
        {

            center_x = grid_cell->center_x;
            center_y = grid_cell->center_y;
            center_z = grid_cell->center_z;

            v = grid_cell->v;
            half_face = grid_cell->half_face_length;

            if (count > 0) 
            {
                if (grid_cell->v > -86.0) 
                {
                    act = true;
                }
            } 
            else 
            {
                act = true;
            }

            if(binary) 
            {
                fwrite (&center_x, sizeof(center_x), 1, output_file);
                fwrite (&center_y, sizeof(center_y), 1, output_file);
                fwrite (&center_z, sizeof(center_z), 1, output_file);
                fwrite (&half_face, sizeof(half_face), 1, output_file);
                fwrite (&v, sizeof(v), 1, output_file);
            }
            else 
            {
                fprintf(output_file, "%g,%g,%g,%g,%g\n", center_x, center_y, center_z, half_face, v);
            }
        }
        grid_cell = grid_cell->next;
    }

    return act;
}