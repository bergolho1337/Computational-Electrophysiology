#include "monodomain_solver.h"
#include "ode_solver.h"
#include <cassert>
#include <cinttypes>
#include <vector>
#include "../utils/logfile_utils.h"
#include "../utils/stop_watch.h"
#include "../purkinje/purkinje.h"
//#include "config/purkinje_config.h"

static inline double ALPHA (double beta, double cm, double dt, double h) 
{
    return (((beta * cm) / dt) * UM2_TO_CM2) * pow (h, 3.0);
    //return (((beta * cm) / dt)) * pow (h, 3.0);
}

struct monodomain_solver *new_monodomain_solver ()
{
    struct monodomain_solver *result = (struct monodomain_solver*)malloc(sizeof(struct monodomain_solver));

    result->beta = 0.14;
    result->cm = 1.2;

    return result;
}

void configure_monodomain_solver_from_options (struct monodomain_solver *the_monodomain_solver,
                                               struct user_options *options) 
{

    assert (the_monodomain_solver);
    assert (options);

    the_monodomain_solver->num_threads = options->num_threads;
    the_monodomain_solver->final_time = options->final_time;

    the_monodomain_solver->dt = options->dt;

    the_monodomain_solver->sigma = options->sigma;
    the_monodomain_solver->beta = options->beta;
    the_monodomain_solver->cm = options->cm;
}

void solve_monodomain (struct monodomain_solver *monodomain_solver, struct ode_solver *ode_solver,\
                        struct grid *grid, struct user_options *configs)
{
    assert (configs);
    assert (grid);
    assert (monodomain_solver);
    assert (ode_solver);

    print_to_stdout_and_file (LOG_LINE_SEPARATOR);

    long ode_total_time = 0, cg_total_time = 0, total_write_time = 0, total_mat_time = 0, total_ref_time = 0,
         total_deref_time = 0, cg_partial, total_config_time = 0;

    uint32_t total_cg_it = 0;

    struct stop_watch solver_time, ode_time, cg_time, part_solver, part_mat, write_time, ref_time, deref_time,
        config_time;

    init_stop_watch (&config_time);

    start_stop_watch (&config_time);

    ///////MAIN CONFIGURATION BEGIN//////////////////
    init_ode_solver_with_cell_model(ode_solver);
    struct stim_config_hash *stimuli_configs = configs->stim_configs;
    struct purkinje_config *pk_config = configs->purkinje_config;

    double last_stimulus_time = -1.0;

    if (stimuli_configs) 
    {
        // Init all stimuli
        STIM_CONFIG_HASH_FOR_EACH_KEY_APPLY_FN_IN_VALUE_AND_KEY (stimuli_configs, init_stim_functions);

        //Find last stimuli
        size_t s_size = stimuli_configs->size;
        double s_end;
        for (int i = 0; i < s_size; i++) 
        {
            for (struct stim_config_elt *e = stimuli_configs->table[i % s_size]; e != 0; e = e->next) 
            {
                s_end = e->value->stim_start + e->value->stim_duration;
                if(s_end > last_stimulus_time) last_stimulus_time = s_end;
            }
        }
    }
    // Configure the functions and set the Purkinje mesh domain
    if (pk_config) 
    {
        set_spatial_purkinje(pk_config,grid);
    } 
    else 
    {
        print_to_stdout_and_file ("No Purkinje configuration provided! Exiting!\n");
        exit (EXIT_FAILURE);
    }

    // Assemble the matrix of the PDE
    Eigen::SparseMatrix<double> A = assembly_matrix(monodomain_solver,grid,pk_config);
    
}

Eigen::SparseMatrix<double> assembly_matrix (struct monodomain_solver *monodomain_solver, struct grid *grid, struct purkinje_config *pk_config) 
{
    struct graph *pk_net = grid->the_purkinje_network;
    int n = pk_net->total_nodes;
    double cm = monodomain_solver->cm;
    double beta = monodomain_solver->beta;
    double sigma = monodomain_solver->sigma;
    double dt = monodomain_solver->dt;
    double h = pk_config->start_h;

    Eigen::SparseMatrix<double> A(n,n);
    std::vector< Eigen::Triplet<double> > coeff;

    double alpha = ALPHA(beta,cm,dt,h);
    
    struct node *ptr = pk_net->list_nodes;
    while (ptr != NULL)
    {
        int u = ptr->id;
        double value = -sigma * h;
        
        struct edge *ptrl = ptr->list_edges;
        while (ptrl != NULL)
        {
            int v = ptrl->dest->id;
            coeff.push_back( Eigen::Triplet<double>(u,v,value) );
            ptrl = ptrl->next;
        }
        value = (ptr->num_edges * sigma * h) + alpha;
        coeff.push_back( Eigen::Triplet<double>(u,u,value) );
        ptr = ptr->next;
    }
    A.setFromTriplets(coeff.begin(),coeff.end());
    A.makeCompressed();

    return A;
}