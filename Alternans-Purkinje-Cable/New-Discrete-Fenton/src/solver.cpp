#include "../include/solver.h"

Solver::Solver (Options *user_options)
{
    nthreads = user_options->num_threads;
    dt = user_options->dt;
    tmax = user_options->tmax;
    print_rate = user_options->print_rate;
    sst_rate = user_options->sst_rate;
    network_filename = user_options->network_filename;
    steady_filename = user_options->sst_filename;
    plot_filename = user_options->plot_filename;

    dx = user_options->start_h / user_options->num_div_cell;
    M = nearbyint(tmax/dt);

    set_plot_cells();

    stim_config = new Stimulus(user_options);
    //stim_config->print();

    
}

void Solver::set_plot_cells ()
{
    // Allocate and initialize the vector 'ids' with the volumes to be plotted from the input file .plt
    plot = (Plot*)malloc(sizeof(Plot));
    FILE *pltFile = fopen(plot_filename.c_str(),"r");
    if (!pltFile)
    {
        cerr << "[-] ERROR! Cannot open PLT file!" << endl;
        exit(EXIT_FAILURE);
    }

    if (!fscanf(pltFile,"%d",&plot->np)) 
    {
        cerr << "[-] ERROR! Reading PLT file!" << endl;
        exit(EXIT_FAILURE);
    }
    plot->ids = (int*)malloc(sizeof(int)*plot->np);

    for (int i = 0; i < plot->np; i++)
    {
        if (!fscanf(pltFile,"%d",&plot->ids[i]))
        {
            cerr << "[-] ERROR! Reading PLT file!" << endl;
            exit(EXIT_FAILURE);
        }
    }

    fclose(pltFile);
}