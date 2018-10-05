#include "user_config.h"

struct user_options* new_user_options (int argc, char *argv[])
{
    struct user_options *config = (struct user_options*)malloc(sizeof(struct user_options));
    
    ifstream in(argv[1]);
    in >> config->dt;
    in >> config->tmax;
    in >> config->start_h;
    in >> config->diameter;
    in >> config->num_div_cell;
    in >> config->sigma_c;
    in >> config->G_gap;
    in >> config->use_steady_state;
    in >> config->pk_network_filename;
    in >> config->plot_filename;
    in.close();

    return config;
}

void free_user_options (struct user_options *configs)
{
    free(configs);
}

void print_user_options (struct user_options *configs)
{
    cout << PRINT_LINE2 << endl;
    cout << "dt = " << configs->dt << endl;
    cout << "tmax = " << configs->tmax << endl;
    cout << "start_h = " << configs->start_h << endl;
    cout << "diameter = " << configs->diameter << endl;
    cout << "num_div_cell = " << configs->num_div_cell << endl;
    cout << "sigma_c = " << configs->sigma_c << endl;
    cout << "G_gap = " << configs->G_gap << endl;
    cout << "use_steady_state = " << configs->use_steady_state << endl;
    cout << "purkinje_network_filename = " << configs->pk_network_filename << endl;
    cout << "plot_filename = " << configs->plot_filename << endl;
    cout << PRINT_LINE2 << endl;
}

void display_usage (const char p_name[])
{
    cout << PRINT_LINE << endl;
    cout << "Usage:> " << p_name << " <input_filename>" << endl;
    cout << PRINT_LINE2 << endl;
    cout << "<input_filename> = Name of the input filename with the parameter values" << endl;
    cout << PRINT_LINE << endl;
}