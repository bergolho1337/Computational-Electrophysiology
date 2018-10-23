#include "../include/options.h"

Options::Options (int argc, char *argv[])
{
    string str;
    ifstream in_file(argv[1]);
    if (!in_file)
    {
        cerr << "[Options] ERROR! Reading input file!" << endl;
        exit(EXIT_FAILURE);
    }

    in_file >> num_threads;
    in_file >> dt;
    in_file >> tmax;

    in_file >> str;
    if (str == "yes")
        steady_state = true;
    else if (str == "no")
        steady_state = false;
    else
    {
        cerr << "[Options] ERROR! Invalid parameter!" << endl;
        exit(EXIT_FAILURE);
    }

    in_file >> print_rate;
    in_file >> sst_rate;
    in_file >> network_filename;
    in_file >> sst_filename;
    in_file >> plot_filename;

    in_file >> start_h;
    in_file >> num_div_cell;
    in_file >> start_diameter;
    in_file >> sigma_c;
    in_file >> G_gap;

    in_file >> stim_current;
    in_file >> stim_start;
    in_file >> stim_duration;
    in_file >> start_period;
    in_file >> end_period;
    in_file >> period_step;
    in_file >> n_cycles;
    in_file >> id_limit;

    in_file.close();

}

void Options::print_user_options ()
{
    cout << "[main]" << endl;
    cout << "num_threads = " << num_threads << endl;
    cout << "dt = " << dt << endl;
    cout << "tmax = " << tmax << endl;
    cout << "Use Steady-State = " << steady_state << endl;
    cout << "Print rate = " << print_rate << endl;
    cout << "Steady-State rate = " << sst_rate << endl;
    cout << "Network filename = " << network_filename << endl;
    cout << "Steady-State filename = " << sst_filename << endl;
    cout << "Plot filename = " << plot_filename << endl << endl;
    
    cout << "[cell]" << endl;
    cout << "start_h = " << start_h << endl;
    cout << "num_div_cell = " << num_div_cell << endl;
    cout << "start diameter = " << start_diameter << endl;
    cout << "sigma_c = " << sigma_c << endl;
    cout << "G_gap = " << G_gap << endl << endl;

    cout << "[stimulus]" << endl;
    cout << "stim_current = " << stim_current << endl;
    cout << "stim_start = " << stim_start << endl;
    cout << "stim_duration = " << stim_duration << endl;
    cout << "start_period = " << start_period << endl;
    cout << "end_period = " << end_period << endl;
    cout << "period_step = " << period_step << endl;
    cout << "n_cycles = " << n_cycles << endl;
    cout << "id_limit = " << id_limit << endl;


}