#include "options.h"

User_Options::User_Options (int argc, char *argv[])
{
    string str;
    ifstream in_file(argv[1]);
    
    in_file >> str;
    if (str == "t")
        steady_state = true;
    else if (str == "s")
        steady_state = false;
    else
    {
        cerr << "[-] ERROR! Invalid parameter!" << endl;
        exit(EXIT_FAILURE);
    }

    in_file >> dt >> tmax;
    in_file >> mesh_filename >> start_h;
    in_file >> sst_filename >> plot_filename;
    in_file >> alfa >> diameter >> sigma_c;

    in_file.close();

}

void User_Options::print_user_options ()
{
    cout << "Steady state = " << steady_state << endl;
    cout << "Dt = " << dt << endl;
    cout << "tmax = " << tmax << endl;
    cout << "Mesh filename = " << mesh_filename << endl;
    cout << "Steady state filename = " << sst_filename << endl;
    cout << "Plot filename = " << plot_filename << endl;
    cout << "alpha = " << alfa << endl;
    cout << "diameter = " << diameter << endl;
    cout << "sigma_c = " << sigma_c << endl;
}