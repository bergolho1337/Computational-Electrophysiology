#include <iostream>
#include <cstdio>

using namespace std;

const int MAX_SIZE = 500;

void write_sst_config_file (const int period, const int start_period, const int step_period)
{
    char filename[MAX_SIZE];
    sprintf(filename,"files/sst_beeler_reuter_8cm_NoGap_%dms.ini",period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=1\n");
    fprintf(file,"dt=0.01\n");
    fprintf(file,"simulation_time=%d\n",20*period);
    if (period == start_period)
        fprintf(file,"use_steady_state=no\n");
    else
        fprintf(file,"use_steady_state=yes\n");
    fprintf(file,"print_rate=1000000\n");
    fprintf(file,"sst_rate=%d\n",2000*period);
    fprintf(file,"network_filename = networks/cable-8cm.vtk\n");
    if (period == start_period)
        fprintf(file,"sst_filename = teste.sst\n");
    else
        fprintf(file,"sst_filename = steady_state/cable-BR-8cm-%dms-NoGap.sst\n",period+step_period);
    fprintf(file,"plot_filename = plot/cable-8cm-BR.plt\n\n");
    fprintf(file,"[cell]\n");
    fprintf(file,"start_h = 0.025\n");
    fprintf(file,"num_div_cell = 1\n");
    fprintf(file,"start_diameter = 0.003\n");
    fprintf(file,"sigma_c = 0.001334\n");
    fprintf(file,"G_gap = 0.628\n");
    fprintf(file,"library_file=shared_libs/libbeeler_reuter_1977.so\n\n");
    fprintf(file,"[stim_feedback]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 1.0\n");
    fprintf(file,"stim_current = 2.0\n");
    fprintf(file,"n_cycles=20\n");
    fprintf(file,"start_period=%d\n",period);
    fprintf(file,"end_period=%d\n",period);
    fprintf(file,"period_step=100\n");
    fprintf(file,"id_limit=5\n");
    fprintf(file,"function=stim_if_id_less_than\n");
    
    fclose(file);
}

void write_simulation_config_file (const int period, const int start_period)
{
    char filename[MAX_SIZE];
    sprintf(filename,"files/simple_beeler_reuter_8cm_NoGap_%dms.ini",period);
    FILE *file = fopen(filename,"w+");

    fprintf(file,"[main]\n");
    fprintf(file,"num_threads=2\n");
    fprintf(file,"dt=0.01\n");
    fprintf(file,"simulation_time=%d\n",2*period);
    fprintf(file,"use_steady_state=yes\n");
    fprintf(file,"print_rate=10\n");
    fprintf(file,"sst_rate=%d\n",2000*period);
    fprintf(file,"network_filename = networks/cable-8cm.vtk\n");
    fprintf(file,"sst_filename = steady_state/cable-BR-8cm-%dms-NoGap.sst\n",period);
    fprintf(file,"plot_filename = plot/cable-8cm-BR.plt\n\n");
    fprintf(file,"[cell]\n");
    fprintf(file,"start_h = 0.025\n");
    fprintf(file,"num_div_cell = 1\n");
    fprintf(file,"start_diameter = 0.003\n");
    fprintf(file,"sigma_c = 0.001334\n");
    fprintf(file,"G_gap = 0.628\n");
    fprintf(file,"library_file=shared_libs/libbeeler_reuter_1977.so\n\n");
    fprintf(file,"[stim_feedback]\n");
    fprintf(file,"stim_start = 0.0\n");
    fprintf(file,"stim_duration = 1.0\n");
    fprintf(file,"stim_current = 2.0\n");
    fprintf(file,"n_cycles=20\n");
    fprintf(file,"start_period=%d\n",period);
    fprintf(file,"end_period=%d\n",period);
    fprintf(file,"period_step=100\n");
    fprintf(file,"id_limit=5\n");
    fprintf(file,"function=stim_if_id_less_than\n");
    
    fclose(file);
}

void write_config_files (const int period, const int start_period, const int step_period)
{
    write_sst_config_file(period,start_period,step_period);
    write_simulation_config_file(period,start_period);
}

int main ()
{
    int start_period = 400;
    int end_period = 200;
    int period_step = 5;
    for (int period = start_period; period >= end_period; period -= period_step)
    {
        printf("Writing %dms files ...\n",period);
        write_config_files(period,start_period,period_step);
    }
    return 0;
}
