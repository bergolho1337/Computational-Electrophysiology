#ifndef _USER_CONFIG_H_
#define _USER_CONFIG_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

using namespace std;

#define PRINT_LINE "--------------------------------------------------------------------------------"
#define PRINT_LINE2 "==================================================================================="
#define MAX_FILENAME_SIZE 100

struct user_options
{
    double start_h;
    double dt;
    double tmax;
    double sigma_c;
    double G_gap;
    double diameter;
    int num_div_cell;
    int use_steady_state;
    char pk_network_filename[MAX_FILENAME_SIZE];
    char plot_filename[MAX_FILENAME_SIZE];
};

struct user_options* new_user_options (int argc, char *argv[]);
void free_user_options (struct user_options *configs);
void print_user_options (struct user_options *configs);

void display_usage (const char p_name[]);

#endif
