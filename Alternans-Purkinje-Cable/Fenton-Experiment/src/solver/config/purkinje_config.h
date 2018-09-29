#ifndef PURKINJE_CONFIG_H_
#define PURKINJE_CONFIG_H_

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

struct purkinje_config
{
    char *network_filename;
    bool network_filename_was_set;
    char *name;
    bool name_was_set;
    double start_h;
    bool start_h_was_set;
};

struct purkinje_config* new_purkinje_config ();

#endif