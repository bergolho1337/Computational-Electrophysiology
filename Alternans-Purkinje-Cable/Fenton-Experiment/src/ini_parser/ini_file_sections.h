//
// Created by sachetto on 12/10/17.
//

#ifndef FENTON_INI_FILE_HEADERS_H
#define FENTON_INI_FILE_HEADERS_H

#define MAIN_SECTION "main"
#define CELL_SECTION "cell"
#define STIM_SECTION "stim"

#define MATCH_SECTION_AND_NAME(s, n) strcmp(section, s) == 0 && strcmp(name, n) == 0
#define MATCH_SECTION(s) strcmp(section, s) == 0
#define MATCH_NAME(v) strcmp(name, v) == 0
#define SECTION_STARTS_WITH(s) strncmp(section, s, strlen(s)) == 0
#define A_STARTS_WITH_B(a ,b) strncmp(a, b, strlen(b)) == 0


#endif // FENTON_INI_FILE_HEADERS_H
