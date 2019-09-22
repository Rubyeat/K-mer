#ifndef NEIGH_UTILITIES_H
#define NEIGH_UTILITIES_H

#include <stddef.h>
#include <vector>
#include "species.h"
#include "dist_matrix.h"

using namespace std;

#define member_size(type, member) sizeof(((type *)0)->member)

char *neigh_strdup(const char *s);

size_t trim_trailing_space(char *s);
size_t filename_copy(const char *path, char *dest, size_t size);
dist_matrix *load_file(vector<Species> species, double **);

#endif //NEIGH_UTILITIES_H
