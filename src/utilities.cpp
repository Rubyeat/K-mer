#include "../headers/utilities.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

char *neigh_strdup(const char *src) {
    char *dst = NULL;

    if (src != NULL) {
        size_t length = strlen(src);

        dst = (char*)malloc(length + 1);

        if (dst != NULL) {
            strcpy(dst, src);
        }
    }

    return dst;
}

size_t trim_trailing_space(char *s) {
    size_t length = strlen(s);
    
    while (length > 0 && isspace(s[length - 1])) {
        length--;
    }
    
    s[length] = '\0';
    
    return length;
}

size_t filename_copy(const char *path, char *dest, size_t size) {
    size_t length = strlen(path);
    size_t end = length;
    
    while (end > 0 && path[end - 1] != '.') {
        end--;
    }
    
    if (end == 0) {
        /* No extension found */
        end = length;
    } else {
        /* Remove the trailing dot */
        end--;
    }
    
    size_t start = end;
    
    while (start > 0 && path[start - 1] != '/') {
        start--;
    }
    
    size_t count = end - start;
    
    if (dest != NULL) {
        size_t n = ((size - 1) < count) ? (size - 1) : count;
        
        memcpy(dest, path + start, n);
        dest[n] = '\0';
    }
    
    return count;
}

dist_matrix *load_file(vector<Species> species, double **ar) {
    int result;
    uint32_t species_count;

    species_count = species.size();
    dist_matrix *dmat = dist_matrix_init(species_count);

    if (!dmat) {
        printf("Unable to create distance matrix");
        return NULL;
    }

    for (uint32_t i = 0; i < species_count; i++) {
        /* species name: up to 30 alphabetic or whitespace characters */
        char species_name[31];
        int a = 0;
        for(; a < species[i].name.length(); a++) {
            species_name[a] = species[i].name[a];
        }
        species_name[a] = 0;

        dist_matrix_set_species_name(dmat, i, species_name);
        dmat->cluster_sizes[i] = 1;

        for (uint32_t j = 0; j < i; j++) {
            double *element = dist_matrix_element(dmat, i, j);
            *element = ar[i][j];

        }
    }
    return dmat;
}
