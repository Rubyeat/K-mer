#include "../headers/phylogenetic_tree.h"
#include "../headers/utilities.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>

string newick_tree;
size_t btree_storage_size(uint32_t nodes_count) {
    return sizeof(btree_storage) + (nodes_count * member_size(btree_storage, nodes[0]));
}

btree_storage *btree_storage_init(uint32_t nodes_count) {
    size_t size = btree_storage_size(nodes_count);

    /* Zero-initialize the storage to ensure that btree_storage_free can correctly
     * free the strings */
    btree_storage *storage = (btree_storage*) calloc(size, 1);

    if (storage != NULL) {
        storage->nodes_count = nodes_count;
    }
    
    return storage;
}

void btree_storage_free(btree_storage *storage) {
    if (storage) {
        for (uint32_t i = 0; i < storage->nodes_count; ++i) {
            free(storage->nodes[i].node_name);
        }

        free(storage);
    }
}

btree_node *btree_storage_fetch(btree_storage *storage) {
    assert(storage->used_nodes < storage->nodes_count);
    
    uint32_t idx = storage->used_nodes;
    btree_node *node = &storage->nodes[idx];
    
    storage->used_nodes++;
    
    return node;
}

char *btree_node_set_name(btree_node *node, const char *name) {
    free(node->node_name);

    node->node_name = neigh_strdup(name);

    return node->node_name;
}

uint32_t btree_get_height(btree_node *root) {
    if (root == NULL) {
        return 0;
    }
    
    uint32_t left = btree_get_height(root->left);
    uint32_t right = btree_get_height(root->right);

    uint32_t height = 1;
    height += (left >= right) ? left : right;
    
    return height;
}

static void btree_print_node(btree_node *root, double distance, uint32_t depth, bool *is_open) {
    if (root == NULL) {
        return;
    }

    btree_print_node(root->right, root->distance_right, depth + 1, is_open);

    if (depth > 0) {
        for (uint32_t i = 1; i < depth; i++) {
            printf("  %c        ", is_open[i] ? '|' : ' ');
        }
    
        printf("  |--%.2lf-- ", distance);
    }

    is_open[depth] = !is_open[depth];

    printf("%s\n", root->node_name);

    btree_print_node(root->left, root->distance_left, depth + 1, is_open);
}

void btree_newick_format(btree_node *root, double distance, uint32_t depth, bool *is_open, bool isLeft, bool isFirstCall) {
    if (root == NULL) {
        return;
    }

    bool is_node = true;
    char* p;
    long converted = strtol(root->node_name, &p, 10);
    if (*p) {
        is_node = false;
    }
    if(is_node) {
        newick_tree += "(" ;
    }
    btree_newick_format(root->right, root->distance_right, depth + 1, is_open,  false, false);

    /*if (depth > 0) {
        for (uint32_t i = 1; i < depth; i++) {
           printf("  %c        ", is_open[i] ? '|' : ' ');
        }

        printf("  |--%.2lf-- ", distance);
    }*/
    is_open[depth] = !is_open[depth];
   /* printf("%s\n", root->node_name);*/

    if(!isLeft && !is_node) {
        newick_tree +=  std::string(root->node_name) + ": " + std::to_string(distance) + ",";
    }
    btree_newick_format(root->left, root->distance_left, depth + 1, is_open, true, false);
    if(isLeft) {
        if (!is_node) {
            newick_tree += std::string(root->node_name) + ":" + std::to_string(distance) ;
        }
    }
    if (is_node ) {
        if(isFirstCall) {
            newick_tree +=  ")";
        } else {
            newick_tree +=  "):" + std::to_string(distance);
        }
        if(!isLeft && !isFirstCall) {
            newick_tree += ",";
        }
    }


}

void btree_print_tree(btree_node *root) {
    assert(root != NULL);
    
    uint32_t height = btree_get_height(root);
    
    bool is_open[height];
    
    for (uint32_t i = 0; i < height; i++) {
        is_open[i] = false;
    }
    
    btree_print_node(root, 0, 0, is_open);
}

string btree_print_newick_tree(btree_node *root) {
    assert(root != NULL);
    newick_tree = "";
    uint32_t height = btree_get_height(root);

    bool is_open[height];

    for (uint32_t i = 0; i < height; i++) {
        is_open[i] = false;
    }

    btree_newick_format(root, 0, 0, is_open, false, true);
    return newick_tree;
}

void btree_print_trees(btree_node **trees, uint32_t tree_count) {
    for (uint32_t i = 0; i < tree_count; i++) {
        btree_print_tree(trees[i]);
        printf("\n");
    }
}
