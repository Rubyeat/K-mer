#include <iostream>
#include <cstdio>
#include "headers/dist_matrix.h"
#include "headers/phylogenetic_tree.h"
#include "headers/neighbour_joining.h"
#include "headers/utilities.h"
using namespace std;

vector<string> speciesName = {"armadillo", "baboon", "cat", "chicken", "chimp", "clouded_leopard", "colobus_monkey", "cow", "dog", "dusky_titi","elephant", "ferret", "galago", "gibbon", "gorilla", "guinea_pig", "heagehog","horse","human", "macaque", "marmoset",  "mouse", "mouse_lemur", "orangatun", "owl_monkey", "pig", "platypus", "rabbit", "rat", "sheep", "shrew","squirrel_monkey", "st_squirrel", "tenrec",  "vervet", "wallaby"};
vector<string> file1Name = {"armadillo_counts_dumps.fa", "baboon_counts_dumps.fa", "cat_counts_dumps.fa", "chicken_counts_dumps.fa", "chimp_counts_dumps.fa", "clouded_leopard_counts_dumps.fa", "colobus_monkey_counts_dumps.fa", "cow_counts_dumps.fa", "dog_counts_dumps.fa", "dusky_titi_counts_dumps.fa","elephant_counts_dumps.fa", "ferret_counts_dumps.fa", "galago_counts_dumps.fa", "gibbon_counts_dumps.fa", "gorilla_counts_dumps.fa", "guinea_pig_counts_dumps.fa", "heagehog_counts_dumps.fa","horse_counts_dumps.fa","human_counts_dumps.fa", "macaque_counts_dumps.fa", "marmoset_counts_dumps.fa", "mouse_counts_dumps.fa", "mouse_lemur_counts_dumps.fa", "orangatun_counts_dumps.fa", "owl_monkey_counts_dumps.fa", "pig_counts_dumps.fa", "platypus_counts_dumps.fa", "rabbit_counts_dumps.fa", "rat_counts_dumps.fa", "sheep_counts_dumps.fa", "shrew_counts_dumps.fa", "squirrel_monkey_counts_dumps.fa", "st_squirrel_counts_dumps.fa", "tenrec_counts_dumps.fa", "vervet_counts_dumps.fa", "wallaby_counts_dumps.fa"};
vector<string> file2Name = {"armadillo.fa", "baboon.fa", "cat.fa", "chicken.fa", "chimp.fa", "clouded_leopard.fa", "colobus_monkey.fa", "cow.fa", "dog.fa", "dusky_titi.fa","elephant.fa", "ferret.fa", "galago.fa", "gibbon.fa", "gorilla.fa", "guinea_pig.fa", "heagehog.fa","horse.fa","human.fa", "macaque.fa", "marmoset.fa", "mouse.fa", "mouse_lemur.fa", "orangatun.fa", "owl_monkey.fa", "pig.fa", "platypus.fa", "rabbit.fa", "rat.fa", "sheep.fa", "shrew.fa", "squirrel_monkey.fa", "st_squirrel.fa", "tenrec.fa", "vervet.fa", "wallaby.fa"};
void printMatrix(vector<Species> species, double **distanceMatrixEucledian, double **distanceMatrixMahalnobis, double **distanceMatrixFactorcount);
int main() {
    vector<Species> species;
    double** distanceMatrixEucledian = new double*[speciesName.size()];
    for (int i = 0; i < speciesName.size(); ++i) {
        distanceMatrixEucledian[i] = new double[speciesName.size()];
    }
    double** distanceMatrixMahalnobis = new double*[speciesName.size()];;
    for (int i = 0; i < speciesName.size(); ++i) {
        distanceMatrixMahalnobis[i] = new double[speciesName.size()];
    }
    double** distanceMatrixFactorcount = new double*[speciesName.size()];;
    for (int i = 0; i < speciesName.size(); ++i) {
        distanceMatrixFactorcount[i] = new double[speciesName.size()];
    }
    for (int i = 0; i < speciesName.size() ; i++) {
        Species s(speciesName[i], 20);
        s.init("../files/Data20/" + file1Name[i], "../files/sequences/" + file2Name[i]);
        species.push_back(s);
    }
    for (int j = 0; j < speciesName.size(); j++) {
        for (int i = 0; i < speciesName.size(); i++) {
            if(j == i){
                distanceMatrixEucledian[j][i] = 0;
                distanceMatrixMahalnobis[j][i] = 0;
                distanceMatrixFactorcount[j][i] = 0;
                continue;
            } else if(i < j) {
                distanceMatrixEucledian[j][i] = distanceMatrixEucledian[i][j];
                distanceMatrixMahalnobis[j][i] = distanceMatrixMahalnobis[i][j];
                distanceMatrixFactorcount[j][i] = distanceMatrixFactorcount[i][j];
            } else {
                vector<double> v = species[j].calculateDistance(species[i]);
                distanceMatrixEucledian[j][i] = v[0];
                distanceMatrixMahalnobis[j][i] = v[1];
                distanceMatrixFactorcount[j][i] = v[2];
            }
        }
    }
    //Construct tree
    dist_matrix *dmatE = load_file(species, distanceMatrixFactorcount);
    double u[dmatE->species_count];
    char cluster_name[2 + 3 * sizeof(dmatE->species_count) + 1];
    uint32_t cluster_id = 1;

    btree_node *partial_trees[dmatE->species_count];
    btree_storage *tree_storage = nj_tree_init(dmatE, partial_trees);

    while (dmatE->species_count >= 2) {
        /* Compute the average distance of each clusters from the others */
        dist_matrix_compute_avg_distances(dmatE, u);

        /* Find the pair of nearest clusters */
        uint32_t c1, c2;
        nj_find_nearest_clusters(dmatE, u, &c1, &c2);

        /* Generate a name for the new cluster */
        unsigned long result = snprintf(cluster_name, sizeof(cluster_name),"%d", cluster_id);
        /* Add a node for the new cluster to the array of partial trees */
        nj_tree_add_node(dmatE, u, tree_storage, partial_trees, cluster_name, c1, c2);

        /* Create a new dist_matrix joining the specified clusters */
        dist_matrix *joined = nj_join_clusters(dmatE, cluster_name, c1, c2);

        if (joined == NULL) {
            /* Error, stop here */
            break;
        }

        /* Release the old distance matrix */
        dist_matrix_free(dmatE);
        dmatE = joined;

        cluster_id++;
    }

    btree_node *phyl_tree = partial_trees[0];

    btree_print_tree(phyl_tree);
    printf("\n\n");

    dist_matrix_free(dmatE);
    btree_storage_free(tree_storage);
    printMatrix(species, distanceMatrixEucledian, distanceMatrixFactorcount, distanceMatrixFactorcount);
    delete[] distanceMatrixEucledian;
    delete[] distanceMatrixMahalnobis;
    delete[] distanceMatrixFactorcount;
    return 0;
}

void printMatrix(vector<Species> species, double **distanceMatrixEucledian, double **distanceMatrixMahalnobis, double **distanceMatrixFactorcount) {
    string x = ("../files/Output/matrixEucledian_" + to_string(species[0].kmer) + ".txt");
    //freopen (x.c_str(),"w",stdout);
    printf("%d\n", speciesName.size());
    for (int k = 0; k < speciesName.size(); ++k) {
        cout << species[k].name << " ";
        for (int i = 0; i <speciesName.size() ; ++i) {
            printf("%.2lf ", distanceMatrixEucledian[k][i]);
        }
        printf("\n");
    }
    //fclose (stdout);
    string xx = ("../files/Output/matrixMahalnobis_" + to_string(species[0].kmer) + ".txt");
    //freopen (xx.c_str(),"w",stdout);
    printf("%d\n", speciesName.size());
    for (int k = 0; k < speciesName.size(); ++k) {
        cout << species[k].name << " ";
        for (int i = 0; i <speciesName.size() ; ++i) {
            printf("%.2lf ", distanceMatrixMahalnobis[k][i]);
        }
        printf("\n");
    }
    //fclose (stdout);
    string xxx = ("../files/Output/matrixFactorCount_" + to_string(species[0].kmer) + ".txt");
    //freopen (xxx.c_str(),"w",stdout);
    printf("%d\n", speciesName.size());
    for (int k = 0; k < speciesName.size(); ++k) {
        cout << species[k].name << " ";
        for (int i = 0; i <speciesName.size() ; ++i) {
            printf("%.2lf ", distanceMatrixFactorcount[k][i]);
        }
        printf("\n");
    }
}