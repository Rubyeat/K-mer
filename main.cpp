#include <iostream>
#include <cstdio>
#include "headers/species.h"
using namespace std;

vector<string> speciesName = {"armadillo", "baboon", "cat", "chicken", "chimp", "clouded_leopard", "colobus_monkey", "cow", "cpbat", "dog", "dusky_titi","elephant", "ferret", "fugu", "galago", "gibbon", "gorilla", "guinea_pig", "heagehog","horse","human", "lemur", "macaque", "marmoset", "monodelphis", "mouse", "mouse_lemur", "muntjak_indian", "oppssum", "orangatun", "owl_monkey", "pig", "platypus", "rabbit", "rat", "rfbat", "sheep", "shrew","squirrel_monkey", "st_squirrel", "tenrec", "tetra", "vervet", "wallaby"};
vector<string> file1Name = {"armadillo_counts_dumps.fa", "baboon_counts_dumps.fa", "cat_counts_dumps.fa", "chicken_counts_dumps.fa", "chimp_counts_dumps.fa", "clouded_leopard_counts_dumps.fa", "colobus_monkey_counts_dumps.fa", "cow_counts_dumps.fa", "cpbat_counts_dumps.fa", "dog_counts_dumps.fa", "dusky_titi_counts_dumps.fa","elephant_counts_dumps.fa", "ferret_counts_dumps.fa", "fugu_counts_dumps.fa", "galago_counts_dumps.fa", "gibbon_counts_dumps.fa", "gorilla_counts_dumps.fa", "guinea_pig_counts_dumps.fa", "heagehog_counts_dumps.fa","horse_counts_dumps.fa","human_counts_dumps.fa", "lemur_counts_dumps.fa", "macaque_counts_dumps.fa", "marmoset_counts_dumps.fa", "monodelphis_counts_dumps.fa", "mouse_counts_dumps.fa", "mouse_lemur_counts_dumps.fa", "muntjak_indian_counts_dumps.fa", "oppssum_counts_dumps.fa", "orangatun_counts_dumps.fa", "owl_monkey_counts_dumps.fa", "pig_counts_dumps.fa", "platypus_counts_dumps.fa", "rabbit_counts_dumps.fa", "rat_counts_dumps.fa", "rfbat_counts_dumps.fa", "sheep_counts_dumps.fa", "shrew_counts_dumps.fa", "squirrel_monkey_counts_dumps.fa", "st_squirrel_counts_dumps.fa", "tenrec_counts_dumps.fa", "tetra_counts_dumps.fa", "vervet_counts_dumps.fa", "wallaby_counts_dumps.fa"};
vector<string> file2Name = {"armadillo.fa", "baboon.fa", "cat.fa", "chicken.fa", "chimp.fa", "clouded_leopard.fa", "colobus_monkey.fa", "cow.fa", "cpbat.fa", "dog.fa", "dusky_titi.fa","elephant.fa", "ferret.fa", "fugu.fa", "galago.fa", "gibbon.fa", "gorilla.fa", "guinea_pig.fa", "heagehog.fa","horse.fa","human.fa", "lemur.fa", "macaque.fa", "marmoset.fa", "monodelphis.fa", "mouse.fa", "mouse_lemur.fa", "muntjak_indian.fa", "oppssum.fa", "orangatun.fa", "owl_monkey.fa", "pig.fa", "platypus.fa", "rabbit.fa", "rat.fa", "rfbat.fa", "sheep.fa", "shrew.fa", "squirrel_monkey.fa", "st_squirrel.fa", "tenrec.fa", "tetra.fa", "vervet.fa", "wallaby.fa"};
int main() {
    vector<Species> species;
    double distanceMatrixEucledian[speciesName.size()][speciesName.size()];
    double distanceMatrixMahalnobis[speciesName.size()][speciesName.size()];
    double distanceMatrixFactorcount[speciesName.size()][speciesName.size()];
    for (int i = 0; i < speciesName.size() ; i++) {
        Species s(speciesName[i], 20);
        s.init("../files/Data20/" + file1Name[i], "../files/sequences/" + file2Name[i]);
        species.push_back(s);
    }
    for (int j = 0; j < speciesName.size(); j++) {
        for (int i = 0; i < speciesName.size(); i++) {
            if(j == i){
                distanceMatrixEucledian[j][i] = 0;
                continue;
            } else if(i < j) {
                distanceMatrixEucledian[j][i] = distanceMatrixEucledian[i][j];
            } else {
                distanceMatrixEucledian[j][i] = species[j].calculateEuclidianDistance(species[i]);
            }
        }
    }
    for (int j = 0; j < speciesName.size(); j++) {
        for (int i = 0; i < speciesName.size(); i++) {
            if(j == i){
                distanceMatrixMahalnobis[j][i] = 0;
                continue;
            } else if(i < j) {
                distanceMatrixMahalnobis[j][i] = distanceMatrixMahalnobis[i][j];
            } else {
                distanceMatrixMahalnobis[j][i] = species[j].calculateMahalnobisDistance(species[i]);
            }
        }
    }
    for (int j = 0; j < speciesName.size(); j++) {
        for (int i = 0; i < speciesName.size(); i++) {
            if(j == i){
                distanceMatrixFactorcount[j][i] = 0;
                continue;
            } else if(i < j) {
                distanceMatrixFactorcount[j][i] = distanceMatrixFactorcount[i][j];
            } else {
                distanceMatrixFactorcount[j][i] = species[j].calculateFactorDistance(species[i]);
            }
        }
    }

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
    //fclose (stdout);

    return 0;
}