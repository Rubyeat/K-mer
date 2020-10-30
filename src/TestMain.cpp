#include <cstdio>
#include "../headers/dist_matrix.h"
#include "../headers/phylogenetic_tree.h"
#include "../headers/neighbour_joining.h"
#include "../headers/utilities.h"
#include "Distance.h"
using namespace std;
/* 20 Species */
/*vector<string> speciesName = {"armadillo", "baboon", "cat", "chicken", "chimp", "clouded_leopard", "colobus_monkey", "cow", "dog", "dusky_titi","elephant", "ferret", "galago", "gibbon", "gorilla", "guinea_pig", "heagehog","horse","human", "orangatun"};
vector<string> file1Name = {"armadillo_counts_dumps.fa", "baboon_counts_dumps.fa", "cat_counts_dumps.fa", "chicken_counts_dumps.fa", "chimp_counts_dumps.fa", "clouded_leopard_counts_dumps.fa", "colobus_monkey_counts_dumps.fa", "cow_counts_dumps.fa", "dog_counts_dumps.fa", "dusky_titi_counts_dumps.fa","elephant_counts_dumps.fa", "ferret_counts_dumps.fa", "galago_counts_dumps.fa", "gibbon_counts_dumps.fa", "gorilla_counts_dumps.fa", "guinea_pig_counts_dumps.fa", "heagehog_counts_dumps.fa","horse_counts_dumps.fa","human_counts_dumps.fa", "orangatun_counts_dumps.fa"};
vector<string> file2Name = {"armadillo.fa", "baboon.fa", "cat.fa", "chicken.fa", "chimp.fa", "clouded_leopard.fa", "colobus_monkey.fa", "cow.fa", "dog.fa", "dusky_titi.fa","elephant.fa", "ferret.fa", "galago.fa", "gibbon.fa", "gorilla.fa", "guinea_pig.fa", "heagehog.fa","horse.fa","human.fa", "orangatun.fa"};
vector<int> sequencCount = {4938920, 5082025, 4746218, 4519823, 4616057, 4578219, 5386352, 5231428, 4369232, 4686137, 4979619, 4965553, 5209548, 5528445, 4607202, 4599354, 4574284, 4643538, 4700560, 5132068, 4641652, 5032268, 5498450, 4887515, 5068389, 4825265, 5202090, 5065741, 4646332};*/

/*ecoli */
/*vector<string> speciesName = {"536a", "APEC01", "ATCC8739", "B4Sb227", "B18BS512", "BW2952", "CB9615", "CFT073", "D1Sd197", "DH10B","E24377A", "E234869", "ED1a", "EDL933", "F2a301", "F2a2457T", "F5b8401","HS","IAI1", "IAI39","MG1655","S88", "Sakai","SE11","SMS35", "SSSs046","UMN026","UTI89", "W3110"};
vector<string> file1Name = {"536a_counts_dumps.fa", "APEC01_counts_dumps.fa", "ATCC8739_counts_dumps.fa", "B4Sb227_counts_dumps.fa", "B18BS512_counts_dumps.fa", "BW2952_counts_dumps.fa", "CB9615_counts_dumps.fa", "CFT073_counts_dumps.fa", "D1Sd197_counts_dumps.fa", "DH10B_counts_dumps.fa","E24377A_counts_dumps.fa", "E234869_counts_dumps.fa", "ED1a_counts_dumps.fa", "EDL933_counts_dumps.fa", "F2a301_counts_dumps.fa", "F2a2457T_counts_dumps.fa", "F5b8401_counts_dumps.fa","HS_counts_dumps.fa","IAI1_counts_dumps.fa", "IAI39_counts_dumps.fa","MG1655_counts_dumps.fa","S88_counts_dumps.fa", "Sakai_counts_dumps.fa","SE11_counts_dumps.fa","SMS35_counts_dumps.fa", "SSSs046_counts_dumps.fa","UMN026_counts_dumps.fa","UTI89_counts_dumps.fa", "W3110_counts_dumps.fa"};
vector<int> sequencCount = {4938920, 5082025, 4746218, 4519823, 4616057, 4578219, 5386352, 5231428, 4369232, 4686137, 4979619, 4965553, 5209548, 5528445, 4607202, 4599354, 4574284, 4643538, 4700560, 5132068, 4641652, 5032268, 5498450, 4887515, 5068389, 4825265, 5202090, 5065741, 4646332};*/

/*Fish*/
/*vector<string> speciesName = {"NC_009057", "NC_009058", "NC_009059", "NC_009060", "NC_009062", "NC_009063", "NC_009064", "NC_009065", "NC_009066", "NC_009067","NC_009459", "NC_010205", "NC_011168", "NC_011169", "NC_011170", "NC_011171", "NC_011177","NC_011179","NC_012055", "NC_013564","NC_013577","NC_013663", "NC_013750","NC_018814","NC_018815"};
vector<string> file1Name = {"NC_009057_counts_dumps.fa", "NC_009058_counts_dumps.fa", "NC_009059_counts_dumps.fa", "NC_009060_counts_dumps.fa", "NC_009062_counts_dumps.fa", "NC_009063_counts_dumps.fa", "NC_009064_counts_dumps.fa", "NC_009065_counts_dumps.fa", "NC_009066_counts_dumps.fa", "NC_009067_counts_dumps.fa","NC_009459_counts_dumps.fa", "NC_010205_counts_dumps.fa", "NC_011168_counts_dumps.fa", "NC_011169_counts_dumps.fa", "NC_011170_counts_dumps.fa", "NC_011171_counts_dumps.fa", "NC_011177_counts_dumps.fa","NC_011179_counts_dumps.fa","NC_012055_counts_dumps.fa", "NC_013564_counts_dumps.fa","NC_013577_counts_dumps.fa","NC_013663_counts_dumps.fa", "NC_013750_counts_dumps.fa","NC_018814_counts_dumps.fa","NC_018815_counts_dumps.fa"};;
vector<int> sequencCount = {16626, 16639, 16603, 16572, 16587, 16598, 16703, 16649, 17045, 16510, 16441, 16793, 16544, 16556, 16543, 16976, 16486, 16457, 16508, 16657, 16798, 16627, 16628, 16588, 16590};*/
//44 species
/*vector<string> speciesName = {"armadillo", "baboon", "cat", "chicken", "chimp", "clouded_leopard", "colobus_monkey", "cow", "dog", "dusky_titi","elephant", "ferret", "galago", "gibbon", "gorilla", "guinea_pig", "heagehog","horse","human", "macaque", "marmoset",  "mouse", "mouse_lemur", "orangatun", "owl_monkey", "pig", "platypus", "rabbit", "rat", "sheep", "shrew","squirrel_monkey", "muntjak_indian", "tenrec",  "vervet", "wallaby"};
vector<string> file1Name = {"armadillo_counts_dumps.fa", "baboon_counts_dumps.fa", "cat_counts_dumps.fa", "chicken_counts_dumps.fa", "chimp_counts_dumps.fa", "clouded_leopard_counts_dumps.fa", "colobus_monkey_counts_dumps.fa", "cow_counts_dumps.fa", "dog_counts_dumps.fa", "dusky_titi_counts_dumps.fa","elephant_counts_dumps.fa", "ferret_counts_dumps.fa", "galago_counts_dumps.fa", "gibbon_counts_dumps.fa", "gorilla_counts_dumps.fa", "guinea_pig_counts_dumps.fa", "heagehog_counts_dumps.fa","horse_counts_dumps.fa","human_counts_dumps.fa", "macaque_counts_dumps.fa", "marmoset_counts_dumps.fa", "mouse_counts_dumps.fa", "mouse_lemur_counts_dumps.fa", "orangatun_counts_dumps.fa", "owl_monkey_counts_dumps.fa", "pig_counts_dumps.fa", "platypus_counts_dumps.fa", "rabbit_counts_dumps.fa", "rat_counts_dumps.fa", "sheep_counts_dumps.fa", "shrew_counts_dumps.fa", "squirrel_monkey_counts_dumps.fa", "muntjak_indian_counts_dumps.fa", "tenrec_counts_dumps.fa", "vervet_counts_dumps.fa", "wallaby_counts_dumps.fa"};
vector<string> file2Name = {"armadillo.fa", "baboon.fa", "cat.fa", "chicken.fa", "chimp.fa", "clouded_leopard.fa", "colobus_monkey.fa", "cow.fa", "dog.fa", "dusky_titi.fa","elephant.fa", "ferret.fa", "galago.fa", "gibbon.fa", "gorilla.fa", "guinea_pig.fa", "heagehog.fa","horse.fa","human.fa", "macaque.fa", "marmoset.fa", "mouse.fa", "mouse_lemur.fa", "orangatun.fa", "owl_monkey.fa", "pig.fa", "platypus.fa", "rabbit.fa", "rat.fa", "sheep.fa", "shrew.fa", "squirrel_monkey.fa", "muntjak_indian.fa", "tenrec.fa", "vervet.fa", "wallaby.fa"};
vector<int> sequencCount = {17187, 20792, 20511, 20295, 18222, 16182, 20792, 20500, 16558, 18524, 20837, 20713, 20641, 20338, 20669, 20627, 20357, 17848, 20863, 20786, 20095, 19391, 19440, 18605, 20842, 17222, 18895, 20714, 20636, 20281, 19052, 17004, 15472, 18440, 18855, 18906};*/

/*7 primates*/
/*vector<string> speciesName = {"Baboon", "Common_ch", "Gibbon", "Gorrila", "Human", "Orangutan", "Pigmy_ch"};
vector<string> file1Name = {"Baboon_counts_dumps.fa", "Common_ch_counts_dumps.fa", "Gibbon_counts_dumps.fa", "Gorrila_counts_dumps.fa", "Human_counts_dumps.fa", "Orangutan_counts_dumps.fa", "Pigmy_ch_counts_dumps.fa"};;
vector<int> sequencCount = {16757, 16800, 16708, 16598, 16806, 16624, 16791};*/

//ChunkMiss
vector<string> speciesName = {"chimp", "gorrila", "human"};
vector<string> file1Name = {"chimp_counts_dumps.fa", "gorrila_counts_dumps.fa", "human_counts_dumps.fa"};
vector<int> sequencCount = {1648526, 1843478, 1906376};

void printMatrix(vector<Species> species, double **distanceMatrixFactorcount);
int main() {
    double maxEntropyKmer = 0;
    double max = 0;
    double* entropy = new double[65];
    for (int iter = 15; iter < 16; ++iter) {
        vector<Species> species;
        for (int i = 0; i < speciesName.size() ; i++) {
            Species s(speciesName[i], iter);
            s.init("/home/yusuf/Ruby/files/ChunkMiss/Data" + std::to_string(iter) + "/" + file1Name[i], sequencCount[i],
                   false);
            species.push_back(s);
        }
        entropy[iter] = Species::calculateEntropy(species);
        if(entropy[iter] > max) {
            max = entropy[iter];
            maxEntropyKmer = iter;
        }
    }
    for (int iter = 15; iter < 16; ++iter) {
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
            Species s(speciesName[i], iter);
            s.init("/home/yusuf/Ruby/files/ChunkMiss/Data" + std::to_string(iter) + "/" + file1Name[i], sequencCount[i], true);
            species.push_back(s);
        }
        for (int j = 0; j < speciesName.size(); j++) {
            for (int i = 0; i < speciesName.size(); i++) {
                if(j == i) {
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

        //printMatrix(species, distanceMatrixMahalnobis);
        //Construct tree and find distance FactorCount
        dist_matrix *dmatE = load_file(species, distanceMatrixMahalnobis);
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
        string newickFormat = btree_print_newick_tree(phyl_tree);
        printf("\n\n");
        btree_print_tree(phyl_tree);
        string s1(newickFormat);
        //cout<< s1<<endl;
        //string s2("((guinea_pig: 0.391644,(armadillo: 0.376809,(chicken: 0.580288,heagehog:0.383956):0.013781):0.006984):0.005294,((elephant: 0.336681,(((cat: 0.106369,clouded_leopard:0.118604):0.144929,(dog: 0.266843,ferret:0.252705):0.018200):0.055911,(cow: 0.321540,horse:0.295120):0.004875):0.015553):0.005461,(galago: 0.341609,(dusky_titi: 0.212400,((baboon: 0.128577,colobus_monkey:0.132463):0.026111,(gibbon: 0.137576,((gorilla: 0.111531,(chimp: 0.112824,human:0.100773):0.006357):0.014026,orangatun:0.136581):0.010076):0.019880):0.033993):0.103275):0.018019):0.005294);");
        //ecoli
        //string s2("((((D1Sd197:1.0,(CB9615:1.0,(Sakai:1.0,EDL933:1.0):1.0):1.0):1.0,(((((SSSs046:1.0,(B18BS512:1.0,B4Sb227:1.0):1.0):1.0,(E24377A:1.0,(IAI1:1.0,SE11:1.0):1.0):1.0):1.0):1.0,(F5b8401:1.0,(F2a2457T:1.0,F2a301:1.0):1.0):1.0):1.0,((ATCC8739:1.0,HS:1.0):1.0,((BW2952:1.0,DH10B:1.0):1.0,(MG1655:1.0,W3110:1.0):1.0):1.0):1.0):1.0):1.0,UMN026:1.0):1.0,((IAI39:1.0,SMS35:1.0):1.0,(E234869:1.0,(536a:1.0,((S88:1.0,(APEC01:1.0,UTI89:1.0):1.0):1.0,(ED1a:1.0,CFT073:1.0):1.0):1.0):1.0):1.0):1.0);");
        //fish
        //string s2("((((((NC_013564:1.0000,NC_013577:1.0000):100.0000,((NC_009459:1.0000,NC_009066:1.0000):1.0000,(NC_010205:1.0000,(NC_012055:1.0000,NC_009067:1.0000):1.0000):1.0000):1.0000):1.0000,((NC_009060:1.0000,NC_009059:1.0000):1.0000,(NC_009064:1.0000,NC_009065:1.0000):1.0000):1.0000):1.0000,(NC_011179:1.0000,NC_011177:1.0000):1.0000):1.0000,(NC_011170:1.0000,NC_011169:1.0000):1.0000):1.0000,(NC_011168:1.0000,NC_009058:1.0000):1.0000,(NC_011171:1.0000,((NC_009057:1.0000,(NC_013750:1.0000,NC_013663:1.0000):1.0000):1.0000,(NC_009062:1.0000,(NC_018814:1.0000,(NC_018815:1.0000,NC_009063:1.0000):1.0000):1.0000):1.0000):1.0000):1.0000);");
        //7 primates
        //string s2("(Baboon:0.0000,(Gibbon:0.0000,(Orangutan:0.0000,(Gorrila:0.0000,(Human:0.0000,(Pigmy_ch:0.0000,Common_ch:0.0000):0.0000):0.0000):0.0000):0.0000):0.0000);");
        //44 species
        //string s2("((((((gibbon:0.0000,(orangatun:0.0000,(gorilla:0.0000,(human:0.0000,chimp:0.0000):0.0000):0.0000):0.0000):0.0000,(colobus_monkey:0.0000,(vervet:0.0000,(baboon:0.0000,macaque:0.0000):0.0000):0.0000):0.0000):0.0000,(dusky_titi:0.0000,(owl_monkey:0.0000,(squirrel_monkey:0.0000,marmoset:0.0000):0.0000):0.0000):0.0000):0.0000,(galago:0.0000,mouse_lemur:0.0000):0.0000):0.0000,(rabbit:0.0000,(guinea_pig:0.0000,(rat:0.0000,mouse:0.0000):0.0000):0.0000):0.0000):0.0000,(((horse:0.0000,(pig:0.0000,(muntjak_indian:0.0000,(sheep:0.0000,cow:0.0000):0.0000):0.0000):0.0000):0.0000,(armadillo:0.0000,(heagehog:0.0000,shrew:0.0000):0.0000):0.0000):0.0000,((cat:0.0000,clouded_leopard:0.0000):0.0000,(ferret:0.0000,dog:0.0000):0.0000):0.0000):0.0000,((wallaby:0.0000,(elephant:0.0000,tenrec:0.0000):0.0000):0.0000,(platypus:0.5000,chicken:0.5000):0.4999):0.0000);");
        //ChunkMiss
        string s2 ("(gorrila:0.0000, (human:0.0000,chimp:0.0000):0.0000);");
        PhyloTree t1 = PhyloTree(s1, false);
        PhyloTree t2 = PhyloTree(s2, false);
        double  dd, d;
        dd = Distance::getRobinsonFouldsDistance(t1, t2, false);
        //d = Distance::getWeightedRobinsonFouldsDistance(t1, t2, false);
        cout << "Entropy " << entropy[iter];
        cout<<"RobinsonFouldDistance Factorcount " << "for kemr " << iter << " is " << dd<<endl;
        //cout<<d<<endl;

        //btree_print_tree(phyl_tree);
       // printMatrix(species, distanceMatrixFactorcount);
        delete[] distanceMatrixFactorcount;
        dist_matrix_free(dmatE);
        btree_storage_free(tree_storage);
        species.clear();
    }


    return 0;
}

void printMatrix(vector<Species> species, double **distanceMatrixFactorcount) {

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