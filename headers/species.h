//
// Created by Yusuf on 9/11/2019.
//

#ifndef KMERS_KMERCOUNTS_H
#define KMERS_KMERCOUNTS_H

#endif //KMERS_KMERCOUNTS_H

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>
#include <algorithm>
using namespace std;

class Species {
public:
    string name;
    int kmer;
    double averageKmerCounts;
    int sequenceLength;
    vector<string> kmers;
    vector<int> kmerCount;
    Species(string n, int kmer);
    int init(string fileName, string fileName2);
    int calculateEuclidianDistance(Species s);
    double calculateMahalnobisDistance(Species s);
    double calculateFactorDistance(Species s);
};


