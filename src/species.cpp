//
// Created by Yusuf on 9/11/2019.
//
#include "../headers/species.h"

Species::Species(string n, int kmer) {
    this->name = n;
    this->kmer = kmer;
}

int Species::init(string fileName, string fileName2) {
    string line;
    string line2;
    ifstream myfile (fileName);
    ifstream myfile2 (fileName2);
    if (myfile.is_open())
    {
        int sum = 0;
        int count = 0;
        while ( getline (myfile,line) )
        {
            getline (myfile,line2);
            string str = line.substr(1, line.size());
            stringstream geek(str);
            int x = 0;
            geek >> x;
            sum += x;
            count++;
            this->kmerCount.push_back(x);
            this->kmers.push_back(line2);
        }
        this->averageKmerCounts = sum / (double) count;
        myfile.close();
    }
    if (myfile2.is_open()) {
        this->sequenceLength = 0;
        while ( getline (myfile2,line2) ) {
            for (int i = 0; i < line2.length(); ++i) {
                if(line2[i] == 'A' || line2[i] == 'C' || line2[i] == 'G' || line2[i] == 'T' || line2[i] == 'N') {
                    this->sequenceLength++;
                }
            }
        }
    }
    return 0;
}

int Species::calculateEuclidianDistance(Species s) {
    int sum = 0;
    for (int i = 0; i < this->kmers.size(); ++i) {
        for (int j = 0; j <s.kmers.size() ; ++j) {
            if(this->kmers[i] == s.kmers[j]) {
                sum += pow((this->kmerCount[i] - s.kmerCount[j]), 2);
                break;
            }
        }
    }
    return sum;
}

double Species::calculateMahalnobisDistance(Species s) {
    double sum = 0;
    for (int i = 0; i < this->kmers.size(); ++i) {
        for (int j = 0; j <s.kmers.size() ; ++j) {
            if(this->kmers[i] == s.kmers[j]) {
                sum += pow(((this->kmerCount[i]/ abs(this->averageKmerCounts - this->kmerCount[i])) - s.kmerCount[j]/ abs(s.averageKmerCounts - s.kmerCount[j])), 2);
                break;
            }
        }
    }
    return sum;
}

double Species::calculateFactorDistance(Species s) {
    double sum = 0;
    double temp = (double)(min(this->sequenceLength, s.sequenceLength) - this->kmer +1);
    for (int i = 0; i < this->kmers.size(); ++i) {
        for (int j = 0; j <s.kmers.size() ; ++j) {
            if(this->kmers[i] == s.kmers[j]) {
                sum += (min(this->kmerCount[i], s.kmerCount[j])) / temp;
                break;
            }
        }
    }
    double x = abs(log10(0.1 + sum));
    return x;
}
