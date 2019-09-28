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
            this->kmerCount[line2] = x;
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

vector<double> Species::calculateDistance(Species s) {
    vector<double> v;
    double sumEuclidian = 0;
    double sumMahalnobis = 0;
    double temp = (double)(min(this->sequenceLength, s.sequenceLength) - this->kmer +1);
    double sumFactor = 0;
    for (int i = 0; i < this->kmers.size(); ++i) {
        sumEuclidian += pow((this->kmerCount[this->kmers[i]] - s.kmerCount[this->kmers[i]]), 2);
        sumMahalnobis += pow(((this->kmerCount[this->kmers[i]]/ abs(this->averageKmerCounts - this->kmerCount[this->kmers[i]])) - s.kmerCount[this->kmers[i]]/ abs(s.averageKmerCounts - s.kmerCount[this->kmers[i]])), 2);
        sumFactor += (min(this->kmerCount[this->kmers[i]], s.kmerCount[this->kmers[i]])) / temp;
    }

    for (int i = 0; i < s.kmers.size(); i++) {
        if(this->kmerCount.count(s.kmers[i]) < 1) {
            sumEuclidian += pow((s.kmerCount[s.kmers[i]]), 2);
            sumMahalnobis += pow((s.kmerCount[s.kmers[i]]/ abs(s.averageKmerCounts - s.kmerCount[s.kmers[i]])), 2);
        }
    }
    v.push_back(sumEuclidian);
    v.push_back(sumMahalnobis);
    v.push_back(abs(log10(0.1 + sumFactor)));
    return v;
}
