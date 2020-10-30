//
// Created by Yusuf on 9/11/2019.
//
#include "../headers/species.h"

map<string, int> Species::totalCount ;

Species::Species(string n, int kmer) {
    this->name = n;
    this->kmer = kmer;
}

int Species::init(string fileName, int sCount, bool calculateSd) {
    this->sequenceLength = sCount;
    string line;
    string line2;
    ifstream myfile (fileName);
    /*ifstream myfile2 (fileName2);*/
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
            if(Species::totalCount.count(line2) < 1) {
                Species::totalCount[line2] = x;
            } else {
                Species::totalCount[line2] += x;
            }
            this->kmers.push_back(line2);
            if(calculateSd) {
                this->calculateStandardDeviation(line2);
            }
        }
        this->averageKmerCounts = sum / (double) count;

        myfile.close();
    } else {
        cout<< "cannot find file" << fileName;
        exit(0);
    }
    /*if (myfile2.is_open()) {
        this->sequenceLength = 0;
        while ( getline (myfile2,line2) ) {
            for (int i = 0; i < line2.length(); ++i) {
                if(line2[i] == 'A' || line2[i] == 'C' || line2[i] == 'G' || line2[i] == 'T' || line2[i] == 'N') {
                    this->sequenceLength++;
                }
            }
        }
    }*/

    return 0;
}

vector<double> Species::calculateDistance(Species s) {
    vector<double> v;
    double sumEuclidian = 0;
    double temp = (double)(min(this->sequenceLength, s.sequenceLength) - this->kmer +1);
    double milhanobis = 0;
    double sumFactor = 0;
    for (int i = 0; i < this->kmers.size(); ++i) {
        sumEuclidian += pow((this->kmerCount[this->kmers[i]] - s.kmerCount[this->kmers[i]]), 2);
        sumFactor += (min(this->kmerCount[this->kmers[i]], s.kmerCount[this->kmers[i]])) / temp;
        if(s.standardDeviation.count(this->kmers[i]) > 0) {
            milhanobis  += pow(( (this->kmerCount[this->kmers[i]] / this->standardDeviation[this->kmers[i]] ) - (s.kmerCount[this->kmers[i]] / s.standardDeviation[this->kmers[i]])), 2);
        } else {
            milhanobis  += pow((this->kmerCount[this->kmers[i]] / this->standardDeviation[this->kmers[i]] ), 2);
        }
        //cout<<milhanobis<<endl;
    }
    for (int i = 0; i < s.kmers.size(); i++) {
        if(this->kmerCount.count(s.kmers[i]) < 1) {
            sumEuclidian += pow((s.kmerCount[s.kmers[i]]), 2);
            milhanobis += pow((s.kmerCount[s.kmers[i]] /  s.standardDeviation[s.kmers[i]]), 2);
        }
    }

    v.push_back(sqrt(sumEuclidian));
    v.push_back(milhanobis);
    v.push_back(abs(log10(0.1 + sumFactor)));
    return v;
}

double Species::calculateEntropy(vector<Species> species) {
    double entropy = 0.0;
    for (int i = 0; i < species.size(); ++i) {
        Species s = species[i];
        for (int j = 0; j < s.kmers.size(); ++j) {
            string kmer = s.kmers[j];
            double p = (double) s.kmerCount[kmer]/ Species::totalCount[kmer];
            entropy += (-p  *  log10(p));
        }
    }
    return entropy;
}

void Species::calculateStandardDeviation(string subsequence) {
    int n;
    int L = this->kmer;
    int M = this->sequenceLength;
    n = M - L + 1;
    int i,j=0,k, flag;
    double p = pow (0.25, L), sum=0.0;
    int Q[64];
    double var_x, expec_x;
    for (i=0; subsequence[i]!='\0'; i++){
        for(k=0; k<=i; k++){
            if(subsequence[k] == subsequence[k + L -1-i]){
                flag=1;
            }
            else{
                flag=0;
                break;
            }
        }
        if(flag==1){
            Q[i]=1;
        }
        else {
            Q[i]=0;
        }
    }
    int minimum= min(L-1, n-1);
    for (k=0; k< minimum;k++){
        sum += (n-(k+1)) * Q[L-(k+1)] * (pow(0.25, k+1));
    }
    expec_x= n* p;
    var_x= (expec_x*(1-expec_x)) + (pow(p,2) * (n-L) * (n-L +1)) + (2 *p*sum);
    double sd = sqrt(var_x);
    this->standardDeviation[subsequence] = sd;
}
