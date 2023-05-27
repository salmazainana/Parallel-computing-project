#include <iostream>
#include <fstream>
#include <set>
#include "outils.h"
#include <string>
#include <vector>
#include <utility> // for std::pair
#include <cmath>

#include "outils.cpp"

int main(){
    
    // Read the input sequences
    std::string A = "AGGACCGCACAACCTTGCAGCTCAGCGACTCGTGGGGTCACACACACTGG";
    std::string B = "GATCACACATAACTGTGGTGTCATGCATTTGGTATCTTTTAATTTTTAGG";

    int m = A.size();
    int n = B.size();
    /*Introduce three dynamic programming tables each one of size (m+1)*(n+1):
        -T1: ai must match bj, ie (i,j) of type 1
        -T2: ai must match a gap, ie (i,j) of type 2
        -T3: bj must match a gap, ie (i,j) of type 3
    */ 
    std::vector<std::vector<double>> T1;
    std::vector<std::vector<double>> T2;
    std::vector<std::vector<double>> T3;


    //Initialize the tables
    double start_type = 1;
    initializeTables(T1, T2, T3, m, n, 1, 1, start_type);

    return 0;
    
}