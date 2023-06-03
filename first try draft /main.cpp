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
    std::string Aprime = "AGGACCGCACAACCTTGCAGCTCAGCGACTCGTGGGGTCACACACACTGG";
    std::string Bprime = "GATCACACATAACTGTGGTGTCATGCATTTGGTATCTTTTAATTTTTAGG";

    int m = A.size();
    int n = B.size();

    // Dynamic programming tables
    std::vector<std::vector<double>> T1(m+1, std::vector<double>(n+1, 0));
    std::vector<std::vector<double>> T2(m+1, std::vector<double>(n+1, 0));
    std::vector<std::vector<double>> T3(m+1, std::vector<double>(n+1, 0));


    num_threads = 10; // p is the number of processors 
    std::vector<std::thread> processors(num_threads);


    /*We define special columns of the dynamic programing tables T1, T2, T3 
        to be columns (0, n/num_threads, 2n/num_threads, ..., n)
        and rows (0, m/num_threads, 2m/num_threads, ..., m)
    */

    std::vector<double> special_columns = {0, n/num_threads, 2*n/num_threads, ..., n};
    std::vector<double> special_rows = {0, m/num_threads, 2*m/num_threads, ..., m};
    
    for (int k=0; k<num_threads-1; k++){

        // we want to give each thread substrings from Aprime to Bprime 
        std::string subA = Aprime.substr(k*m/num_threads +1,(k+1)*m/num_threads);
        std::string subB = Bprime.substr(k*n/num_threads +1,(k+1)*n/num_threads);


       processors[k] = 



    }
       

    return 0;
    
}