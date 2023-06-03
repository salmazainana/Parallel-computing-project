#include <iostream>
#include <fstream>
#include <set>
#include "outils.h"
#include <string>
#include <vector>
#include <utility> // for std::pair
#include <cmath>
#include <algorithm>

//CONSTANTE H: penalty of the first gap in the alignment
const double h = -2;
const double g = -1;

double standard_score(char a, char b) {
    std::set<char> sigma = {'A', 'C', 'G', 'T'};
    if (sigma.find(a) == sigma.end() || sigma.find(b) == sigma.end()) {
        std::cerr << "Error: invalid character in input" << std::endl;
        exit(1);
    }
    if (a == b) return 2;
    else return 0;
}

double gap_penalty(int k) {
    //  not the same penalty as the one introfuced in the paper 
    return -k;
}

double sc(std::vector<std::pair<int,int>> alignment) {
    int score = 0; 
    for (auto it = alignment.begin(); it != alignment.end(); ++it) {
        if (it->first == 0 || it->second == 0) {
            score += gap_penalty(1);
        }
        else{
            score += standard_score(it->first, it->second);
        }
    }
    return score;
}

int type(std::pair<int, int> element){
    // return the type of the element(i,j)
    // if i = 0 return 2
    // if j = 0 return 3
    // else 0 
    if (element.first == 0 ){
        return 2;
    }
    else if (element.second == 0){
        return 3;
    }
    else{
        return 1;
    }
}

double esc_s(std::vector<std::pair<int, int>> alignment, double start_type){
    // if start_type equal == -2 or -3 and the first element of the alignment is of type abs(start_type)
    // then we return Constant H 
    // else if strat_type >0 and type and the first element of the alignment is not of type abs(start_type)
    // then we return -infinity 
    // else we return 0
    if (start_type == -2 || start_type == -3){
        if (type(alignment[0]) == abs(start_type)){
            return h;
        }
        else{
            return 0;
        }
    }
    if (start_type > 0){
        if (type(alignment[0]) != abs(start_type)){
            return -INFINITY;
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }

}

double esc_e(std::vector<std::pair<int, int>> alignment, double end_type){
    // if end_type equal == -2 or -3 and the last element of the alignment is of type abs(end_type)
    // then we return Constant h 
    // else if end_type>0 and type and the last element of the alignment is not of type abs(end_type)
    // then we return -infinity 
    // else we return 0
    if (end_type == -2 || end_type == -3){
        if (type(alignment[alignment.size()-1]) == abs(end_type)){
            return h;
        }
        else{
            return 0;
        }
    }
    if (end_type > 0){
        if (type(alignment[alignment.size()-1]) != abs(end_type)){
            return -INFINITY;
        }
        else{
            return 0;
        }
    }
    else{
        return 0;
    }
}

double extended_score(std::vector<std::pair<int, int>> alignment, double start_type, double end_type) {
    if (start_type == end_type && start_type ==-1){
        return sc(alignment);
    }
    return sc(alignment)+ esc_s(alignment, start_type) + esc_e(alignment, end_type);
}

void initializeTables(std::vector<std::vector<double>>& T1, std::vector<std::vector<double>>& T2, std::vector<std::vector<double>>& T3,
                      int iprime, int jprime, int isec, int jsec, double start_type) {
    // Initialize tables with appropriate sizes
    T1[iprime - 1][jprime - 1] = -INFINITY;
    T2[iprime - 1][jprime - 1] = -INFINITY;
    T3[iprime - 1][jprime - 1] = -INFINITY;

    //leftmost column 
    for (int i = iprime; i < isec + 1; i++) {
        T1[i][jprime - 1] = -INFINITY;
        T2[i][jprime - 1] = -INFINITY;
        if (start_type == -3){
            T3[i][jprime - 1] = -g*(i-iprime+1) ; 
        }
        else if ( (start_type== 3) || (start_type== -1) || (start_type== -2) ) {
            T3[i][jprime - 1] = -h -g*(i - iprime+1);
        }
        else{    
            T3[i][jprime - 1] =-INFINITY;
        }
    }
    //topmost column 
    for (int j = jprime; j < jsec + 1; j++) {
        T1[iprime - 1][j] = -INFINITY;
        T3[iprime - 1][j] = -INFINITY;
        if (start_type == -2){
            T2[iprime -1 ][j] = - g*(j-jprime+1) ; 
        }
        else if ( (start_type== -3) || (start_type== -1) || (start_type== 2) ){
            T2[iprime -1 ][j] = - h - g*(j- jprime + 1);
        }
        else{    
            T2[iprime -1 ][j] =-INFINITY;
        }
    }

    // Set initial values in t1, t2, and t3
    if (start_type <=1){
        //start_type in {-3,-2,-1,1,2,3}
        if (start_type == 1 || start_type == -1){
            T1[iprime - 1][jprime - 1] = 0;
        }
        if (start_type == -2){
            T2[iprime - 1][jprime - 1] = 0;
        }
        if (start_type == -3){
            T3[iprime - 1][jprime - 1] = 0;
        }
    }
    
    if (start_type > 0) {
        if (start_type == 2) {
            T2[iprime - 1][jprime] = -(g + h);
        }
        else{
            T2[iprime - 1][jprime] =-INFINITY;
        }
        if (start_type == 3) {
            T3[iprime][jprime - 1] = -(g + h); 
        }
        else{
            T3[iprime][jprime - 1] =-INFINITY;
        }
    }

    // Fill the remaining cells of t1, t2, and t3
    for (int i = iprime; i < isec+1; i++) {
        for (int j = jprime; j < jsec +1; j++) {
            T1[i][j] =  sc({{i, j}}) + std::max({T1[i - 1][j - 1], T2[i - 1][j - 1], T3[i - 1][j - 1]});
            T2[i][j] = std::max({T1[i - 1][j] - (g + h), T2[i - 1][j] - g, T3[i - 1][j] - (g + h)});
            T3[i][j] = std::max({T1[i - 1][j] - (g + h), T2[i-1][j] - (g + h), T3[i-1][j] - g});
        }

    }
}

double hprime(int k, double end_type){
    if (end_type == -2 || end_type == -3){
        if (k == end_type){
            return h;
        }
        else{
            return 0;
        }
    }
    return 0;
}


// Optimal score function for the alignment of two subsequences
// if end_type > 0 then T|end_type|[isec,jsec] is the optimal score
// if end_type < 0 then max(T1[isec,jsec],T2[isec,jsec]+hprime(-2), T3[isec,jsec]+hprime(-3))is the optimal score
// such that if k==end_type and end_type in {-2,-3} then hprime(k)=h else hprime(k) == 0

double optimal_score(std::vector<std::vector<double>>& T1, std::vector<std::vector<double>>& T2, std::vector<std::vector<double>>& T3,
                  int iprime, int jprime, int isec, int jsec, double end_type) {
    if (end_type > 0) {
        if (end_type == 1) {
            return T1[isec][jsec];
        }
        else if (end_type == 2) {
            return T2[isec][jsec];
        }
        else {
            return T3[isec][jsec];
        }
    }
    else {
        double temp = T2[isec][jsec]+ hprime(-2, end_type);
        double temp1 = std::max( T1[isec][jsec],temp);
        double temp2 = T3[isec][jsec] + hprime(-3, end_type) ;
        return std::max( temp1,temp2);
    }
}


//Traceback function for the alignment of two subsequences
std::vector<std::pair<int, int>> traceback( std::vector<std::vector<double>>& T1, std::vector<std::vector<double>>& T2, std::vector<std::vector<double>>& T3,
                  int iprime, int jprime, int isec, int jsec, double start_type){
    std::vector<std::pair<int, int>> alignment;
    std::pair<int, int> s = {iprime, jprime};
    alignment.push_back(s);
    if (end_type>0){
        
    }



}

/*  We want to find elements of a partial balanced partition 
    which is a list of cells where the subdivision occurs 
    by using the tables T1 T2 T3. 
*/

// if [i,j]k is on the path of the solution C through T
// then the original problem of finding alignment between A and B with start_type = s, end_type = e
// is reduced to finding alignment between A[1,i] and B[1,j], start_type = s, end_type = k 
// and A[i+1,m] and B[j+1,n], start_type = -k, end_type = e

// To find cells that lie on the solution we use the method in paper [10]

/*Reverse method:
 We consider revA and revB as the reverse of A and B respectively.
    We use the same tables TR1 TR2 TR3 to find the cells that lie on the solution.
    We use the same method as above but with the following changes:
    1. start_type = e and end_type = k
    2. We start from [isec,jsec] and end at [iprime,jprime]
    The optimal score is found in [iprime,jprime] and the cells that lie on the solution
     are found by backtracking from [iprime,jprime] to [isec,jsec]
*/ 

// This function finds the cells that lie on the solution


void ReversedTables(std::vector<std::vector<double>>& TR1, std::vector<std::vector<double>>& TR2,
                     std::vector<std::vector<double>>& TR3,
                      int isec, int jsec, int iprime, int jprime, double end_type) {
    // NEED TO REDO 
    initializeTables(TR1, TR2, TR3, isec, jsec, iprime, jprime, end_type);}


double opt(int i, int j,std::vector<std::vector<double>>&  T1,std::vector<std::vector<double>>&  T2, 
        std::vector<std::vector<double>>& T3,std::vector<std::vector<double>>&  TR1, std::vector<std::vector<double>>& TR2,
        std::vector<std::vector<double>>& TR3 ){
            // This function finds the optimal score of the alignment between A and B
            // that passes through the cell [i,j]

            double Tmax = std::max(std::max(T1[i][j], T2[i][j]), T3[i][j]);
            double TRmax = std::max(std::max(TR1[i][j], TR2[i][j]), TR3[i][j]);
            double temp = std::max(Tmax+TRmax, T2[i][j] + TR2[i][j] + h);
            return std::max(temp, T3[i][j] + TR3[i][j] + h);
        }

// A solution passes through the cell [i,j] if and only if opt(i,j) is equal to the score of an optimal alighment 
// between A and B ie opt(i,j) == T1[m][n] where m is the length of A and n is the length of B