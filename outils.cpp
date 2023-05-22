#include <iostream>
#include <fstream>
#include <set>
#include "outils.h"
#include <string>
#include <vector>
#include <utility> // for std::pair
#include <cmath>

//CONSTANTE H: penalty of the first gap in the alignment
const int H = -2;

int standard_score(char a, char b) {
    std::set<char> sigma = {'A', 'C', 'G', 'T'};
    if (sigma.find(a) == sigma.end() || sigma.find(b) == sigma.end()) {
        std::cerr << "Error: invalid character in input" << std::endl;
        exit(1);
    }
    if (a == b) return 2;
    else return 0;
}

int gap_penalty(int k) {
    //  not the same penalty as the one introfuced in the paper 
    return -k;
}

int sc(std::vector<std::pair<int, int>> alignment) {
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

int esc_s(std::vector<std::pair<int, int>> alignment, double start_type){
    // if start_type equal == -2 or -3 and the first element of the alignment is of type abs(start_type)
    // then we return Constant H 
    // else if strat_type >0 and type and the first element of the alignment is not of type abs(start_type)
    // then we return -infinity 
    // else we return 0
    if (start_type == -2 || start_type == -3){
        if (type(alignment[0]) == abs(start_type)){
            return H;
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

int esc_e(std::vector<std::pair<int, int>> alignment, double end_type){
    // if end_type equal == -2 or -3 and the last element of the alignment is of type abs(end_type)
    // then we return Constant H 
    // else if end_type>0 and type and the last element of the alignment is not of type abs(end_type)
    // then we return -infinity 
    // else we return 0
    if (end_type == -2 || end_type == -3){
        if (type(alignment[alignment.size()-1]) == abs(end_type)){
            return H;
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

int extended_score(std::vector<std::pair<int, int>> alignment, double start_type, double end_type) {
    if (start_type == end_type && start_type ==-1){
        return sc(alignment);
    }
    return sc(alignment)+ esc_s(alignment, start_type) + esc_e(alignment, end_type);
}

