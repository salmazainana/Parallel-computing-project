#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <utility> // for std::pair

double standard_score(char a, char b);
double gap_penalty(int k);
double sc(std::vector<std::pair<int,int>> alignment);
int type (std::pair<int, int> element);
double esc_s(std::vector<std::pair<int,int>> alignment, double start_type);
double esc_e(std::vector<std::pair<int,int>> alignment, double end_type);
int extended_score(std::vector<std::pair<int,int>> alignment);
void initializeTables(std::vector<std::vector<int>>& T1, std::vector<std::vector<int>>& T2, std::vector<std::vector<int>>& T3,
                      int isec, int jsec, int iprime, int jprime, double start_type);
double optimal_score(std::vector<std::vector<double>>& T1, std::vector<std::vector<double>>& T2, std::vector<std::vector<double>>& T3,
                  int isec, int jsec, int iprime, int jprime, double end_type);

