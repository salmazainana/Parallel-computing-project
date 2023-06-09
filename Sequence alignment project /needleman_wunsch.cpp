#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <mutex>

//TO RUN
//g++ -std=c++11 -o nw needleman_wunsch.cpp
// then ./nw

// Setup the penalties
int match_penalty = 1;
int mismatch_penalty = -1;
int gap_penalty = -1;

// Calculate a cell
void calculate_cell(int i, int j, std::vector<std::vector<int>>& M, const std::string& seq1, const std::string& seq2) {
    int score;
    if(seq1[i - 1] == seq2[j - 1])
        score = match_penalty;
    else
        score = mismatch_penalty;

    M[i][j] = std::max({M[i-1][j-1] + score, M[i-1][j] + gap_penalty, M[i][j-1] + gap_penalty});
}

// Calculate the matrix
void calculate_scores(std::vector<std::vector<int>>& M, const std::string& seq1, const std::string& seq2) {
    // CHANGE NUM OF THREADS
    int num_threads = 4;
    int m = seq1.length();
    int n = seq2.length();

    for(int sum = 2; sum <= m + n; ++sum) {
        int z1 = std::max(1, sum - n);
        int z2 = std::min(sum - 1, m);

        std::vector<std::thread> threads;
        for(int z = z1; z <= z2; z += num_threads) {
            for(int t = z; t < std::min(z + num_threads, z2 + 1); ++t) {
                threads.push_back(std::thread(calculate_cell, t, sum - t, std::ref(M), std::ref(seq1), std::ref(seq2)));
            }

            for(auto& thread : threads) {
                thread.join();
            }

            threads.clear();
        }
    }
}

// Prints the Aligned sequences
void output(const std::vector<std::vector<int>>& M, const std::string& seq1, const std::string& seq2) {
    int i = seq1.length();
    int j = seq2.length();

    std::string alignmentA, alignmentB;

    while(i > 0 && j > 0) {
        int score = M[i][j];
        int score_diag = M[i-1][j-1];
        int score_up = M[i][j-1];
        int score_left = M[i-1][j];

        if(score == score_diag + ((seq1[i-1] == seq2[j-1]) ? match_penalty : mismatch_penalty)) {
            alignmentA += seq1[i-1];
            alignmentB += seq2[j-1];
            --i;
            --j;
        } else if(score == score_up + gap_penalty) {
            alignmentA += '-';
            alignmentB += seq2[j-1];
            --j;
        } else {
            alignmentA += seq1[i-1];
            alignmentB += '-';
            --i;
        }
    }

    while(i > 0) {
        alignmentA += seq1[i-1];
        alignmentB += '-';
        --i;
    }

    while(j > 0) {
        alignmentA += '-';
        alignmentB += seq2[j-1];
        --j;
    }

    std::reverse(alignmentA.begin(), alignmentA.end());
    std::reverse(alignmentB.begin(), alignmentB.end());

    std::cout << "Alignment 1: " << alignmentA << "\n";
    std::cout << "Alignment 2: " << alignmentB << "\n";
}

// Prints the matrix
void table(const std::vector<std::vector<int>>& M) {
    for (const auto& row : M){
        for (int elem : row){
            std::cout << elem << " ";
        }
        std::cout << "\n";
    }
}

void NW(const std::string& seq1, const std::string& seq2) {
    int m = seq1.length();
    int n = seq2.length();
    std::vector<std::vector<int>> M(m+1, std::vector<int>(n+1, 0));

    for(int i = 0; i <= m; i++) M[i][0] = i * gap_penalty;
    for(int j = 0; j <= n; j++) M[0][j] = j * gap_penalty;

    // Calculates the matrix
    calculate_scores(M, seq1, seq2);

    // Outputs table
    table(M);

    // Outputs alligned sequence
    output(M, seq1, seq2);
}

int main() {
    std::string seq1 = "";
    std::string seq2 = "";

    std::cout << "L1: " << seq1.length() << "\n";
    std::cout << "L2: " << seq2.length() << "\n";

    auto start = std::chrono::high_resolution_clock::now(); // Start the clock

    NW(seq1, seq2);

    auto end = std::chrono::high_resolution_clock::now(); // Stop the clock
    std::chrono::duration<double> diff = end-start; // Find the elapsed time
    std::cout << "Execution time: " << diff.count() << " s\n";

    return 0;
}
