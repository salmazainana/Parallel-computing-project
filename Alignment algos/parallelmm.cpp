#include <iostream>
#include <vector>
#include <thread>

using namespace std;

// Function to perform sequence alignment for a specific range of rows using the Myers-Miller algorithm
void alignRows(const string& seq1, const string& seq2, vector<vector<int>>& score, int startRow, int endRow) {
    int n = seq2.length();

    for (int i = startRow; i <= endRow; i++) {
        for (int j = 1; j <= n; j++) {
            int match = score[i - 1][j - 1] + (seq1[i - 1] != seq2[j - 1]);
            int deleteOp = score[i - 1][j] + 1;
            int insertOp = score[i][j - 1] + 1;
            score[i][j] = min(match, min(deleteOp, insertOp));
        }
    }
}


// Function to perform sequence alignment using the Myers-Miller algorithm with parallelization
void sequenceAlignment(const string& seq1, const string& seq2) {
    int m = seq1.length();
    int n = seq2.length();

    // Calculate the maximum possible score
    int maxScore = m + n;

    // Create the score matrix
    vector<vector<int>> score(m + 1, vector<int>(n + 1, 0));

    // Fill the first row of the score matrix
    for (int j = 0; j <= n; j++) {
        score[0][j] = j;
    }

    // Fill the first column of the score matrix
    for (int j = 0; j <= n; j++) {
        score[j][0] = j;
    }

    // Define the number of threads
    int numThreads = thread::hardware_concurrency();

    // Calculate the chunk size for each thread
    int chunkSize = (m + numThreads - 1) / numThreads;

    // Create a vector to store the threads
    vector<thread> threads;

    // Perform sequence alignment in parallel
    for (int i = 1; i <= m; i += chunkSize) {
        int startRow = i;
        int endRow = min(i + chunkSize - 1, m);
        threads.emplace_back(alignRows, cref(seq1), cref(seq2), ref(score), startRow, endRow);
    }

    // Wait for all threads to finish
    for (auto& t : threads) {
        t.join();
    }

    // Trace back the alignment
    int i = m;
    int j = n;
    string alignedSeq1 = "";
    string alignedSeq2 = "";

    while (i > 0 && j > 0) {
        int currentScore = score[i][j];
        int diagonalScore = score[i - 1][j - 1];
        int leftScore = score[i][j - 1];
        int upScore = score[i - 1][j];

        if (currentScore == diagonalScore + (seq1[i - 1] != seq2[j - 1])) {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            i--;
            j--;
        } else if (currentScore == upScore + 1) {
            alignedSeq1 = seq1[i - 1] + alignedSeq1;
            alignedSeq2 = '-' + alignedSeq2;
            i--;
        } else {
            alignedSeq1 = '-' + alignedSeq1;
            alignedSeq2 = seq2[j - 1] + alignedSeq2;
            j--;
        }
    }

    // Handle the remaining characters, if any
    while (i > 0) {
        alignedSeq1 = seq1[i - 1] + alignedSeq1;
        alignedSeq2 = '-' + alignedSeq2;
        i--;
    }

    while (j > 0) {
        alignedSeq1 = '-' + alignedSeq1;
        alignedSeq2 = seq2[j - 1] + alignedSeq2;
        j--;
    }

    // Print the alignment
    cout << "Aligned sequence 1: " << alignedSeq1 << endl;
    cout << "Aligned sequence 2: " << alignedSeq2 << endl;
    cout << "Alignment score: " << score[m][n] << "/" << maxScore << endl;

}

int main() {
    string seq1 = "AGTACGCA";
    string seq2 = "TATGC";

    sequenceAlignment(seq1, seq2);

    return 0;
}