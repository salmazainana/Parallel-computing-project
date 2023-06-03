#include <iostream>
#include <vector>
#include <algorithm>

const int GAP_OPEN = -10;      // Gap opening penalty
const int GAP_EXTEND = -1;    // Gap extension penalty
const int MATCH = 2;          // Match score
const int MISMATCH = -1;      // Mismatch score

// Function to perform sequence alignment using Gotoh algorithm
void gotohAlignment(const std::string& seq1, const std::string& seq2) {
    int n = seq1.length();
    int m = seq2.length();

    // Initialize matrices
    std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));
    std::vector<std::vector<int>> gapH(n + 1, std::vector<int>(m + 1, 0));
    std::vector<std::vector<int>> gapV(n + 1, std::vector<int>(m + 1, 0));

    // Initialize gap extension penalties
    for (int i = 1; i <= n; i++) {
        gapH[i][0] = GAP_OPEN + (i - 1) * GAP_EXTEND;
        dp[i][0] = gapH[i][0];
    }
    for (int j = 1; j <= m; j++) {
        gapV[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND;
        dp[0][j] = gapV[0][j];
    }

    // Fill in the matrices
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            // Calculate match/mismatch score
            int matchScore = (seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH;

            // Calculate gap scores
            gapH[i][j] = std::max(gapH[i][j - 1] + GAP_EXTEND, dp[i][j - 1] + GAP_OPEN);
            gapV[i][j] = std::max(gapV[i - 1][j] + GAP_EXTEND, dp[i - 1][j] + GAP_OPEN);

            // Calculate current cell score
            dp[i][j] = std::max({ dp[i - 1][j - 1] + matchScore, gapH[i][j], gapV[i][j] });
        }
    }

    // Traceback to find the alignment
    std::string align1, align2;
    int i = n, j = m;

    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && dp[i][j] == dp[i - 1][j - 1] + ((seq1[i - 1] == seq2[j - 1]) ? MATCH : MISMATCH)) {
            align1 = seq1[i - 1] + align1;
            align2 = seq2[j - 1] + align2;
            i--;
            j--;
        } else if (j > 0 && dp[i][j] == gapH[i][j]) {
            align1 = '-' + align1;
            align2 = seq2[j - 1] + align2;
            j--;
        } else {
            align1 = seq1[i - 1] + align1;
            align2 = '-' + align2;
            i--;
        }
    }

    // Print the alignment
    std::cout << "Alignment:" << std::endl;
    std::cout << align1 << std::endl;
    std::cout << align2 << std::endl;
}

int main() {
    std::string seq1 = "AGTACGCA";
    std::string seq2 = "TATGC";

    gotohAlignment(seq1, seq2);

    return 0;
}