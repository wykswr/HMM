//
// Created by Yukai Wang on 2022-11-14.
// Student#: 34271296

#include "AlignmentModelProb.h"

double AlignmentModelProb::forward() {
    if (inferred)
        return probability;
    for (int i = 0; i < seq1.length(); ++i) {
        for (int j = 0; j < seq2.length(); ++j) {
            match[i][j] = prob(i, j) *
                          (trans[1][1] * score(match, i - 1, j - 1) + trans[2][1] * score(x_only, i - 1, j - 1) +
                           trans[3][1] * score(y_only, i - 1, j - 1));
            x_only[i][j] = prob(i, -1) * (trans[1][2] * score(match, i - 1, j) + trans[2][2] * score(x_only, i - 1, j));
            y_only[i][j] = prob(-1, j) * (trans[1][3] * score(match, i, j - 1) + trans[3][3] * score(y_only, i, j - 1));
        }
    }
    int n = static_cast<int>(seq1.length() - 1);
    int m = static_cast<int>(seq2.length() - 1);
    probability =
            trans[1][4] * score(match, n, m) + trans[2][4] * score(x_only, n, m) + trans[3][4] * score(y_only, n, m);
    inferred = true;
    return probability;
}

AlignmentModelProb::AlignmentModelProb(const std::string &seq1, const std::string &seq2) : AlignmentModel(seq1, seq2),
                                                                                           probability(-1) {}
