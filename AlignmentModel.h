//
// Created by Yukai Wang on 2022-11-07.
// Student#: 34271296

#ifndef HMM_ALIGNMENTMODEL_H
#define HMM_ALIGNMENTMODEL_H

#include <string>
#include <vector>
#include <tuple>


class AlignmentModel {
protected:
    // Transition matrix
    std::vector<std::vector<double>> trans;
    const std::string seq1, seq2;
    // Likelihood matrix
    std::vector<std::vector<double>> match, x_only, y_only;
    bool inferred;

    // The P_{x_i, y_j}, if y_j is "-", you can set j to -1.
    double prob(int, int) const;

    // Getter of likelihood matrix, embedding base case.
    static double score(const std::vector<std::vector<double>> &, int, int) ;

    // The Vierbi dynamic programming algorithm.
    void infer();

    // Util for r_hidden_chain.
    void trace_back(int, int, int, std::vector<int> &) const;

    // Trace back the hidden states chain Z_i via likelihood matrix inferred.
    std::vector<int> r_hidden_chain();

public:
    // Initialize transition matrix and likelihood matrix.
    AlignmentModel(const std::string &, const std::string &);

    // Reverse r_hidden_chain
    std::vector<int> hidden_chain();

    // Using r_hidden_chain to generate the alignment.
    std::tuple<std::string, std::string> alignment();
};

#endif //HMM_ALIGNMENTMODEL_H
