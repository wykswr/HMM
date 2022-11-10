//
// Created by Yukai Wang on 2022-11-07.
// Student#: 34271296

#ifndef HMM_VITERBI_H
#define HMM_VITERBI_H

#include <string>
#include <vector>
#include <tuple>


class Viterbi {
private:
    std::vector<std::vector<double>> trans;
    const std::string seq1, seq2;
    std::vector<std::vector<double>> match, x_only, y_only;
    bool inferred;

    double prob(int, int) const;

    double score(const std::vector<std::vector<double>> &, int, int) const;

    void infer();

    void trace_back(int, int, int, std::vector<int> &) const;

    std::vector<int> r_hidden_chain();

public:
    Viterbi(const std::string &, const std::string &);

    std::vector<int> hidden_chain();

    std::tuple<std::string, std::string> alignment();
};

#endif //HMM_VITERBI_H
