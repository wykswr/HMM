//
// Created by Yukai Wang on 2022-11-07.
// Student#: 34271296

#include "AlignmentModel.h"
#include <algorithm>


AlignmentModel::AlignmentModel(const std::string &seq1, const std::string &seq2)
        : match(std::vector<std::vector<double>>(seq1.length() + 1)),
          x_only(std::vector<std::vector<double>>(seq1.length() + 1)),
          y_only(std::vector<std::vector<double>>(seq1.length() + 1)),
          seq1(seq1), seq2(seq2), trans(4),
          inferred(false) {
    for (int i = 0; i < seq1.length() + 1; ++i) {
        match[i] = std::vector<double>(seq2.length() + 1, -1);
        x_only[i] = std::vector<double>(seq2.length() + 1, -1);
        y_only[i] = std::vector<double>(seq2.length() + 1, -1);
    }
    match[0][0] = 1;
    x_only[0][0] = 0;
    y_only[0][0] = 0;
    trans[0] = std::vector<double>{0, 1. / 3, 1. / 3, 1. / 3, 0};
    trans[1] = std::vector<double>{0, 7. / 10, 1. / 10, 1. / 10, 1. / 10};
    trans[2] = std::vector<double>{0, 1. / 4, 13. / 20, 0, 1. / 10};
    trans[3] = std::vector<double>{0, 1. / 4, 0, 13. / 20, 1. / 10};
}

double AlignmentModel::prob(int i, int j) const {
    if (i <= 0 || j <= 0)
        return 1. / 4;
    auto a = seq1[i - 1];
    auto b = seq2[j - 1];
    if (a == b)
        return 4. / 28;
    return 1. / 28;
}

double AlignmentModel::score(const std::vector<std::vector<double>> &arr, int i, int j) {
    if (i == -1 || j == -1)
        return 0;
    return arr[i][j];
}

void AlignmentModel::infer() {
    for (int i = 0; i < seq1.length() + 1; ++i) {
        for (int j = 0; j < seq2.length() + 1; ++j) {
            match[i][j] = match[i][j] != -1 ? match[i][j] : prob(i, j) *
                                                            std::max({trans[1][1] * score(match, i - 1, j - 1),
                                                                      trans[2][1] * score(x_only, i - 1, j - 1),
                                                                      trans[3][1] * score(y_only, i - 1, j - 1)});
            x_only[i][j] =
                    x_only[i][j] != -1 ? x_only[i][j] : prob(i, -1) * std::max({trans[1][2] * score(match, i - 1, j),
                                                                                trans[2][2] * score(x_only, i - 1, j)});
            y_only[i][j] =
                    y_only[i][j] != -1 ? y_only[i][j] : prob(-1, j) * std::max({trans[1][3] * score(match, i, j - 1),
                                                                                trans[3][3] * score(y_only, i, j - 1)});
        }
    }
    inferred = true;
}

void AlignmentModel::trace_back(int i, int j, int state, std::vector<int> &record) const {
    if (i < 0 || j < 0)
        return;
    if (i == 0 && j == 0)
        return;
    record.push_back(state);
    double next_match, next_x, next_y;
    int next_state;
    switch (state) {
        case 1:
            next_match = trans[1][1] * score(match, i - 1, j - 1);
            next_x = trans[2][1] * score(x_only, i - 1, j - 1);
            next_y = trans[3][1] * score(y_only, i - 1, j - 1);
            if (next_match > next_x)
                next_state = next_match > next_y ? 1 : 3;
            else
                next_state = next_x > next_y ? 2 : 3;
            trace_back(i - 1, j - 1, next_state, record);
            break;
        case 2:
            next_match = trans[1][2] * score(match, i - 1, j);
            next_x = trans[2][2] * score(x_only, i - 1, j);
            next_state = next_match > next_x ? 1 : 2;
            trace_back(i - 1, j, next_state, record);
            break;
        case 3:
            next_match = trans[1][3] * score(match, i, j - 1);
            next_y = trans[3][3] * score(y_only, i, j - 1);
            next_state = next_match > next_y ? 1 : 3;
            trace_back(i, j - 1, next_state, record);
            break;
    }
}

std::vector<int> AlignmentModel::r_hidden_chain() {
    if (!inferred)
        infer();
    std::vector<int> chain;
    int state, n, m;
    n = static_cast<int>(seq1.length());
    m = static_cast<int>(seq2.length());
    if (score(match, n, m) > score(x_only, n, m))
        state = score(match, n, m) > score(y_only, n, m) ? 1 : 3;
    else
        state = score(x_only, n, m) > score(y_only, n, m) ? 2 : 3;
    trace_back(n, m, state, chain);
    return chain;
}

std::vector<int> AlignmentModel::hidden_chain() {
    auto rhc = r_hidden_chain();
    return {rhc.rbegin(), rhc.rend()};
}

std::tuple<std::string, std::string> AlignmentModel::alignment() {
    auto chain = hidden_chain();
    auto p1 = seq1.cbegin();
    auto p2 = seq2.cbegin();
    std::vector<char> s1(hidden_chain().size());
    std::vector<char> s2(hidden_chain().size());
    int idx = 0;
    for (auto r: chain) {
        switch (r) {
            case 1:
                s1[idx] = *p1++;
                s2[idx] = *p2++;
                break;
            case 2:
                s1[idx] = *p1++;
                s2[idx] = '-';
                break;
            case 3:
                s1[idx] = '-';
                s2[idx] = *p2++;
                break;
        }
        idx++;
    }
    return std::make_tuple(std::string(s1.begin(), s1.end()), std::string(s2.begin(), s2.end()));
}
