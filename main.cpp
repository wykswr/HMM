// @author Yukai Wang, student#: 34271296
#include <iostream>
#include "Viterbi.h"
#include <fstream>
#include <tuple>
#include <sstream>

using namespace std;

tuple<vector<string>, vector<string>> read_file(const string &path) {
    ifstream my_file;
    my_file.open(path);
    string line;
    bool indicator = false;
    vector<string> s1, s2;
    int cnt = 0;
    if (my_file.is_open()) {
        while (getline(my_file, line)) {
            if (indicator) {
                indicator = false;
                if (cnt++ % 2 == 0)
                    s1.push_back(line);
                else
                    s2.push_back(line);
            }
            indicator = line.rfind('>', 0) == 0 || line.rfind('@', 0) == 0;
        }
        my_file.close();
    } else
        throw runtime_error("Could not open the file!");
    return make_tuple(s1, s2);
}

tuple<vector<string>, vector<string>> align(const tuple<vector<string>, vector<string>> &sequences) {
    auto &seq1 = get<0>(sequences);
    auto &seq2 = get<1>(sequences);
    vector<string> states, alignments;
    if (seq1.size() != seq2.size())
        throw runtime_error("Pairs' number doesn't match!");
    for (int i = 0; i < seq1.size(); ++i) {
        auto hmm = Viterbi(seq1[i], seq2[i]);
        auto state = hmm.hidden_chain();
        stringstream st;
        copy(state.begin(), state.end(), std::ostream_iterator<int>(st, ""));
        auto alignment = hmm.alignment();
        states.push_back(st.str());
        alignments.push_back(get<0>(alignment));
        alignments.push_back(get<1>(alignment));
    }
    return make_tuple(states, alignments);
}

int main(int argc, char **argv) {
    string fastq, chain_out, alignment_out;
    vector<string> args(argv + 1, argv + argc);
    string usage = "Syntax: HMM -i <infile> -a <path output> -b <alignment output>";
    if (argc == 1) {
        cout << usage << endl;
        return 0;
    }
    for (auto p = args.cbegin(); p != args.cend(); p++) {
        if (*p == "-h" || *p == "--help") {
            cout << usage << endl;
            return 0;
        }
        if (*p == "-i") {
            fastq = *++p;
            continue;
        }
        if (*p == "-a") {
            chain_out = *++p;
            continue;
        }
        if (*p == "-b") {
            alignment_out = *++p;
            continue;
        }
    }
    ofstream chains, alignments;
    chains.open(chain_out);
    alignments.open(alignment_out);
    int cnt = 0;
    auto result = align(read_file(fastq));
    for (auto &r: get<0>(result))
        chains << 0 << r << 4 << endl;
    for (auto &r: get<1>(result)) {
        alignments << r << endl;
        if (cnt++ % 2 == 1)
            alignments << endl;
    }
    chains.close();
    alignments.close();
    return 0;
}
