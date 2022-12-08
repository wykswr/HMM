/*
 * @author Yukai Wang
 */
#include <iterator>
#include <iostream>
#include "AlignmentModel.h"
#include "AlignmentModelProb.h"
#include <fstream>
#include <tuple>
#include <sstream>

using namespace std;

tuple<vector<string>, vector<string>> read_file(const string &path) {
    /*
     * Read a paired fastq file into two vectors, if the sequence number is odd, it'll throw an exception
     */
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
    /*
     * Read 2 vectors with same length, do sequence alignment at corresponding location,
     * output the matching result of each 2 sequences to 2 vectors respectively.
     */
    auto &seq1 = get<0>(sequences);
    auto &seq2 = get<1>(sequences);
    if (seq1.size() != seq2.size())
        throw runtime_error("Pairs' number doesn't match!");
    vector<string> states, alignments;
    for (int i = 0; i < seq1.size(); ++i) {
        auto hmm = AlignmentModel(seq1[i], seq2[i]);
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

vector<double> forward(const tuple<vector<string>, vector<string>> &sequences) {
    /*
     *  Read 2 vectors with same length, generate the matching likelihood into a vector.
     */
    auto &seq1 = get<0>(sequences);
    auto &seq2 = get<1>(sequences);
    if (seq1.size() != seq2.size())
        throw runtime_error("Pairs' number doesn't match!");
    vector<double> probabilities(seq1.size());
    for (int i = 0; i < seq1.size(); ++i) {
        auto hmm = AlignmentModelProb(seq1[i], seq2[i]);
        probabilities[i] = hmm.forward();
    }
    return probabilities;
}

int main(int argc, char **argv) {
    string fastq, chain_out, alignment_out;
    bool forward_mode = false;
    // Parser the args
    vector<string> args(argv + 1, argv + argc);
    string usage = "Syntax: HMM [-f] -i <infile> -a <path | probability> -b <alignment>";
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
        if (*p == "-f") {
            forward_mode = true;
            continue;
        }
    }
    // Forward only
    if (forward_mode) {
        ofstream probabilities;
        probabilities.open(chain_out);
        auto result = forward(read_file(fastq));
        for (auto &r: result)
            probabilities << r << endl;
        probabilities.close();
        return 0;
    }
    // Alignment mode
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
