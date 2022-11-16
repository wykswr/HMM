//
// Created by Yukai Wang on 2022-11-14.
// Student#: 34271296

#ifndef HMM_ALIGNMENTMODELPROB_H
#define HMM_ALIGNMENTMODELPROB_H

#include "AlignmentModel.h"


class AlignmentModelProb : private AlignmentModel {
private:
    double probability;

public:
    // Using AlignmentModel's likelihood matrix and DP to get the matching likelihood.
    double forward();

    AlignmentModelProb(const std::string &, const std::string &);

};


#endif //HMM_ALIGNMENTMODELPROB_H
