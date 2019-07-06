#ifndef __SolutionMerger_H__
#define __SolutionMerger_H__

#include <math.h>
#include <map>
#include <vector>
#include <set>
#include <queue>
#include "SettingsHandler.h"
#include <set_obj.h>
#include <merge_sol.h>
#include <STPGModel.h>
#include "QOL/CpuTimer.h"

class SolutionMerger{
  private:
    // int temp;
    SettingsHandler sh;
    // STPGModel *probModel;

  public:
    // SolutionMerger(const SettingsHandler sh, const STPGModel *probModel) : sh(sh),probModel(probModel) {std::cout << "solution merger initialised" << std::endl;};
    SolutionMerger(const SettingsHandler sh){std::cout << "solution merger initialised" << std::endl;};

    void merge(const std::vector<set_obj>&sols, merge_sol &merged_sols);
    void test(){std::cout << "test" << std::endl;};

    ~SolutionMerger(){};

};

#endif