#ifndef __SolutionMerger_H__
#define __SolutionMerger_H__

#include <math.h>
#include <map>
#include <vector>
#include <set>
#include <queue>
#include "SettingsHandler.h"
#include <set_obj.h>
#include <STPGModel.h>
#include "QOL/CpuTimer.h"

class SolutionMerger{
  private:
    SettingsHandler sh;
    STPGModel *probModel;

  public:
    SolutionMerger(const SettingsHandler sh, const STPGModel *probModel) : sh(sh),probModel(probModel) {std::cout << "solution merger initialised" << std::endl;};

    void merge(const std::vector<set_obj>&sols, std::vector<std::vector<int> > &groups, std::vector<int> &group_map);

    ~SolutionMerger(){};

};

#endif