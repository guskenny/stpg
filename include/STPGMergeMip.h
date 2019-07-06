#ifndef __STPGMergeMip_H__
#define __STPGMergeMip_H__

#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
#include "QOL/CpuTimer.h"
#include "merge_sol.h"
#include <math.h>
#include <map>
#include <vector>
#include <set>
#include "SettingsHandler.h"
#include <STPGModel.h>
#include <string>
#include <sstream>
#include <util.h>

class STPGMergeMip{
  private:
    STPGModel *probModel;
    SettingsHandler sh;
    set_obj best_sol;
    merge_sol groups;
    std::vector<qol::Variable> z;

  public:
    bool rel_gap;

    STPGMergeMip(const SettingsHandler sh, STPGModel *base_model,merge_sol &groups, set_obj &best_sol);

    void set_groups(merge_sol &new_groups){
        groups = new_groups;
    };

    void solve(set_obj &sol);

    void initMIPModel(qol::MIPSolver &mip);

    ~STPGMergeMip(){};

};

class STPGMergeCallback : public qol::Callback{
    STPGModel *probModel;
    const std::vector<qol::Variable> &z;
    int *countPtr;
    qol::MIPSolver *mipPtr;
    merge_sol groups;
    public:
    STPGMergeCallback(STPGModel *_probModel, const std::vector<qol::Variable> &_z,int *_countPtr,qol::MIPSolver *_mipPtr,const merge_sol &_groups):probModel(_probModel),z(_z),countPtr(_countPtr),mipPtr(_mipPtr),groups(_groups){};
    ~STPGMergeCallback(){};
    Status callback(Progress where);
};

#endif
