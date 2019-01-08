#ifndef __STPGMip_H__
#define __STPGMip_H__

#include <boost/format.hpp>
#include "QOL/QolMIP.h"
#include "QOL/QolColFormulation.h"
#include "QOL/CplexFormulation.h"
#include "QOL/GurobiFormulation.h"
#include "QOL/CpuTimer.h"
#include <math.h>
#include <map>
#include <vector>
#include <set>
#include "SettingsHandler.h"
#include <STPGModel.h>
#include <string>
#include <sstream>
#include <util.h>

class STPGMip{
  private:
    STPGModel *probModel;
    SettingsHandler sh;
    set_obj best_sol;
    std::vector<qol::Variable> x;
    std::vector<qol::Variable> y;

  public:
    bool rel_gap;

    STPGMip(const SettingsHandler sh, STPGModel *base_model, set_obj &best_sol);

    void solve(set_obj &sol);

    void initMIPModel(qol::MIPSolver &mip);

    ~STPGMip(){};

};

class STPGCallback : public qol::Callback{
    STPGModel *probModel;
    const std::vector<qol::Variable> &x;
    const std::vector<qol::Variable> &y;
    int *countPtr;
    qol::MIPSolver *mipPtr;
    public:
    STPGCallback(STPGModel *_probModel,const std::vector<qol::Variable> &_x, const std::vector<qol::Variable> &_y,int *_countPtr,qol::MIPSolver *_mipPtr):probModel(_probModel),x(_x),y(_y),countPtr(_countPtr),mipPtr(_mipPtr){};
    ~STPGCallback(){};
    Status callback(Progress where);
};

#endif
