#ifndef __QOL_SOLUTION_PRIMAL__
#define __QOL_SOLUTION_PRIMAL__

#include "QolData.h"
#include "QolSolutionBase.h"

namespace qol {
  // This solution stores:
  // (a) primal feasible set of values
  // (b) primal infeasible (relaxed/lower-bound) set of values
  class SolutionPrimal : public SolutionBase, public DblVec {
  public:
    virtual void writeSolutionTxt(std::ostream& os) const {
      DblVec::const_iterator it;
      for (it=(*this).begin(); it!=(*this).end(); ++it) { 
        os << (*it);
      }
    };
    virtual void writeSolutionXML(std::ostream& os) const {};
    virtual void writeSolutionCVS(std::ostream& os) const {};
  };
}
#endif
