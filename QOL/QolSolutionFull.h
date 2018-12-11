#ifndef __QOL_SOLUTION_FULL__
#define __QOL_SOLUTION_FULL__

#include "QolData.h"
#include "QolSolutionBase.h"

namespace qol {
  // This solution stores:
  // (a) primal and dual feasible set of values (plus the variable name)
  // (b) primal and dual infeasible (relaxed/lower-bound) set of values
  class SolutionFull : public SolutionBase, public ExtendedVarVec {
  public:
    virtual void writeSolutionTxt(std::ostream& os) const {
      ExtendedVarVec::const_iterator it;
      for (it=(*this).begin(); it!=(*this).end(); ++it) { 
        os << (*it) << std::endl;
      }
    };
    virtual void writeSolutionXML(std::ostream& os) const {};
    virtual void writeSolutionCVS(std::ostream& os) const {};
  };
}
#endif
