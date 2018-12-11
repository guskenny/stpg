#ifndef __QOL_SOLUTION__
#define __QOL_SOLUTION__

#include "QolData.h"

namespace qol {
  // a solution may represent a:
  // (a) primal feasible set of values
  // (b) primal infeasible (relaxed/lower-bound) set of values
  class SolutionBase{
  private:
    double obj;				 // objective
  public:
    double getObjective () const {return obj;}
	void setObjective (double _o) { obj = _o;}
	virtual ~SolutionBase() {} // just to declare destructor virtual  

    SolutionBase &operator <= (SolutionBase &s){ return this->obj<=s.obj ? *this : s;}

    virtual void writeSolutionTxt(std::ostream& os) const {};
    virtual void writeSolutionXML(std::ostream& os) const {};
    virtual void writeSolutionCVS(std::ostream& os) const {};
  };
}
#endif
