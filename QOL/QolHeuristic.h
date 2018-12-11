#ifndef __QOL_HEURISTIC__
#define __QOL_HEURISTIC__

#include "QolSolutionBase.h"
#include "QolUtil.h"

namespace qol {
  class Heuristic {
  protected:
    // NeighbourhoodSearch - need a primal feasible solution
    // Repair - need a primal, not necessarily feasible solution
    // Constructive - just need some weights (costs) for each variable
    HeuristicType type;
  public:

    virtual Status operator() (const qol::DblVec &input, qol::SolutionBase &soln) = 0;
    
    virtual HeuristicType getType() const {return type;}
    virtual void setType (HeuristicType _type) { type = _type; }

    Heuristic (HeuristicType _type) : type (_type) {}
    virtual ~Heuristic (){}

  };
}

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
