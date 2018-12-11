// So far this is an empty class. It will implement the ACO algorithm. Just added to help with the 
// organisation of the whole library.

#ifndef __ACO_OPTIMISER__
#define __ACO_OPTIMISER__

#include "QolOptimiser.h"

namespace qol {
  class ACO : public Optimiser {
  public:
    ACO (qol::ModelType type, qol::FormulationType form): Optimiser(type, form){}
    virtual ~ACO(void) {}
  protected:
    /// interaction with solver
    virtual qol::Status AlgorithmSolve () ; // do heuristic/primal solving
  };

} //namespace qol

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
