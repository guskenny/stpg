#include "QolOptimiser.h"
///TODO
//
namespace qol {
  void Optimiser::writeLP(const std::string filename) {
    for (Index idx = 0; idx<model->nSubprob(); ++idx) {
      std::stringstream fname;
      fname << filename << "_" << idx;
      model->subproblem (idx).writeLP (fname.str().c_str());
    }
  }
}

