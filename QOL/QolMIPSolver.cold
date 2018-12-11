#include "QolMIPSolver.h"
namespace qol {
void MIPSolver::extractSolution () {
    if (sln) delete(sln);
    if (_type==qol::FULL)
      sln = new qol::SolutionFull();
    else
      sln = new qol::SolutionPrimal();
    sln->setObjective (getObjective ());
    const Index numVars = nVar();
    for (Index v=0; v<numVars; ++v) {
      if (_type==qol::PRIMAL) {
        ((SolutionPrimal*)sln)->push_back (getPrimal (qol::Variable(v)));
      }else {
        qol::ExtendedVariable variable (qol::Variable(v), getVarName(qol::Variable(v)),
          getPrimal(qol::Variable(v)), getDual(qol::Variable(v)));
        ((SolutionFull*)sln)->push_back (variable);
      }
    }
  }
};
