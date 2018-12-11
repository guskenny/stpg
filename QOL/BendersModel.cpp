#include "BendersModel.h"

namespace qol {

  Index BendersModel::nVar() const {
    return 0;
  }

  Index BendersModel::nConstr()const {
    return 0;
  } 

  qol::ConstraintMIP BendersModel::getConstraintRow(Constraint c) {
    return qol::ConstraintMIP ();
  }

  ///SolveRelaxed
  Status BendersModel::solveRelaxed() {
      return qol::HEURISTIC;
  }

  Status BendersModel::solveExact(){  // do branch and bound as well
    return qol::FAILED;
  }

  void BendersModel::setVarBounds(Variable var,double lb,double ub){ // fix/unfix
    /*Index problem = InvalidIndex;
    qol::Variable v = modelVar (var, problem);
    if (problem==InvalidIndex) return;
    qol::MIPSolver &prob = subproblem (problem);

    prob.setVarLB(v,lb);
    prob.setVarUB(v,ub);
    */
  }

  /// scheme or maybe to do some perturbation in a Benders scheme.
  void BendersModel::setVarCost(Variable var,double obj, Index problem){ // objective coefficient
  }
}
