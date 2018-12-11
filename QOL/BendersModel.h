#ifndef __BENDERS_MODEL_H__
#define __BENDERS_MODEL_H__

#include "QolModel.h"

namespace qol {
  class BendersModel: public Model {

  private:
	  DblVec viol;
  public:
    BendersModel (): Model (){}
    virtual ~BendersModel (){}

    virtual ModelType type() const {return BENDERS;};

    virtual Index nConstr() const ; // explicit/linking constraints

    virtual qol::ConstraintMIP getConstraintRow(Constraint c);

    virtual Status solveRelaxed();

    virtual Status solveExact();

    virtual Status solve () { return solveRelaxed(); }

    virtual void setVarBounds(Variable var,double lb,double ub); // fix/unfix

    /// setting variable cost may be sensible as part of Lagrangian
    /// scheme or maybe to do some perturbation in a Benders scheme.
    virtual void setVarCost(Variable var,double obj, Index problem); // objective coefficient

    virtual Index nVar() const;

	//virtual std::vector <double> getLimits () const { return viol; }

    virtual void setNumVar(Index n) {};

    /// <summary>
    /// Given problem variable, tell us which model variable corresponds
    /// to it. Not sure if this is needed but could be used to test if
    /// modelVar(x).valid() and hence whether the subproblem interacts
    /// with the global problem variablle x
    /// </summary>
    virtual qol::Variable modelVar(const qol::Variable probVar, Index problem) const  {
      return qol::Variable();
    }

    virtual qol::Variable problemVar(Index prob, qol::Variable modVar) const {
      return qol::Variable();
    }
  };
}

#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
