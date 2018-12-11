#ifndef __LAGRANGIAN_MODEL_H__
#define __LAGRANGIAN_MODEL_H__

#include "QolModel.h"
#include <assert.h>

namespace qol {
  class LagrangianModel: public Model {

    std::vector < qol::ConstraintRow > lnkCnstrs;
    DblVec primal;

    ///<summary>
    /// Array of ranges for the linking constraints associated to a given subproblem.
    /// The var [i] will be at problem varAtProblem[i].first, with index varAtProblem[i].second
    ///</summary>
    //std::vector < std::pair <Index, qol::Variable > > varAtProblem;
    ///<summary>
    /// Array to know the subproblem where the constraint 'i' (at cnstrAtProblem[i])
    /// can be found. The problem is cnstrAtProblem[i].first while the position at the
    /// problem is cnstrAtProblem[i].second
    std::vector < std::pair <Index, qol::Constraint > > cnstrAtProblem;

  private:
    /// The dual lower bound (0 for >= constraints otherwise -INF)
    /// The dual upper bound (0 for <= constraints otherwise +INF)
		DblVec dualLB, dualUB;
    DblVec lagMulti;
  public:
    LagrangianModel (): Model() {
      vars.clear();
      lagMulti.clear();
      varAtProblem.clear();
      cnstrAtProblem.clear();
      primal.clear();
      dualLB.clear();
      dualUB.clear();
      redCost.clear();
      lnkCnstrs.clear();
      params.addParam (qol::LAGITER, "lagrangian-iterations", 100, 0, 1000, 100);
      params.addParam (qol::LAGTIMER, "lagrangian-timer", 200, 0, 1000, 200);
    }

    virtual ~LagrangianModel (){}

    /// <summary>
    /// returns the model type
    /// </summary>
    virtual ModelType type() const {return LAGRANGIAN;}

    // number of variables in this model
    virtual Index nVar() const;

    //still to implement
    virtual void setNumVar(Index n) {};

    /// <summary>
    /// returns the number of constraints
    /// </summary>
    virtual Index nConstr()const;

    virtual qol::ConstraintMIP getConstraintRow (Constraint c);

    virtual Status solveRelaxed(); // solve LP relaxation (root node) only
    virtual Status solveExact() { return qol::FAILED; };  // do branch and bound as well

    virtual Status solve ();

    virtual void setVarBounds(Variable var,double lb,double ub); // fix/unfix

    /// setting variable cost may be sensible as part of Lagrangian
    /// scheme or maybe to do some perturbation in a Benders scheme.
    virtual void setVarCost(Variable var,double obj, Index problem); // objective coefficient

    double getObjectiveBound ();

    /// <summary>
    /// Mapping between variables in the problem and the model formulation
    /// problemVar tells us for each variable in the subproblem the
    /// corresponding variable in prob (if it exists)
    /// Typically expect that the problme variables are identical to the
    /// model variables (as assumed here) or a subset of the model
    /// variables. However more complicated mappings are possible and
    /// allowed for here. In particular if the Model is for a subproblem
    /// of the original problem we may only be setting a subset of the
    /// variables that is actually determined by the solutiono to this
    /// subproblem
    /// </summary>
    virtual qol::Variable problemVar(Index prob, qol::Variable modVar) const {
      for (Index i=0; i<problemsVars[prob].size(); ++i)
        if (problemsVars[prob][i] == Index(modVar))
          return qol::Variable(i);
      return qol::Variable (InvalidIndex);
    }

    /// <summary>
    /// Given problem variable, tell us which model variable corresponds
    /// to it. Not sure if this is needed but could be used to test if
    /// modelVar(x).valid() and hence whether the subproblem interacts
    /// with the global problem variablle x
    /// </summary>
    virtual qol::Variable modelVar(const qol::Variable probVar, Index problem) const  {
      if (Index(probVar)>=problemsVars[problem].size())
        return qol::Variable (InvalidIndex);
      return problemsVars [problem][probVar];
    }

    /// <summary>
    /// Add a subproblem to the model:
    ///     - Add the constraints of the problem to cnstrAtProblem
    ///     - Add the variables of the problem to varAtProblem
    /// </summary>
    virtual void addProblem (qol::MIPSolver* p) ;

    void reducedCost();

    qol::ConstraintRow getLnkConstraint (int idx) const { return lnkCnstrs[idx]; }

    /// <summary>
    /// set the linking constraints row
    /// </summary>
    /// <param name="r">
    /// vector of linking constraints. Each index 'i' contains the
    /// linking contraints of the subproblem 'i'
    /// </param>
    void setLnkConstraints (std::vector < qol::ConstraintRow > r) {
      dualLB.resize (r.size(), -inf);
      dualUB.resize (r.size(), inf);
      lnkCnstrs = r;
    }

    /// <summary>
    /// set the linking constraints row
    /// </summary>
    /// <param name="r">
    /// vector of linking constraints. Each index 'i' contains the
    /// linking contraints of the subproblem 'i'
    /// </param>
    void addLnkConstraints (std::vector < qol::ConstraintRow > r) {
      dualLB.resize (dualLB.size()+r.size(), -inf);
      dualUB.resize (dualUB.size()+r.size(), inf);
      for (int i=0; i<r.size(); ++i) {
        lnkCnstrs.push_back (r[i]);
   //     for (Index p = 0; p<nSubprob (); ++p) {
   //       problems[p]->addConstraint (Row (r[i]));
   //     }
      }
    }

    protected:
        void calculateViolations (Index p, DblVec& viol);
  };
}
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
