/// Model is a way of specifying the problem to be solved but potentially
/// only implicitly - that is we may just provide a method for solving
/// some of the subproblems rather than explicit constraints.
#ifndef __QOL_MODEL__
#define __QOL_MODEL__

#include <assert.h>
#include <list>
#include <vector>
#include <string>

#include "QolUtil.h"
#include "QolData.h"
#include "QolHeuristic.h"
#include "QolMIPSolver.h"

#include "CpuTimer.h"

namespace qol {
  /// Model is an abstract base class for that corresponds both to a formulation
  /// (which may only be implicit encoded in solvers rather than explicit) and
  /// a basic bit of solver code for solving subproblems.  It may
  /// (a) be inspected for number of subproblems, linking constraints etc
  /// (b) have its subproblems solved either exactly or LP relaxtion only
  /// (c) Potentially be modified eg by setting lagrangian vector or value
  ///     of linking variables in a Benders problem
  /// Depending on the type of model a "solve" may produce a lower bound
  /// (Lagrangian Model or solveRelaxed for MIPformulation) or an  upper bound
  /// (Benders Model, solveExact for MIPformulation).
  /// Note that the Model class may also be used to hold a subproblem of the
  /// original problem in which case the solve method(s) may implement a
  /// specialised solution algorithm (eg shortest path) for the simpler
  /// subproblem structure. Further nesting (having subsubproblems) is not
  /// explicitly disallowed but in the first instance will not really be
  /// catered for.


  class Model {
  protected:
    // current solution
    DblVec dual;		  // lagrangian/dual vector (length nConstr)
    DblVec redCost;		// reduced cost vector for dual
   // DblVec limits;    // for the violations

    double UB;			//upper bound
    double LB;			//lower bound

    /// each problem is represented by its formulation
    std::vector <qol::MIPSolver*> problems;

    std::vector<qol::Heuristic *>heuristic; // model specific heuristics
    std::vector<qol::Variable> vars;
    qol::Parameters params;
    CpuTimer timer;

    std::vector < Index > varAtProblem; //varAtProblem[i] = Problem where i is used
    std::vector < std::vector < Index > > problemsVars; //problemsVars[i][j] -> global index of the variable j in problem i
  public:
    Model (){
      problems.clear();
      heuristic.clear();
      dual.clear();
      redCost.clear();
      vars.clear();
      varAtProblem.clear();
      problemsVars.clear();
     // limits.clear();
      UB = inf;
      LB = -inf;
    }

    virtual ~Model() {}

    /// <summary>
    /// returns the model type
    /// </summary>
    virtual ModelType type() const { return QOLMODEL;}

    /// <summary>
    /// returns the number of subproblems
    /// </summary>
    virtual inline Index nSubprob() const { return problems.size();}

    virtual Index nConstr() const = 0; // explicit/linking constraints

    /// <summary>
    /// Adds a new subproblem to the model
    /// </summary>
    /// <param name="p">Subproblem</param>
    virtual void addProblem (qol::MIPSolver* p) {problems.push_back (p);}

    // number of variables in this model
    virtual Index nVar() const = 0;

    /*
    virtual DblVec getLimits () const { return limits; }
    virtual void setLimits (DblVec l) { limits = l; }
    virtual void setLimits (std::vector<double> l) { limits = l; }
    */

    // inspect the constraints: warning this is generally an expensive
    // operation as it may copy all of the constraint information from
    // elsewhere
    virtual qol::ConstraintMIP getConstraintRow (qol::Constraint c) = 0;

    virtual qol::MIPSolver &subproblem(Index s=0) {
      qol::MIPSolver* p; 
      (s<nSubprob()) ? (p = problems[s]) : (p = NULL) ; 
      return *p;
    }

    // is there a difference between solveRelaxed & solveExact ?
    // for unimodular matrices the LP relaxation is integer so no difference
    // also no difference if all variables are continuous
    //virtual bool relaxedIsExact() const = 0;

    // solveRelaxed should set reducedCost (and maybe dual unless dual is
    // an input as in Lagrangian relaxation)
    virtual Status solveRelaxed() = 0; // solve LP relaxation (root node) only

    virtual Status solveExact() = 0; // do branch and bound as well

    virtual Status solve () = 0; //this method should call either solveExact or solveRelaxed depending on the method

    // model modification / updates for model variables:
    virtual void setNumVar(Index n) = 0;	///< resize to given number of variables

    /// Variable bounds may be set due to global solution information that
    /// may have fixed the value in some way or at least reduced the
    /// domains. Alternatively it may be used as part of a say a branching
    /// or heuristic scheme to search for a solution and then "unfix" the
    /// bounds again.
    virtual void setVarBounds(qol::Variable var,double lb,double ub) = 0; // fix/unfix

    /// setting variable cost may be sensible as part of Lagrangian
    /// scheme or maybe to do some perturbation in a Benders scheme.
    /// Not sure it makes sense to change variable costs in a MIP formulation
    virtual void setVarCost(qol::Variable var,double obj, Index problem) = 0; // objective coefficient

    DblVec getReducedCost () const { return redCost; }

    void setReducedCost (DblVec rc) { redCost = rc; }

    void setUB (double _ub) { UB = _ub; }
    double getUB () const { return UB; }

    void setLB (double _lb) { LB = _lb; }
    double getLB () const { return LB; }

    void addVarAtProblem (qol::Variable v, Index p) {
      Index idx = v;
      if (varAtProblem.size()<= idx)
        varAtProblem.resize (idx+1);
      varAtProblem[idx] = p;
      if (problemsVars.size()<=p)
        problemsVars.resize (p+1);
      problemsVars[p].push_back (v);

      //varAtProblem[idx].second = v;
    }

  protected:
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
    virtual Variable problemVar(Index prob, qol::Variable modVar) const = 0;

    /// Given problem variable, tell us which model variable corresponds
    /// to it. Not sure if this is needed but could be used to test if
    /// modelVar(x).valid() and hence whether the subproblem interacts
    /// with the global problem variablle x
    virtual Variable modelVar(const qol::Variable probVar, Index problem) const  = 0;

  };
}	// end qol namespace
#endif
