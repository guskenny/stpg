#include "LagrangianModel.h"

namespace qol {

  Index LagrangianModel::nVar() const {
    return (vars.size());
  }

  Index LagrangianModel::nConstr()const {
    return (cnstrAtProblem.size());
  }  // explicit/linking constraints

  qol::ConstraintMIP LagrangianModel::getConstraintRow(Constraint c) {
    return (problems[cnstrAtProblem[c].first]->getConstr (
      cnstrAtProblem[c].second));
  }

  void LagrangianModel::addProblem (qol::MIPSolver* p) {
    const Index numVars = p->nVar();

    Index idxProb = problems.size();
    for (Index i = 0; i<numVars ; ++i) {
      qol::Variable v = p->getVar (i);
    //  varAtProblem.push_back (std::make_pair (idxProb, i));
      vars.push_back (v);
    }
    cnstrAtProblem.reserve (nConstr()+p->nConstr());
    for (Index i = 0; i<p->nConstr(); ++i) {
      qol::ConstraintMIP c = p->getConstr(i);
      cnstrAtProblem.push_back(std::make_pair (idxProb, i));
    }
    problems.push_back (p);
  }

  ///SolveRelaxed
  Status LagrangianModel::solveRelaxed() {
    std::vector <Status> status (nSubprob());
    DblVec lb (nSubprob(), 0.0);

    for (Index p = 0; p<nSubprob (); ++p) {
      qol::MIPSolver * sprob = problems[p];
      ////
      //Solve each sproblem using CPLEX, Gurobi,...
      ////
      status[p] = sprob->solveRelaxed();

      lb[p] = problems[p]->getObjective();
    }

    LB = lb.sum();
    return qol::HEURISTIC;
  }

  ///Solve
  Status LagrangianModel::solve() {
    lagMulti.clear();
    lagMulti.resize (lnkCnstrs.size(), 0.0);
    DblVec viol (lnkCnstrs.size(), 0.0); // = limits;
    bool isFeasible = false;
    double modViol;
    //double lagCost;
    double step = 0.3;  //
    double multi = 0.9; //
    int noProgressCnt = 0;

    std::vector <Status> status (nSubprob());
    DblVec lb (nSubprob(), 0.0);

    dual.resize (nSubprob(), 0.0);
    primal.resize (nVar());

    for (Index i=0; i<lnkCnstrs.size(); ++i) {
      if (lnkCnstrs[i].rhs<dualLB[i]) {
        lnkCnstrs[i].rhs=dualLB[i];
      }
      if (lnkCnstrs[i].rhs>dualUB[i]) {
        lnkCnstrs[i].rhs=dualUB[i];
      }
    }

    const Index maxIter = Index(params.getParamValue (LAGITER));
    timer.setTimeLimit (params.getParamValue (LAGTIMER));

    for (Index it = 0; it<maxIter && !timer.timeLimitReached() && step>1e-8; ++it) {
      //backup
      std::vector < qol::ConstraintRow > copyLnkCnstrs = lnkCnstrs;

      for (Index i=0; i<lnkCnstrs.size(); ++i) {
        for (Index ii=0; ii<lnkCnstrs[i].size(); ++ii) {
          lnkCnstrs[i][ii].first-=lagMulti[i];
        }
      }
      for (Index p = 0; p<nSubprob (); ++p) {
        qol::MIPSolver  *sprob = problems[p];
        ////
        //Solve each sproblem using CPLEX, Gurobi,...
        ////
        for (Index i=0; i<lnkCnstrs.size(); ++i) {
          qol::Row row;
          for (Index ii=0; ii<lnkCnstrs[i].size(); ++ii) {
            qol::Variable v = problemVar (p, copyLnkCnstrs[i][ii].second);
            if (!v.valid()) continue;
            
            row.lhs.push_back (CoeffIdxPair(copyLnkCnstrs[i][ii].first, v));
          }
          if (row.lhs.size()==0) continue;
          row.op = copyLnkCnstrs[i].sense;
          row.rhs = copyLnkCnstrs[i].rhs;
          sprob->addConstraint (row);
        }
        
        status[p] = sprob->solveRelaxed();

        calculateViolations (p, viol);

        lb[p] = problems[p]->getObjective();
      }
      //restore backup
      lnkCnstrs = copyLnkCnstrs;

      //reducedCost ();
      double lastLB = lb.sum();
      if(lastLB < LB-0.001){
			  if(++noProgressCnt > 15){
          return qol::INFEASIBLE;
			  }
			  multi *= 0.7;
		  }

      LB = std::max (LB, lastLB);

      modViol = 0.0;
      for (Index i=0; i<viol.size(); ++i) {
        if (( viol[i] > 0.0 && lagMulti[i] < 0.0 ) || ( viol[i] < 0.0 ))
				  modViol += viol[i] * viol[i];
      }

      if (modViol>0.0) {
		    modViol = sqrt(modViol);
        if (UB!=inf)
		      step = multi * (UB-LB)/modViol;
        else
          step*=multi;
      }

      multi *= 0.99;

      for (Index i=0; i<lagMulti.size(); ++i) {
        lagMulti[i] = std::max (0.0, lagMulti[i]-step*viol[i]);
      }

      for (Index l=0; l<lnkCnstrs.size(); ++l) {
        lnkCnstrs[l].rhs = std::max(dualLB[l], 
          std::min(dualUB[l],lnkCnstrs[l].rhs));
      }

      if (modViol==0.0) {
        return qol::OPTIMAL;
      }
    }

    return qol::HEURISTIC;
  }

  // fix/unfix variable bounds
  void LagrangianModel::setVarBounds(Variable var,double lb,double ub){
    Index problem = InvalidIndex;
    qol::Variable v = modelVar (var, problem);

    if (problem==InvalidIndex) return;
    qol::MIPSolver &prob = subproblem (problem);

    prob.setVarLB(v,lb);
    prob.setVarUB(v,ub);
  }

  /// setting variable cost may be sensible as part of Lagrangian
  /// scheme or maybe to do some perturbation in a Benders scheme.
  void LagrangianModel::setVarCost(Variable var,double obj, Index problem){ // objective coefficient
    if (problem==InvalidIndex) return;
    Variable v = modelVar (var, problem);
    qol::MIPSolver &prob = subproblem (problem);
    prob.setObjCoeff (v, obj); 
  }

  /// <summary>
  /// Calculate the reduced cost for each variable given the duals.
  /// Note Lagrangian is c x + lambda (b- A x) so need to
  /// subtract duals!
  /// Reduced cost(j,t) = c_j,t - \sum_k=1^duration lambda_(t-k) * R_j
  /// </summary>
  void LagrangianModel::reducedCost() {
    redCost.clear();

    redCost.resize (nVar(), 0.0);

    for (Index p=0; p<nSubprob (); ++p) {
      qol::MIPSolver * sprob = problems[p];
      Index var = 0;

      const Index nvar = sprob->nVar(); //sprob->matrix.size() ;

      for (Index j = 0; j< nvar; ++j) {

        const Index v = modelVar (sprob->getVar (j), p);

        redCost [v] = 0.0;

        for (int t = 0; t<sprob->matrix[j].size(); ++t) {
          //redCost[v] += sprob->getObjCoeff (j)*std::max 
          //                  (t-sprob->getVarUB(j), 0.0);
          double accRedCost = 0.0;

          #pragma omp parallel for reduction(-:accRedCost)
          for (int i = 0; i<(int)lnkCnstrs.size(); ++i) {
            for (Index ii=0; ii<lnkCnstrs[i].size(); ++ii) { 
              if (lnkCnstrs[i][ii].second== Variable(v))
                accRedCost -= sprob->getDual (qol::Constraint 
                              (sprob->matrix[j][t].second))*lnkCnstrs[i].rhs;
            }
          }
          redCost [v] += accRedCost;
          ++var;
        }
      }          
    }
  }

  void LagrangianModel::calculateViolations (Index p, DblVec& viol) {
    const qol::MIPSolver * sprob = problems[p];
    //Calculate violations
    const Index nvar = sprob->matrix.size(); //sprob->nVar();
    //const Index nconst = sprob->nConstr();
    //	  const double tol = problems[p]->getTolerance(); //LOOK
    const double tol = 0.0; //LOOK
    Index prob = p;

    for (Index i=0; i<lnkCnstrs.size(); ++i) {
      double val = 0.0;
      for (Index ii=0; ii<lnkCnstrs[i].size(); ++ii) {
        Variable v = problemVar (prob, lnkCnstrs[i][ii].second); //, prob);
        if (!v.valid()) continue;
        val += lnkCnstrs[i][ii].first*sprob->getPrimal (v);        
      }

      switch (lnkCnstrs[i].sense) {
        case 'L': if (val>lnkCnstrs[i].rhs) viol[i] = -1;
          break;
        case 'E': if (val!=lnkCnstrs[i].rhs) viol[i] = -1;
          break;
        default: if (val<lnkCnstrs[i].rhs) viol[i] = -1;
          break;
      }
    }
  }

  double LagrangianModel::getObjectiveBound () {
    double lb = 0.0;
    for (Index p = 0; p<nSubprob (); ++p) {
      lb+= problems[p]->getObjectiveBound();
    }
    return lb;
  }
}

