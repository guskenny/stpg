#include "PSO.h"

namespace qol {
  void Particle::initialise(size_t numVar,size_t numConstr, int nParticles) {
    dual.resize(numConstr,0.0);
    perturb.resize(numVar,0.0);
    x.resize(numVar,0);
    viol.resize(numConstr,0.0);
    rc.resize(numVar,0.0);
    isFeasible = false;
    dVel.clear();
    // set up first particle with zero multiple, second mostly zero
    // etc through to last one with lots of non-zero multipliers
    for (int i = 0; i<numConstr; ++i) {
      dVel.push_back (rand->uniform() < idx/double (nParticles)
        ? 0.0 : -rand->uniform(0,maxCost));
    }
    lbCost = -inf;
    ubCost = inf;
    pVel.clear();
    for (int i = 0; i<numVar; ++i) {
      pVel.push_back (rand->uniform() < idx/double (nParticles)
        ? 0.0 : -rand->uniform(0,maxCost));
    }
    //	model->setLimits (viol);
  }


  qol::Status Particle::heuristic (int iter, int heurFreq) {
    // only run heuristic if not feasible
    if ( heurFreq> 0 && iter % heurFreq == 0 && ! isFeasible) {
    }
    return qol::HEURISTIC;
  }

  qol::Status Particle::solveInitial () {
    qol::Status status;

    if (!model) return qol::FAILED;

    status = model->solveRelaxed();
    lbCost = model->getLB();

    status = model->solve();

    lbCost = model->getLB();
    ubCost = model->getUB();

    if (status!=HEURISTIC) {
      log->warning ("Unsolved problem");
      return status;
    }

    if (lbCost>bestLbCost) {
      bestLbCost = lbCost;
      bestIter = 0;
    }
    else { //This should never happen
      char str [250];
      sprintf (str, "Best lower bound not improved in the first iteration!  %d  Particle \n", int(idx));
      log->message (str);
    }

    if (ubCost < bestUbCost) {
      bestUbCost = ubCost;
    }

    pVel = 0;
    dVel = 0;

    return qol::HEURISTIC;
  }

  ///<summary>
  /// Main method. Calls the lagrangian or any other model and solves it
  /// If this is the first iteration (-1), it just calls solveInitial
  ///</summary>
  qol::Status Particle::solve (int lbcheckfreq, double & factor, int iter) {
    qol::Status status;

    if (iter==-1) return solveInitial();

    if( iter % lbCheckFreq != 0){
      rc += perturb; 
      status = model->solve ();
      if (status==FAILED) return status;
    }
    else {
      DblVec zero (perturb.size(), 0.0);
      std::swap (perturb, zero);
      status = model->solve ();
      std::swap (perturb, zero);
      if ((status==FAILED) || (status==INFEASIBLE)) {
        isFeasible = false;
        return status;
      }
      if (status==HEURISTIC) isFeasible = true;
    }

    lbCost = model->getLB();
    ubCost = model->getUB();

    if (ubCost < bestUbCost) {
      bestUbCost = ubCost;
    }

    if( iter % lbcheckfreq== 0){
      if (lbCost > bestLbCost) {
        bestLbCost = lbCost;
        bestIter = iter;
      }
      else {
        factor = std::max(subgradFactorMin,subgradFactorMult*factor);
      }
    }

    return qol::HEURISTIC;
  }

  /** Find the "best" solution that is feasible with respect to the
  *  current constraint. Here best means low cost in terms of the
  *  reduced cost. The output is a set of variables and values that
  *  make the constraint feasible. Used to determine the perturbation
  *  direction.
  *  @param constraint (IN) the constraint index
  *  @param feas (OUT) a set of (<index>,<0/1>) pairs for all variables
  *         involved in the constraint that are feasible w.r.t.
  *         this one constraint.
  */
  void Particle::fixConstraint (int col, SparseVec & feas, DblVec x) {
    feas.clear();
    std::vector <IdxVal> cnstrs;
    IntVec cnstr_idx (model->nConstr(), -1);
    qol::LagrangianModel* tmodel = (qol::LagrangianModel*)model;

    for (Index p=0; p<model->nSubprob (); ++p) {
      qol::MIPSolver & sprob = tmodel->subproblem (p);
      for (int j = 0; j<int (sprob.matrix.size()); ++j) {
        int firstcol = j+sprob.matrix[j].size()+col+1;
        int lastcol = std::min (double(firstcol+sprob.getVarLB (j)),
            double(j*sprob.matrix[j].size()+sprob.matrix[j].size()-1));
        for (int i = firstcol; i<lastcol ; ++i) {
          if (x[i]<=0) continue;

          IdxVal temp (sprob.getVarUB (j), j);
          cnstrs.push_back (temp);
          cnstr_idx[j] = i;
          break;
        }
      }
    }

    if (cnstrs.size() > tmodel->nSubprob()) {
      char str [100];
      sprintf (str, "*ERROR: Primal solution not feasible \n");
      log->error (str);
    }

    std::sort (cnstrs.begin(), cnstrs.end());
    double violVal = viol [col];
    while (violVal<0) {
      int max = cnstrs.back().first;
      int var = cnstrs.back().second;
      violVal += tmodel->getLnkConstraint (var).rhs;

      feas.add(cnstr_idx[var],0);
      cnstrs.pop_back();
    }

    while (!cnstrs.empty ()) {
      int var = cnstrs.back().second;
      feas.add (cnstr_idx[var], 1);
      cnstrs.pop_back();
    }
  }

  //Update the vector with perturbation direction
  void Particle::perturbationDirection (double eps, DblVec dualLB,
    DblVec dualUB, DblVec x) {
      Index dSize = dual.size();
      perturbDir.clear();
      perturbDir.resize (dSize, 0.0);

      for ( int i=0; i<dSize; ++i ) {
        if (((viol[i]>eps) && (dualLB[i]<-eps)) ||
          ((viol[i]<-eps) && (dualUB[i]>eps))) {
            SparseVec feas;
            fixConstraint (i, feas, x);
            for( SparseVec::iterator it=feas.begin(); it!=feas.end(); ++it ){
              if(it->second >= 1) perturbDir[it->first] -= 1;
              if(it->second <= 0) perturbDir[it->first] += 1;
            }
        }
      }
  }

  // Update the velocity of the particle
  void Particle::updateVelocity (double eps, double subgradFactor,
    double globalFactor, double bestUB,
    double velFactor, DblVec bestDual,
    DblVec dualLB, DblVec dualUB,
    DblVec bestPerturb, DblVec bestx) {
      double norm = 0.0;
      Index dSize = dVel.size();
      Index pSize = pVel.size();

      perturbationDirection (eps, dualLB, dualUB, bestx);

      for(int i=0;i<dSize;++i) norm += viol[i]*viol[i];

      double rAsc = rand->uniform () ;  // 0,1 uniform
      double rGlobal = rand->uniform () * globalFactor;

      double stepSize;
      if (norm!=0) {
        stepSize = rand->uniform( subgradFactor/2.0,subgradFactor) *
          (bestUB-lbCost)  / norm;
      }
      else {
        stepSize = eps;
      }

      for(int i=0;i<dSize;++i) {
        dVel [i] = velFactor * dVel[i] + stepSize*viol[i]+
          rGlobal*(bestDual[i]-dual[i]);
      }

      for (int i=0; i<pSize; ++i) {
        pVel[i] = pVel[i]*velFactor +
          perturbDir[i]*rAsc*perturbFactor+
          perturbFactor*rGlobal*(1-2*bestx[i])+
          rGlobal*(bestPerturb[i]-perturb[i]);
      }
  }

  //Make a step
  void Particle::makeStep (DblVec dualLB, DblVec dualUB) {
    Index dSize = dVel.size();

    for(int i=0;i<dSize;++i){
      dual[i] += dVel [i];
      dual[i] = std::max(dualLB[i],std::min(dualUB[i],dual[i]));
    }

    perturb*=0.5;
    perturb+=pVel;        // add step in perturbation velocity
  }

  // This method needs to be called right after creating the PSO algorithm
  // In general, any algorithm needs to implement this method to
  // read the specific params of each algorithm

  void PSO::initialise(int argc,char **argv) {
    int currentarg = 0;
    while (currentarg<argc) {
      if (!strcmp (argv[currentarg], "--alg"))
      {
        currentarg++;
        if (!strcmp (argv[currentarg], "PSO")) {
          currentarg++;
          if (strcmp (argv[currentarg], "--params")) return;
          currentarg++;
          while ((currentarg+2<argc) /*&& ()*/) {
            if (!strcmp (argv[currentarg], "-iterations"))
              params.setParamVal (qol::ITERATIONS, atoi (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-subgradfactor"))
                params.setParamVal (qol::SUBGRADFACTOR, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-velocity"))
              params.setParamVal (qol::VELOCITY, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-globalfactor"))
              params.setParamVal (qol::GLOBALFACTOR, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-heurfreq"))
              params.setParamVal (qol::HEURISTIC_FREQ, atoi (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-lbcheckfreq"))
              params.setParamVal (qol::LBCHECKFREQ, atoi (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-perturbfactor"))
              params.setParamVal (qol::PERTURBFACTOR, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-subgradfmin"))
              params.setParamVal (qol::SUBGRADFMIN, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-subgradfmult"))
              params.setParamVal (qol::SUBGRADFMULT, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-particles"))
              params.setParamVal (qol::NPARTICLES, atoi (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-absgap"))
              params.setParamVal (qol::ABSGAP, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-relgap"))
              params.setParamVal (qol::RELGAP, atof (argv[currentarg+1]));
            if (!strcmp (argv[currentarg], "-eps"))
              params.setParamVal (qol::EPS, atof (argv[currentarg+1]));

            currentarg+=2;
          }
        }
        else {
          log->error ("PSO called when PSO was not specified!\n");
        }
        break;
      }
      currentarg++;
    }

    swarm.clear();

    double subgradFactorMin = params.getParamValue (SUBGRADFMIN);
    double subgradFactorMult = params.getParamValue (SUBGRADFMULT);

    int nParticles = int (params.getParamValue (qol::NPARTICLES));
    for (int i=0; i<nParticles; ++i) {
      swarm.push_back (new Particle(i, 0.0, model, type, nParticles, 0.0,
        subgradFactorMin, subgradFactorMult, log));
    }

    bestUB = inf;
    bestLB = -inf;
    bestDual.clear();
    bestViol.clear();
    dualUB.clear();
    dualLB.clear();
    bestPerturb.clear();
    bestx.clear();

    bestDual.resize (model->nConstr(), inf);
    bestViol.resize (model->nConstr(), 0.0);
    bestPerturb.resize (model->nVar(), 0.0);
    bestx.resize (model->nVar(), inf);
    dualUB.resize (model->nConstr(), inf);
    dualLB.resize (model->nConstr(), -inf);

    absGap = params.getParamValue (ABSGAP);
    relGap = params.getParamValue (RELGAP);
    globalFactor = params.getParamValue (GLOBALFACTOR);
    velFactor = params.getParamValue (VELOCITY);
    perturbFactorVal = params.getParamValue (PERTURBFACTOR);
    subgradFactorVal = params.getParamValue (SUBGRADFACTOR);
    subgradFactorMin = params.getParamValue (SUBGRADFMIN);
    subgradFactorMult = params.getParamValue (SUBGRADFMULT);
    lbcheckfreq = int (params.getParamValue (LBCHECKFREQ));
    nParticles = int(params.getParamValue (NPARTICLES));
    heurFreq = int (params.getParamValue (HEURISTIC_FREQ));
    status = qol::FAILED;

    initialised = true;
  }

  ///<summary>
  /// PSO main method. Call the solve method of each particle,
  /// updates the velocity of each particle, make the steps,...
  ///</summary>
  qol::Status PSO::AlgorithmSolve () {

    for (swarmIt=swarm.begin(); swarmIt!=swarm.end();++swarmIt) {
      (*swarmIt)->setLBCheckFreq (lbcheckfreq);
    }

    //---------- initialisation -------------------------
    for (swarmIt=swarm.begin(); swarmIt!=swarm.end();++swarmIt) {
      double temp = 0.0;
      status = (*swarmIt)->solve(0, temp, -1);
      if (status == qol::FAILED)
        return status;
    }

    updateBest();

    DblVec subgradFactor(nParticles,subgradFactorVal);
    DblVec perturbFactor(nParticles,perturbFactorVal);

    //Iterations
    int numIter = int (params.getParamValue (ITERATIONS));

    for ( int iter = 0 ; (iter < numIter) && (bestLB + absGap <= bestUB)
      && (fabs(bestUB-bestLB/bestUB) > relGap) ; ++iter ) {

        for (swarmIt=swarm.begin(); swarmIt!=swarm.end();++swarmIt) {

          Index idx = std::distance (swarm.begin(), swarmIt);

          // --------- calculate step size/direction ----------
          (*swarmIt)->updateVelocity (eps, subgradFactor[idx], globalFactor,
                    bestUB, velFactor, bestDual, dualLB,
                    dualUB, bestPerturb, bestx);

          //---------- make a step ----------------------------
          (*swarmIt)->makeStep (dualLB, dualUB);

          //---------- solve subproblems ----------------------
          status = (*swarmIt)->solve (lbcheckfreq, subgradFactor[idx], iter);

          if (status==qol::FAILED) continue;

          //---------- heuristic ------------------------------
          status = (*swarmIt)->heuristic (iter, heurFreq);
          if (status==qol::FAILED) continue;

          //---------- print info -----------------------------
          printIterStatus (iter, (*swarmIt));
        } //swarmIt
        updateBest ();
    } //iter

    return qol::HEURISTIC;
  }

  ///<summary>
  /// Method to store the best solution found so far
  /// Iterates over all of the particles asking for the best UB and the best LB
  ///</summary>
  void PSO::updateBest () {
    for (swarmIt=swarm.begin(); swarmIt!=swarm.end();++swarmIt) {
      if ((*swarmIt)->getLB()>bestLB) {
        bestLB = (*swarmIt)->getLB();
        bestDual = (*swarmIt)->getDual();
      }

      if ((*swarmIt)->getUB()<bestUB) {
        bestUB = (*swarmIt)->getUB();
        bestPerturb = (*swarmIt)->getPerturb();
      }
    }
  }

  void PSO::setViolations (DblVec _viol) {
    //	viol = _viol;
    for (swarmIt=swarm.begin(); swarmIt!=swarm.end();++swarmIt) {
      (*swarmIt)->setViolations(_viol);
    }
  }

  ///<summary>
  /// Print stats at the end of each iteration for each particle
  ///</summary>
  void PSO::printIterStatus (int iter, Particle *p) {
    int printLevel = int (params.getParamValue (VERBOSITY));
    int printFreq = int (params.getParamValue (PRINTFREQ));

    if( printLevel > 1 && iter % printFreq == 0){
      char str[200];
      sprintf(str,"\tp%02d: LB=%g UB=%g feas=%d minViol=%g\n",
        p->getIdx (), p->getLB (), p->getUB(), p->getIsFeasible(),
        p->getViolations().min());
      log->message (str);

      if(printLevel > 2){
        //printf("\t\tstepSize=%g randGlobal=%g\n",
        //	   stepSize,randGlobal);
        sprintf(str,"\t\tRedCst %g - %g\n", p->getReducedCost().min(),
          p->getReducedCost().max());
        //printf("\t\tRedCst %g - %g\n", p->getReducedCost().min(),
        //                               p->getReducedCost().max());
        log->message (str);
        sprintf (str, "\t\tdVel %g - %g\n",p->getDualVel().min(),
          p->getDualVel().max());
        log->message (str);
        sprintf (str, "\t\tpDir %g - %g\n",p->getPerturbDir().min(),
          p->getPerturbDir().max());
        log->message (str);
        sprintf (str, "\t\tpVel %g - %g\n",p->getPerturVel().min(),
          p->getPerturVel().max());
        log->message (str);
        sprintf (str, "\t\tdual %g - %g\n",p->getDual().min(),
          p->getDual().max());
        log->message (str);
        sprintf (str,"\t\tpert %g - %g\n",p->getPerturb().min(),
          p->getPerturb().max());
        log->message (str);
      }
    }
  }
} //namespace qol

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
* Local Variables:
* tab-width: 4
* eval: (c-set-style "stroustrup")
* End:
*/
