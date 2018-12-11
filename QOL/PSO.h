#ifndef __PSO_OPTIMISER__
#define __PSO_OPTIMISER__

#include "QolOptimiser.h"
#include "QolRandom.h"

namespace qol {

  /** A particle defines a current solution in both primal and dual space.
  * It might be possible to remove this class and work directly with the model from 
  * the swarm, but I think conceptually is clearer this way
  */ 
  class Particle {
  protected:
    Index idx;
    double maxCost;
    int lbCheckFreq;
    qol::Model* model;
    qol::ModelType type;
    Log * log;

    /** multiplier for perturbations (should be small) */
    double perturbFactor;

    DblVec perturbDir;

    DblVec dVel;		        ///< dual velocity
    DblVec pVel;		        ///< perturbation velocity
    DblVec dual;	          ///< lagrange multiplier for each relaxed constraint
    DblVec perturb;         ///< perturbation for each variable
    DblVec viol;            /// The viol vector contains the violation from the primal solution
    /// calculated as b-Ax for the relaxed constraints

  private: 
    void initialise(size_t numVar,size_t numConstr, int nParticles);

    qol::Status solveInitial ();

  public:
    IntVec x;	                  ///< primal solution for dual & perturbation

    DblVec  rc;                 /// Reduced cost vector (for current dual & perturbation)
    bool    isFeasible;	        /// Is primal solution feasible?
    double  ubCost;		          /// Upper bound
    double  lbCost;		          /// Lower bound
    double  bestLbCost;         /// Best LB found so far
    double  bestUbCost;         /// Best LB found so far
    int     bestIter;           /// Best iteration index
    double  subgradFactorMin, subgradFactorMult;

    Random * rand;

    /// standard constructor with given size of problem
    Particle (Index _idx, double cost, qol::Model* m, qol::ModelType _t, 
      int nParticles, int seed, double factorMin, 
      double factorMult, Log* _l): idx(_idx), 
      maxCost (cost), model (m), type(_t), log(_l) {

        initialise (model->nVar(), model->nConstr(), nParticles);
        bestLbCost = -inf;
        bestUbCost = inf;
        rand = new Random ();
        rand->seed (seed);
        perturbFactor = 0.0;
        lbCheckFreq = 1;
        bestIter = -1;
        subgradFactorMin = factorMin;
        subgradFactorMult = factorMult;
    }

    virtual ~Particle() {}

    void fixConstraint (int constraint, SparseVec & feas, DblVec x);

    void perturbationDirection (double eps, DblVec dualLB, DblVec dualUB, 
      DblVec x);

    void updateVelocity (double eps, double subgradFactor, double globalFactor,
      double bestUB, double velFactor, DblVec bestDual, 
      DblVec dualLB, DblVec dualUB, DblVec bestPerturb, 
      DblVec bestx);
    void makeStep (DblVec dualLB, DblVec dualUB);

    ///Main method
    qol::Status solve (int lbcheckfreq, double & factor, int iter=-1);
    qol::Status heuristic (int iter, int heurFreq);

    double getLB () const { return lbCost; }
    double getUB () const { return ubCost; }

    int getIdx () const { return int(idx); }

    double getBestLB () const { return bestLbCost; }
    void setBestLB (int _best) { bestLbCost = _best; }

    DblVec getReducedCost () const { return rc; }

    void setpVel (DblVec _vel) { pVel = _vel; }
    DblVec getpVel () const { return pVel; }

    void setPerturVel (DblVec _vel) { pVel = _vel; }
    DblVec getPerturVel () const { return pVel; }

    void setDualVel (DblVec _vel) { dVel = _vel; }
    DblVec getDualVel () const { return dVel; }

    void setDual (DblVec _d) { dual = _d; }
    DblVec getDual () const { return dual; }

    void setPerturb (DblVec _p) { perturb = _p; }
    DblVec getPerturb () const { return perturb; }

    DblVec getPerturbDir () const { return perturbDir; }

    void setPerturbFactor (double _p) { perturbFactor = _p; }
    double getPerturbFactor () const { return perturbFactor; }

    DblVec getViolations () const { return viol; }
    void setViolations (DblVec _viol) { viol = _viol; }

    void setLBCheckFreq (int _freq) { lbCheckFreq = _freq; }

    void setBestIter (int _best) { bestIter = _best; }
    int getBestIter () const { return bestIter; }

    bool getIsFeasible () const { return isFeasible; }
  };

  class PSO: public Optimiser {
  protected:
    double maxCost;       //for lagrangian, the maximum lagrangian cost
    DblVec dualUB, dualLB;

  private:
    std::vector<Particle *> swarm;
    std::vector<Particle *>::iterator swarmIt;
    double bestLB, bestUB;
    DblVec bestDual, bestViol, bestPerturb, bestx;

    double absGap, relGap, globalFactor, velFactor, perturbFactorVal,
      subgradFactorVal, subgradFactorMin, subgradFactorMult;
    int lbcheckfreq, nParticles, heurFreq;

  public:
    PSO (qol::ModelType type, qol::FormulationType form): 
        Optimiser(type, form) {
          swarm.clear();
          bestUB = inf;
          bestLB = -inf;
          dualUB.clear();
          dualLB.clear();
          bestDual.clear();
          bestViol.clear();
          bestPerturb.clear();
          bestx.clear();
          log = new qol::Log(3, "logfile.log", "");

          bestDual.resize (model->nConstr(), inf);
          bestViol.resize (model->nConstr(), 0.0);
          bestPerturb.resize (model->nVar(), 0.0);
          bestx.resize (model->nVar(), inf);
          dualUB.resize (model->nConstr(), inf);
          dualLB.resize (model->nConstr(), -inf);

          absGap = relGap = globalFactor = velFactor = 0.0;
          perturbFactorVal = subgradFactorVal = subgradFactorMin = 0.0;
          subgradFactorMult = 0.0;
          lbcheckfreq = nParticles = 0, heurFreq = 0;
          status = qol::FAILED;

          params.addParam (qol::ABSGAP, "abs-gap", 0.0, 0.0, 1.0, 0.0001);
          params.addParam (qol::RELGAP, "rel-gap", 0.0, 0.0, 1.0, 0.999);
          params.addParam (qol::VERBOSITY, "verbosity", 0, 0, 3, 1);
          params.addParam (qol::PRINTFREQ, "print-freq", 0, 0, 3, 1);
          params.addParam (qol::ITERATIONS, "iterations", 0, 0, 10000, 1000);
          params.addParam (qol::EPS, "eps", 0, 0, 0.1, 1.e-16);

          params.addParam (qol::SUBGRADFACTOR, "subgradfactor", 0.0, 0.0, 1.0, 0.4);
          params.addParam (qol::VELOCITY, "velocityfactor", 0.0, 0.0, 1.0, 0.1);
          params.addParam (qol::GLOBALFACTOR, "globalfactor", 0.0, 0.0, 1.0, 0.1);
          params.addParam (qol::HEURISTIC_FREQ, "heur-freq", 0, 0, 1, 1);
          params.addParam (qol::LBCHECKFREQ, "lbcheckfreq", 0, 0, 100, 10);
          params.addParam (qol::PERTURBFACTOR, "perturbfactor", 0.0, 0.0, 1.0, 1e-3);
          params.addParam (qol::SUBGRADFMIN, "subgradfmin", 0.0, 0.0, 1.0, 0.01);
          params.addParam (qol::SUBGRADFMULT, "subgradfmult", 0.0, 0.0, 1.0, 0.8);

          params.addParam (qol::NPARTICLES, "CPU Number", 2, 0, 20, omp_get_max_threads());
          params.addParam (qol::CPU_NUMBER, "nparticles", 5, 0, 20, omp_get_max_threads());
        }    

        virtual void initialise(int argc,char **argv);

        virtual ~PSO() {
          if (log) delete (log);
        }

        void setViolations (DblVec viol);

  protected:
    /// interaction with solver
    virtual qol::Status AlgorithmSolve () ; // do heuristic/primal solving

  private:
    void updateBest ();
    void printIterStatus (int iter, Particle *p);
  };
}
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
* Local Variables:
* tab-width: 4
* eval: (c-set-style "stroustrup")
* End:
*/
