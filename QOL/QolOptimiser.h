#ifndef __QOL_OPTIMISER__
#define __QOL_OPTIMISER__

#include "CpuTimer.h"
#include "QolUtil.h"
#include "QolData.h"
#include "QolModel.h"
#include "QolColFormulation.h"
#include "LagrangianModel.h"
#include "VolumeModel.h"
#include "BendersModel.h"
#include "QolSolutionPrimal.h"
#include "QolMIPSolver.h"
#include "QolParams.h"
#include "QolLog.h"
#include <list>
#include <iostream>
#include <sstream>

namespace qol {

  /// Optimiser is the abstract base class for an optimisation
  /// method or non-trivial heuristic. The main idea is that any
  /// optimiser may be given some resources (number of threads or
  /// processes & CPU/elapsed time to run with) and will try to
  /// generate some solutions. Typcially an Optimiser will require
  /// a particular kind of model (eg MIPformulation, LagrangianModel)
  /// However that should be set-up in the constructor
  class Optimiser {
  protected:
    /// Parameters
    qol::Parameters params;
    qol::Status status;

    //flag to make sure the algorithm has been initialised
    bool initialised;

    qol::FormulationType form;
    qol::MIP * formulation;
    qol::Model * model;
    qol::ModelType type;
    qol::Log * log;

    // Timer is a class for timing CPU & wall usage
    // it also allows us to set limits on these
    CpuTimer timer;

    std::list<SolutionPrimal *> primal;     // primal feasible solutions, (best first)
    std::list<SolutionPrimal *> bound;	    // lower bound solutions

  public:
    Optimiser (qol::ModelType t, qol::FormulationType f){
      type = t;
      switch (type) {
        case BENDERS: model = new qol::BendersModel();
          break;
        case LAGRANGIAN: model = new qol::LagrangianModel();
          break;
/*        case VOLUME: model = new qol::VolumeModel;
          break;
*/
        default: return;
          break;
      }
      form = f;
      switch (form){
        case QOLFORM: formulation = new qol::QolColFormulation();
          break;
          /*******
          *****/
      }
      primal.clear();
      bound.clear();
      params.clear();
      timer.reset();
      initialised = false;
    }

    virtual ~Optimiser (){
      if (model) delete(model);
      if (formulation) delete(formulation);
    }

    virtual void initialise (int argc,char **argv) = 0;

    qol::FormulationType getFormulationType () const {return form;}

    virtual void setThreads(int nThreads) {params.addParam (qol::THREADS, "nthreads", 0, 0,32, nThreads); } // how many threads we may use

    virtual int getThreads() { return int (params.getParamValue (qol::THREADS)); };
    //virtual int getThreads() = 0;

    /// Params
    // Allow setting of optimiser specific parameters in a generic way
    Index numParam() const { return params.size(); } // how many are there

    // interaction with solver
    virtual qol::Status solve() {
      qol::Status status = FAILED;
      if (initialised) 
        status = AlgorithmSolve();
      else
        log->error ("The algorithm hasn't been initialised\n");
      initialised = false;
      return status;
    }

    virtual void writeLP(const std::string filename) ;

    void addProblem (qol::MIPSolver* p) {model->addProblem (p);}

    //Return the model. Used to interact directly with the model from main files.
    qol::Model * getModel () { return model; }
  
  protected:
    virtual qol::Status AlgorithmSolve () = 0; // do heuristic/primal solving
  };
}
#endif

/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
