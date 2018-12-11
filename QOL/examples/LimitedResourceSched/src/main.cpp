#include <iostream>
#include <sstream>

#ifdef OPENMP
#include <omp.h>
#else
#define  omp_get_max_threads() 1
#endif

#include "data.h"

#include "PSO.h"
#include "QolColFormulation.h"
#include "GurobiFormulation.h"
#ifdef _USE_GUROBI_
#  include <gurobi_c++.h> // for exception definition
#endif
#include "CplexFormulation.h"


qol::AlgorithmType alg;
qol::FormulationType form;
qol::ModelType modtype;


void usage(char **argv) {
  std::cerr << "USAGE: " << argv[0] << " "
    << "[--alg val] { [--params ] [-param val] } [--model mod] [--solver solv] <input_file>\n"
    << "Solve Resource Constrained Scheduling problems \n"
    << "Arguments give file text data format (As for ACO algorithm)\n"
    << "Optionally set parameters (eg --maxIter 100 --maxCPU 600). \n"
    << "The optional argument 'vol' or 'v' activates the volume\n"
    << "algorithm instead of lagragian particle swarm optimisation"
    << std::endl;
  exit(1);
}

void parse_args(int argc, char **argv) { 
    int currentarg = 0;

    alg = qol::QOLPSO;
    modtype = qol::LAGRANGIAN;
    form = qol::CPLEX;

    while (currentarg<argc) {
      if (!strcmp (argv[currentarg], "--alg"))
      {
        currentarg++;
        if (!strcmp (argv[currentarg], "PSO")) {
          alg = qol::QOLPSO;          
        }
        if (!strcmp (argv[currentarg], "ACO")) {
          alg = qol::QOLACO;
        }
      }
      if (!strcmp (argv[currentarg], "--model")) {
         currentarg++;
        if (!strcmp (argv[currentarg], "LAG")) {
          modtype = qol::LAGRANGIAN;
        }
        if (!strcmp (argv[currentarg], "BENDERS")) {
          modtype = qol::BENDERS;
        }
      }
      if (!strcmp (argv[currentarg], "--solver")) {
         currentarg++;
        if (!strcmp (argv[currentarg], "CPLEX")) {
          form = qol::CPLEX;
        }
        if (!strcmp (argv[currentarg], "GUROBI")) {
          form = qol::GUROBI;
        }
      }
      currentarg++;
    }
}

void createModel (DATA data,  qol::Optimiser *optimiser) {
  qol::QolColFormulation *model = new qol::QolColFormulation();
  std::vector < std::vector <qol::Variable> > xjt (data.total_jobs), zjt(data.total_jobs);

  /// put all the variables in the model (general)

  for(int j=0; j < data.total_jobs; ++j){
    xjt[j] = model->addVarVec (data.max_time);
    zjt[j] = model->addVarVec (data.max_time);
    for(int t=0; t < data.max_time; t++){
      stringstream x_name, z_name;
      x_name << "x(" << j << ")(" << t << ")";
      z_name << "z(" << j << ")(" << t << ")";

      model->setVarName (xjt[j][t], x_name.str());
      model->setVarName (zjt[j][t], z_name.str());

      model->setVarLB (xjt[j][t], 0.0);
      model->setVarLB (zjt[j][t], 0.0);

      model->setVarUB (xjt[j][t], 1.0);
      model->setVarUB (zjt[j][t], 1.0);

      model->setObjCoeff (xjt[j][t], 0.0);
      model->setObjCoeff (zjt[j][t], 0.0);

      model->setVarType (xjt[j][t], qol::Variable::INTEGER);
      model->setVarType (zjt[j][t], qol::Variable::INTEGER);

      if(t < data.rel[j]+data.dur[j]){
        model->setVarUB (xjt[j][t], 0.0);
        model->setVarUB (zjt[j][t], 0.0);
      }
    }
    // job must get completed by end of the horizon
    model->setVarLB (zjt[j][data.max_time-1], 1.0);
  }

  qol::Expression obj;
  qol::Index currentProb = 0;

  qol::MIPSolver *problem =
# ifdef _USE_GUROBI_
	  new qol::GurobiFormulation();
# else
	  new qol::CplexFormulation ();
# endif
  for (int m=0; m<data.machines; ++m) {
    qol::MIPSolver *problem;

    if (optimiser->getFormulationType()==qol::GUROBI) {
#ifdef _USE_GUROBI_
      try {
        problem = new qol::GurobiFormulation ();
      }
      catch(GRBException e) {
		// can't do this unless compiling agains C++ libraries
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
        exit (-1);
      }
#endif
    }
    else {
#ifdef _USE_CPLEX_
      if (optimiser->getFormulationType()==qol::CPLEX) {
        problem = new qol::CplexFormulation ();
      }
#endif
    }

    qol::Constraint temp;

    std::vector < std::vector <qol::Variable> > subxjt (data.jobs_in_machines[m]), subzjt(data.jobs_in_machines[m]);
    problem->setVarSize (data.max_time);
    for (int j = data.first_job_in_machine[m]; j<data.first_job_in_machine[m]+data.jobs_in_machines[m]; ++j){
      int idx = j - data.first_job_in_machine[m];
      subxjt [idx] = problem->addVarVec (data.max_time);
      subzjt [idx] = problem->addVarVec (data.max_time);
      for (int t = 0; t<data.max_time; ++t) {
        ///extract the variables from the model and include them in the subproblems
        qol::Variable var = model->getVar (xjt[j][t]);
        qol::Index temp = qol::Index (var);
        optimiser->getModel()->addVarAtProblem (var, currentProb);
      //  model-> addVarAtProblem (var, optimiser->nSubprob()-1);
        problem->setVarLB (subxjt[idx][t], model->getVarLB (var));
        problem->setVarUB (subxjt[idx][t], model->getVarUB (var));
        problem->setVarName (subxjt[idx][t], model->getVarName (var));
        problem->setObjCoeff (subxjt[idx][t], model->getObjCoeff (var));
      }
      for (int t = 0; t<data.max_time; ++t) {
        qol::Variable var = model->getVar (zjt[j][t]);
        optimiser->getModel()->addVarAtProblem (var, currentProb);
        problem->setVarLB (subzjt[idx][t], model->getVarLB (var));
        problem->setVarUB (subzjt[idx][t], model->getVarUB (var));
        problem->setVarName (subzjt[idx][t], model->getVarName (var));
        problem->setObjCoeff (subzjt[idx][t], model->getObjCoeff (var));
      }
    }

    for (int t = 0; t<data.max_time; ++t) {
      qol::Expression sum1;
      for (int j = data.first_job_in_machine[m]; j<data.first_job_in_machine[m]+data.jobs_in_machines[m]; ++j){
        int idx = j - data.first_job_in_machine[m];
        if (t>0) {
          problem->addConstraint (subzjt[idx][t] - subzjt[idx][t-1] >= 0);// once the job is completed stays completed
          problem->addConstraint (subxjt[idx][t] - (subzjt[idx][t]-subzjt[idx][t-1]) == 0);// relation between x and z variables
        }
        int t1 = std::min(t+data.dur[j], data.max_time-1);
        sum1+=subzjt[idx][t1];
        sum1-=subzjt[idx][t];
      }

      problem->addConstraint (sum1<=1);

      for(int d=0; d < data.precedences.size(); d++){
        int pred = data.precedences[d][0]-1;
        int succ = data.precedences[d][1]-1;
        if(data.job_belongs_to_machine[pred] != m) continue;
        if(data.job_belongs_to_machine[succ] != m) continue;
        if(t>=data.dur[succ]) {
          problem->addConstraint (subzjt[succ-data.first_job_in_machine[m]][t] - subzjt[pred-data.first_job_in_machine[m]][t-data.dur[succ]]<=0); // precedence relations
        }
      }
    }

    int i = data.first_job_in_machine[m]*data.max_time;
    for(int j=data.first_job_in_machine[m]; j < data.first_job_in_machine[m]+data.jobs_in_machines[m]; j++){
      int idx = j-data.first_job_in_machine[m];
      for(int t=0; t < data.max_time; t++){
        problem->setObjCoeff (subxjt[idx][t], data.weight[j]); ///obj += 1.0*xjt[j][t];
        //problem->setObjCoeff (subxjt[idx][t], 1.0); ///obj += 1.0*xjt[j][t];
        i++;
        if(t < data.rel[j]+data.dur[j]){
          qol::Expression temp;
          temp = subxjt[idx][t];
          problem->addConstraint (temp == 0); // job cannot be completed before rel+proc time;
          temp = subzjt[idx][t];
          problem->addConstraint (temp == 0);
        }
      }
    }

    std::vector <qol::ConstraintRow> lnkConstraints;
    for(int j=data.first_job_in_machine[m]; j < data.first_job_in_machine[m]+data.jobs_in_machines[m]; j++){
      int idx = j-data.first_job_in_machine[m];
      qol::ConstraintRow power;
      for(int t=0; t < data.max_time; t++){
        qol::Variable var = model->getVar (zjt[j][t]);
        qol::CoeffVarPair pair (1.0, var /*subzjt[idx][t]*/);
        power.push_back (pair);
      }
      power.rhs = data.power[j];
      power.sense = 'L';
      model->addConstraint (power);
      //problem->addConstraint (power);

      lnkConstraints.push_back (power);
    }
    ((qol::LagrangianModel *)optimiser->getModel())->addLnkConstraints (lnkConstraints);

    std::stringstream fname;
    fname << "out_" << m << ".lp";
    problem->writeLP (fname.str().c_str());
    optimiser->addProblem (problem);
    currentProb++;
  }
}


int main(int argc,char **argv) {
  if(argc < 2)
    usage(argv);

  parse_args (argc, argv);

  //this should go into parse_args
  const char *filename = argv[argc-1];
  if( ! fopen(filename,"r") ){
    std::cerr << "Cannot open " << filename << " for reading\n";
    usage(argv);
  }

  // override defaults with command line arguments

  DATA data;
  std::ifstream in(filename);

  data.ReadFile(in);

#ifdef _USE_CPLEX_
  qol::Optimiser * optimiser;
  if (alg == qol::QOLPSO) {
    optimiser = new qol::PSO(modtype, form);
  }
#else
  qol::PSO * optimiser = new qol::PSO(qol::LAGRANGIAN, qol::GUROBI);
#endif

  createModel(data, optimiser);
  optimiser->initialise (argc, argv);
  optimiser->solve();

  return 0;
}
