#include "QolCommData.h"

namespace qol {
  CommData::CommData () {
    //Initialize variables
    me = nprocs = -1;
    if (vars) vars->clear();
    if (cnstrs) cnstrs->clear();

    //Initialize MPI
//    MPI_Comm_rank(MPI_COMM_WORLD,&me);
//    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  }

  CommData::~CommData () {
    if (vars) delete(vars);
    if (cnstrs) delete(cnstrs);
  }

  void CommData::send_variables () {
  }

  void CommData::send_constraints() {
  }

  void CommData::recv_variables () {
    std::vector <Variable> recv_vars;
    //problem->setvariables(recv_vars);
    //vars = problem->variables;
  }

  void CommData::recv_constraints() {
    std::vector <Constraint> recv_cnstr;
    //problem->setconstraints(recv_cnstr);
    //cnstr = problem->constraints;
  }

  void CommData::perform_swap_problem () {
    ///TODO
  }
}
