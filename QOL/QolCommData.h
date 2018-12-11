#ifndef __QOL_COMMDATA_H__
#define __QOL_COMMDATA_H__

//#include "mpi.h"
#include "QolData.h"
#include "QolUtil.h"
#include <vector>

namespace qol {
  class CommData {

  private:
    int me,nprocs;

    //not sure about how to put this
    std:: vector <Variable> *vars;
    std::vector <Constraint> *cnstrs;

    void send_variables ();
    void send_constraints();

    void recv_variables ();
    void recv_constraints();

  public:
    CommData ();
    ~CommData ();

    //Sends and receives variables/constraints... something. Not sure yet.
    //It uses the private methods of this class.
    void perform_swap_problem ();
  };
}
#endif
