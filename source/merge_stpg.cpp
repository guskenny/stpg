#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <STPGSolver.h> 
//#include <math.h>

int main(int argc, const char** argv){
//  qol::MIPSolver *mipPtr=0; // uncomment for original


  if(argc <= 1){
    std::cout<<"Usage : ./merge_stpg [OPTION] PATH "<<std::endl;
    std::cout<<"       Options"<<std::endl;
    //std::cout<<"          -p    PCPSP"<<std::endl;
    std::cout<<"Sample "<<std::endl;
    std::cout<<"      ./mining -g Data/zuck_small"<<std::endl;
    std::cout<<"Options: g/c solve with gurobi/cplex\n";
    std::cout<<"-r solve relaxed LP only\n";
    exit(-1);
  }
  try{

    STPGSolver solver(argc,argv);
    solver.solve();

  } // end try statement
  catch (...) { std::cerr << "Unknown Error Occured" << std::endl;}
  PE("\n********* END OF PROGRAM *********")
  return 0;
}
