// WARNING: this file adapted for comparison purposes from the ilobendersatsp.cpp
// file that comes with CPLEX so copyright belongs to IBM.

// Example ilobendersatsp.cpp solves a flow MILP model for an
// Asymmetric Traveling Salesman Problem (ATSP) instance
// through Benders decomposition.
//
// The arc costs of an ATSP instance are read from an input file.
// The flow MILP model is decomposed into a master ILP and a worker LP.
//
// The master ILP is then solved by adding Benders' cuts during
// the branch-and-cut process via the cut callback functions.
// The cut callback functions add to the master ILP violated Benders' cuts
// that are found by solving the worker LP.
//
//
// To run this example, command line arguments are required:
//     ilobendersatsp.cpp [filename]
//
//     filename  Is the name of the file containing the ATSP instance (arc costs).
//               If filename is not specified, the instance
//               ../../../examples/data/atsp.dat is read
//
//
// ATSP instance defined on a directed graph G = (V, A)
// - V = {0, ..., n-1}, V0 = V \ {0}
// - A = {(i,j) : i in V, j in V, i != j }
// - forall i in V: delta+(i) = {(i,j) in A : j in V}
// - forall i in V: delta-(i) = {(j,i) in A : j in V}
// - c(i,j) = traveling cost associated with (i,j) in A
//
// Flow MILP model
//
// Modeling variables:
// forall (i,j) in A:
//    x(i,j) = 1, if arc (i,j) is selected
//           = 0, otherwise
// forall k in V0, forall (i,j) in A:
//    y(k,i,j) = flow of the commodity k through arc (i,j)
//
// Objective:
// minimize sum((i,j) in A) c(i,j) * x(i,j)
//
// Degree constraints:
// forall i in V: sum((i,j) in delta+(i)) x(i,j) = 1
// forall i in V: sum((j,i) in delta-(i)) x(j,i) = 1
//
// Binary constraints on arc variables:
// forall (i,j) in A: x(i,j) in {0, 1}
//
// Flow constraints:
// forall k in V0, forall i in V:
//    sum((i,j) in delta+(i)) y(k,i,j) - sum((j,i) in delta-(i)) y(k,j,i) = q(k,i)
//    where q(k,i) =  1, if i = 0
//                 = -1, if k == i
//                 =  0, otherwise
//
// Capacity constraints:
// forall k in V0, for all (i,j) in A: y(k,i,j) <= x(i,j)
//
// Nonnegativity of flow variables:
// forall k in V0, for all (i,j) in A: y(k,i,j) >= 0
//

#include "CplexFormulation.h"
#include "GurobiFormulation.h"
#include "QolUtil.h"
#include <string>
#include <fstream>

typedef std::vector<qol::DblVec> DblMatrix; // vector<vector<double>>
typedef std::vector<std::vector<qol::Variable> > VarMatrix;

// Implementation class for the user-defined lazy constraint callback.
// The function BendersLazyCallback allows to add Benders' cuts as lazy constraints.
class BendersCallback : public qol::Callback
{
	qol::MIPSolver *DPptr;
	VarMatrix u;		 // dual variables (flow constraints)
	std::vector<VarMatrix> v; // dual vars, capacity constraints
	const VarMatrix &x;			  // variables from the master
	const size_t numNodes;
public:
	bool separateFracSols;
	BendersCallback(qol::FormulationType solver,
					size_t n,const VarMatrix &xMaster) :
		numNodes(n),x(xMaster),separateFracSols(true)
		{ createWorkerLP(solver); }
	~BendersCallback() { delete DPptr; }
	void createWorkerLP(qol::FormulationType solver);
	Status callback(Progress where);
};


// This routine creates the master ILP (arc variables x and degree constraints).
//
// Modeling variables:
// forall (i,j) in A:
//    x(i,j) = 1, if arc (i,j) is selected
//           = 0, otherwise
//
// Objective:
// minimize sum((i,j) in A) c(i,j) * x(i,j)
//
// Degree constraints:
// forall i in V: sum((i,j) in delta+(i)) x(i,j) = 1
// forall i in V: sum((j,i) in delta-(i)) x(j,i) = 1
//
// Binary constraints on arc variables:
// forall (i,j) in A: x(i,j) in {0, 1}
//
void createMasterILP(qol::MIP &mip,VarMatrix &x,const DblMatrix &arcCost)
{
   size_t i, j;
   size_t numNodes = arcCost.size();
   char varName[100];
   // Create variables x(i,j) for (i,j) in A 
   // For simplicity, also dummy variables x(i,i) are created.
   // Those variables are fixed to 0 and do not partecipate to 
   // the constraints.
   x.reserve(numNodes);
   for (i = 0; i < numNodes; ++i) {
	   x.push_back(mip.addVarVec(numNodes));
	   for (j = 0; j < numNodes; ++j) {
		   sprintf(varName, "x.%d.%d", (int) i, (int) j);
		   mip.setVarType(x[i][j],qol::Variable::BINARY);
		   mip.setVarName(x[i][j],varName);
		   mip.setVarLB(x[i][i],0); 
		   if(j!=i) mip.setVarUB(x[i][j],1);
		   else mip.setVarUB(x[i][i],0); 
		   // Create objective function: minimize sum((i,j) in A ) c(i,j) * x(i,j)
		   mip.setObjCoeff(x[i][j],arcCost[i][j]);
      }
   }

   // Add the out degree constraints.
   // forall i in V: sum((i,j) in delta+(i)) x(i,j) = 1

   for (i = 0; i < numNodes; ++i) {
	   qol::Expression expr;
	   for (j = 0;   j < numNodes; ++j)  if(i!=j) expr += x[i][j];
	   mip.addConstraint(expr == 1);
   }

   // Add the in degree constraints.
   // forall i in V: sum((j,i) in delta-(i)) x(j,i) = 1
   for (i = 0; i < numNodes; i++) {
	   qol::Expression expr;
	   for (j = 0;   j < numNodes; j++) if(j!=i) expr += x[j][i];
	   mip.addConstraint(expr == 1);
   }

}// END createMasterILP

// This routine set up the formulation algorithm to solve the worker LP, and
// creates the worker LP (i.e., the dual of flow constraints and
// capacity constraints of the flow MILP)
//
// Modeling variables:
// forall k in V0, i in V:
//    u(k,i) = dual variable associated with flow constraint (k,i)
//
// forall k in V0, forall (i,j) in A:
//    v(k,i,j) = dual variable associated with capacity constraint (k,i,j)
//
// Objective:
// minimize sum(k in V0) sum((i,j) in A) x(i,j) * v(k,i,j)
//          - sum(k in V0) u(k,0) + sum(k in V0) u(k,k)
//
// Constraints:
// forall k in V0, forall (i,j) in A: u(k,i) - u(k,j) <= v(k,i,j)
//
// Nonnegativity on variables v(k,i,j)
// forall k in V0, forall (i,j) in A: v(k,i,j) >= 0
//
void BendersCallback::createWorkerLP(qol::FormulationType solver)
{

    size_t i, j, k;
	if(solver == qol::CPLEX){
		qol::CplexFormulation *cpx = new qol::CplexFormulation;
		DPptr = cpx;  
		// Turn off the presolve reductions and set the CPLEX optimizer
		// to solve the worker LP with primal simplex method.
		// VERY IMPORTANT!!!!
		// if we don't turn this off the presolve might deduce unboundedness and
		// not give us an extremal ray!
		CPXsetintparam(cpx->env,CPX_PARAM_REDUCE,0); // VERY IMPORTANT !!!!
		
		CPXsetintparam(cpx->env,CPX_PARAM_LPMETHOD,CPX_ALG_PRIMAL);
		//cplex.setParam(IloCplex::Reduce, 0);
		//cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);
	}else if(solver == qol::GUROBI){
		qol::GurobiFormulation *grb = new qol::GurobiFormulation;
		DPptr = grb;  
		GRBsetintparam(grb->env,GRB_INT_PAR_DUALREDUCTIONS,0);
		// need the following parameter to get unbounded ray information
		GRBsetintparam(grb->env,GRB_INT_PAR_INFUNBDINFO,1);
		GRBsetintparam(grb->env,GRB_INT_PAR_METHOD,GRB_METHOD_PRIMAL);
		// don't show information about sub-problems
		GRBsetintparam(grb->env,GRB_INT_PAR_LOGTOCONSOLE,0);
	}
	qol::MIPSolver &DP=*DPptr; // dual problem
   
	// Create variables v(k,i,j) forall k in V0, (i,j) in A
	// For simplicity, also dummy variables v(k,i,i) are created.
	// Those variables are fixed to 0 and do not partecipate to 
	// the constraints.
	v.resize(numNodes);
	for (k = 1; k < numNodes; ++k) {
		v[k].resize(numNodes);
		for (i = 0; i < numNodes; ++i) {
			v[k][i] = DP.addVarVec(numNodes);
			for(j = 0; j < numNodes; ++j) {
				char varName[100];
				sprintf(varName, "v.%d.%d.%d", (int) k, (int) i, (int) j); 
				DP.setVar(v[k][i][j],0,qol::inf,0,
						  qol::Variable::CONTINUOUS,varName);
			}
			DP.setVarUB(v[k][i][i],0);
		}
	}
   

	// Create variables u(k,i) forall k in V0, i in V
	u.resize(numNodes);
	// Set names for variables u(k,i)
	for (k = 1; k < numNodes; ++k) {
		u[k] = DP.addVarVec(numNodes);
		for(i = 0; i < numNodes; ++i) {
			char varName[100];
			sprintf(varName, "u.%d.%d", (int) k, (int) i);
			DP.setVar(u[k][i],-qol::inf,qol::inf,0.0,
					  qol::Variable::CONTINUOUS,varName);
		}
	}

	// Initial objective function is empty
	//obj.setSense(IloObjective::Minimize);

	// Add constraints:
	// forall k in V0, forall (i,j) in A: u(k,i) - u(k,j) <= v(k,i,j)

	for (k = 1; k < numNodes; ++k) {
		for(i = 0; i < numNodes; ++i) {
			for(j = 0; j < numNodes; ++j) {
				if ( i != j ) {
					DP.addConstraint(u[k][i] - u[k][j] - v[k][i][j] <= 0);
				}
			}
		}
	}

}// END createWorkerLP



// This routine separates Benders' cuts violated by the current x solution.
// Violated cuts are found by solving the worker LP
//
qol::Callback::Status BendersCallback::callback(Progress where)
{
	std::cout << "callback where " << where << ": LB="
			  << getBound() << " - UB=" << getBest()
			  << std::endl;
	switch(where){		
		case PRESOLVE: case SIMPLEX: return OK;
		default:
			if( ! separateFracSols && where != MIP_SOL)
				return OK;		// only do integer cuts
			break;
	}
   size_t i, j, k;
   qol::MIPSolver &DP=*DPptr; // dual problem

   // Update the objective function in the worker LP:
   // minimize sum(k in V0) sum((i,j) in A) x(i,j) * v(k,i,j)
   //          - sum(k in V0) u(k,0) + sum(k in V0) u(k,k)
   
   for (k = 1; k < numNodes; ++k) {
      for (i = 0; i < numNodes; ++i) {
		  for (j = 0; j < numNodes; ++j) {
			  DP.setObjCoeff(v[k][i][j],primalSolution[x[i][j]]);
		  }
      }
   }
   for (k = 1; k < numNodes; ++k) {
	   DP.setObjCoeff(u[k][k],1);
	   DP.setObjCoeff(u[k][0],-1);
   }

   // Solve the worker LP
   qol::Status  status = DP.solveRelaxed(); // solve as LP

   // A violated cut is available iff the solution status is Unbounded
   if ( status == qol::UNBOUNDED ) {
	   std::cout << "\tUnbounded - ";
	   qol::DblVec ray(DP.nVar(),-1);
	   // Get the violated cut as an unbounded ray of the worker LP
	   DP.getRay(ray);
	   double r=0;
	   for(size_t ri=0;ri<ray.size();++ri) r+=ray[ri]*ray[ri];
	   std::cout << "Ray length " << ray.size()
				 << " expect " << (numNodes-1)*numNodes*(1+numNodes)
				 << " norm^2 = " << r << std::endl;
      // Compute the cut from the unbounded ray. The cut is:
      // sum((i,j) in A) (sum(k in V0) v(k,i,j)) * x(i,j) >=
      // sum(k in V0) u(k,0) - u(k,k)
	   qol::Expression cutLHS;
	   double lhsVal=0;
	   for(i=0;i<numNodes;++i)
		   for(j=0;j<numNodes;++j){
			   double sum=0;	
			   for(k=1;k<numNodes;++k)	
				   sum += ray[v[k][i][j]];
			   if(sum != 0.0)
				   std::cout << "+ " << sum  << "*x"<<i<<","<<j<<" ";
			   lhsVal += sum*primalSolution[x[i][j]];
			   cutLHS += sum*x[i][j];
		   }
	   double cutRHS=0;
	   for(k=1;k<numNodes;++k)
		   cutRHS += ray[u[k][0]] - ray[u[k][k]];
	   std::cout << ">= " << cutRHS
				 << "  viol. = " << (cutRHS-lhsVal)
				 << std::endl;
	   addCut( cutLHS >= cutRHS);
   }else{
	   std::cout << "\terror status: " << DP.getErrorMessage() << std::endl;
	   std::cout << "\tsolution status: " << status << std::endl;
	   DP.writeLP("bounded.lp");
	   std::cout << "\tWrote bounded.lp\n";
   }

   return OK;

} // END separate

void usage (char *progname)
{  using namespace std;
   cerr << "Usage:     " << progname << " {0|1} {c|g} [filename]"    << endl;
   cerr << " 0:        Benders' cuts only used as lazy constraints," << endl;
   cerr << "           to separate integer infeasible solutions."    << endl;
   cerr << " 1:        Benders' cuts also used as user cuts,"        << endl;
   cerr << "           to separate fractional infeasible solutions." << endl;
   cerr << " c:        Use CPLEX as solver (default)" 	 			 << endl;
   cerr << " g:        Use Gurobi as solver"  			 			 << endl;
   cerr << " filename: ATSP instance file name (default atsp.dat)."  << endl;
   exit(1);
} // END usage

int main(int argc, char **argv)
{  using namespace std;
   qol::MIPSolver *masterPtr = 0;
   try {
      const char* fileName = "atsp.dat";

	  // Check the command line arguments
	  qol::FormulationType solver = qol::CPLEX;
	  bool separateFracSols=true;		
	  for(int a=1;a<argc;++a)
		  if( argv[a][1] == '\0'){
			  switch(argv[a][0]){
				  case 'c': case 'C':
					  solver = qol::CPLEX;
					  masterPtr = new qol::CplexFormulation;
					  break;
				  case 'g': case 'G':
					  solver = qol::GUROBI;
					  masterPtr = new qol::GurobiFormulation;
					  break;
				  case '0': separateFracSols = false; break;
				  case '1': separateFracSols = true;  break;
				  default:
					  usage(argv[0]);
			  }
		  }else
			  fileName = argv[a];

	  std::cout << "Benders' cuts separated to cut off: ";
      if ( separateFracSols ) {
         std::cout << "Integer and fractional infeasible solutions." << endl;
      } else {
         std::cout << "Only integer infeasible solutions." << endl;
      }
	  std::cout << "Reading data from " << fileName << std::endl;

      // Read arc_costs from data file (17 city problem)
	  DblMatrix arcCost;
	  size_t numNodes;
	  ifstream data(fileName);
      if ( !data ) throw(-1);
	  data >> numNodes;
	  arcCost.resize(numNodes);
	  for(size_t i=0;i<numNodes;++i){
		  arcCost[i].resize(numNodes);
		  for(size_t j=0;j<numNodes;++j)
			  data >> arcCost[i][j];
	  }
      data.close();

      // create master ILP
	  if(masterPtr == 0) masterPtr = new qol::CplexFormulation;	  
	  qol::MIPSolver &master = *masterPtr; 
	  
      VarMatrix x;		  
      createMasterILP(master, x,arcCost);

      // Create worker IloCplex algorithm and worker LP for Benders' cuts separation
	  BendersCallback benders(solver,numNodes,x);  // calls createWorkerLP();
	  benders.separateFracSols = separateFracSols;
	  
      // Set up the cut callback to be used for separating Benders' cuts
      master.setCallback(&benders);
	  if(solver == qol::CPLEX){
		  CPXENVptr env=(dynamic_cast<qol::CplexFormulation *>(masterPtr))->env;
		  CPXsetintparam(env,CPX_PARAM_PREIND,0);
		  //masterCplex.setParam(IloCplex::PreInd, IloFalse); 
		  // Set the maximum number of threads to 1. 
		  // This instruction is redundant: If MIP control callbacks are registered, 
		  // then by default CPLEX uses 1 (one) thread only.
		  // Note that the current example may not work properly if more than 1 threads 
		  // are used, because the callback functions modify shared global data.
		  // We refer the user to the documentation to see how to deal with multi-thread 
		  // runs in presence of MIP control callbacks. 

		  //masterCplex.setParam(IloCplex::Threads, 1); 
		  CPXsetintparam(env,CPX_PARAM_THREADS,1);
	  
		  CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
		  // Turn on traditional search for use with control callbacks
		  //masterCplex.setParam(IloCplex::MIPSearch, IloCplex::Traditional);
		  CPXsetintparam(env,CPX_PARAM_MIPSEARCH,CPX_MIPSEARCH_TRADITIONAL);
	  }else if(solver == qol::GUROBI){					// gurobi
		  GRBenv *env = (dynamic_cast<qol::GurobiFormulation *>(masterPtr))->env;
		  GRBsetintparam(env, GRB_INT_PAR_PRESOLVE,0);
		  // precrush=1 means cuts are transformed based on presolve reductions
		  GRBsetintparam(env, GRB_INT_PAR_PRECRUSH,1); 
		  // Need to specifically enable lazy constraints
		  // (not sure if this applies is required for cuts as well?)
		  GRBsetintparam(env, GRB_INT_PAR_LAZYCONSTRAINTS,1);
		  GRBsetintparam(env, GRB_INT_PAR_LOGTOCONSOLE,1);
		  GRBsetintparam(env, GRB_INT_PAR_THREADS, 1); // is this necessary?
	  }else{
		  std::cout << "ERROR: solver type " << solver << " not handled\n";
		  return -1;
	  }
      //if ( separateFracSols )
      //   masterCplex.use(BendersUserCallback(masterEnv, x, workerCplex, v, u, workerObj));

      // Solve the model and write out the solution
	  qol::Status status = master.solveExact() ;
	  switch(status ) {
		  case qol::FAILED:  case qol::INFEASIBLE: case qol::UNBOUNDED:
			  std::cout << "Problem could not be solved: " << status
						<< endl;
			  break;
		  default:
			  std::cout << endl << "Solution status: " << status << endl;
			  std::cout << "Objective value: "
						<< master.getObjective() << endl;
	  }
	  master.writeLP("master.lp");
	  std::cout << "Wrote final master to master.lp\n";
	  if ( status == qol::OPTIMAL ) {
		  // Write out the optimal tour
		  size_t i,j;
		  std::vector<int> succ(numNodes,-1);
		  for (i = 0; i < numNodes; i++) {
			  for(j = 0; j < numNodes; j++) {
                  if ( master.getPrimal(x[i][j]) > 1e-03 ){
					  succ[i] = j;
					  break;
				  }
               }
            }
            std::cout << "Optimal tour:" << endl;
            i = 0;
            while ( succ[i] != 0 ) {
               std::cout << i << ", ";
               i = succ[i];
            }
            std::cout << i << endl;
	  } else { 
		  std::cout << "Solution status is not Optimal" << endl;
		  std::cout << "No solution available" << endl;
      }
   } catch (qol::Exception error){
	   cerr << "QOL error: " << error.what() << std::endl;
	   return 2;
   } catch (std::exception error){
	   cerr << "std error: " << error.what() << std::endl;
	   return 3;
   } catch (...) {
      cerr << "Unknown exception caught!" << endl;
   }
   if(masterPtr) delete masterPtr;
   return 0;

} // END main
