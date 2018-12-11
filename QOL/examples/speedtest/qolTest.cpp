// mipTest.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"

#include <iostream>
#include <boost/format.hpp>
#include "QolMIP.h"
#include "QolColFormulation.h"
#include "CplexFormulation.h"
#include "GurobiFormulation.h"
#include "CpuTimer.h"

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>
#include <gurobi_c++.h>

//using namespace qol;

/// \code
/// using namespace qol
/// MIP::Formulation m;
/// MIP::Variable x(m),y(m),z(m);
/// x >= 0; y >= 0; z >= 0; // set lower bounds
/// x <= 1; y <= 2; z <= 1; // set upper bounds
/// x.setObjCoeff(2); y.setObjCoeff(-1.0);
/// m.addRow( x + 2*y == 2.5 )
/// Expression sum;
/// for(i=0;i<3;++i) sum += m.getVar(i)
/// m.addRow( sum <= 2 );
/// m.loadProblem();
/// \endcode
/// This will solve a problem of the form
/// \f[ \min 2 x - y \f]
/// Subject to
/// \f$  x+2\times y = 2.5\f$ \f[ \sum_{var\in\{x,y,z\}} var \le 2 \f]
/* the following doxygen markup only works in LaTeX not in html
  \f{eqnarray*}{
    x+2y &=& 2.5\\ \sum_{var\in\{x,y,z\}} var &\le& 2 
  \f}
*/


void createLargeMIP(qol::MIP &mip, int size)
// create problem with size*size variables and 2*size constraints
{
	qol::CpuTimer timer,totalTime;
	std::vector<std::vector<qol::Variable> > x;
	x.reserve(size);
	for(int i=0;i<size;++i){
		x.push_back(mip.addVarVec(size));
		for(qol::Index j=0;j<size;++j){
			mip.setVarLB(x[i][j],0.0);
			mip.setVarUB(x[i][j],1.0);
			mip.setVarName(x[i][j],boost::str(boost::format(
												  "x_%03d_%03d")%i%j));
			mip.setVarType(x[i][j],qol::Variable::BINARY);
			mip.setObjCoeff(x[i][j],i*j/double(size));
		}
	}
	std::cout << "Created " << mip.nVar() << " variables in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " wall-time" << std::endl;
	timer.reset();
	qol::CpuTimer tmr;
	for(int cnt=0;cnt < 10;++cnt)
		for(qol::Index i=0;i<size;++i){
			qol::Expression sum;
			for(qol::Index j=0;j<size;++j)
				sum += x[i][j];
			mip.addConstraint( sum == 1.0).setName("LHS_%d_%d",i,cnt);
		}
	std::cout << size << "\tsimple sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt)
		for(qol::Index j=0;j<size;++j) {
			qol::Expression sum;
			for(qol::Index i=0;i<size;++i)
				sum += 1.0*x[i][j];
			mip.addConstraint( sum == 1.0).setName("RHS_%d_%d",j,cnt);
		}
	std::cout << size << "\tcoeff. sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		qol::Expression sumAll;
		for(qol::Index i=0;i<size;++i){
			qol::Expression sum1,sum2;
			for(qol::Index j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(qol::Index j=size/3;j<size;++j)
				sum2 += x[i][j];
			sumAll = sum1 + sum2;
		}
		mip.addConstraint( sumAll >= cnt ).setName("all_%d",cnt);
	}
	std::cout << 10 << "\texpression sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";	
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		for(qol::Index i=0;i<size;++i){
			qol::Expression sum1,sum2;
			for(qol::Index j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(qol::Index j=size/3;j<size;++j)
				sum2 += x[i][j];
			mip.addConstraint( sum1 + 10 >= sum2+5 ).setName("lhs_rhs_%d",cnt);
		}
	}
	std::cout << size << "\texpression on both sides constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";	
	std::cout << "Created " << mip.nConstr() << " constraints in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	try{
		tmr.reset();
		qol::MIPSolver &solver = dynamic_cast<qol::MIPSolver&>(mip);
		solver.load();
		std::cout << "Loaded in " << tmr.elapsedSeconds() << " sec (CPU) "
				  << tmr.elapsedWallTime() << " sec wall "
				  << timer.elapsedSeconds() << " total CPU" << std::endl;
	} catch (...){
		std::cout << "Cannot load - column formulation only\n";
	}
		
 	timer.reset();
  	mip.writeLP("xxx.lp");
  	std::cout << "Wrote xxx.lp in "
  			  << timer.elapsedSeconds() << " sec (CPU) "
  			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	std::cout << "Total time: "
			  << totalTime.elapsedSeconds() << " sec (CPU) "
			  << totalTime.elapsedWallTime() << " sec wall-time" << std::endl;
}


void createLargeMIP_CPLEX(IloEnv &env,IloCplex &cplex,int size)
// create problem with size*size variables and 2*size constraints
{
	qol::CpuTimer timer,totalTime;
	IloModel model(env);
	IloArray<IloIntVarArray> x(env,size);
	IloExpr obj(env);
	for(int i=0;i<size;++i){
		x[i] = IloIntVarArray(env,size,0,1);
		for(IloInt j=0;j<size;++j){
			x[i][j].setLB(0);
			x[i][j].setUB(1);
			x[i][j].setName(boost::str(boost::format(
										   "x_%03d_%03d")%i%j).c_str());
			//mip.setVarType(x[i][j],qol::Variable::BINARY);
			obj += x[i][j] * (i*j/double(size));
		}
	}
	model.add(IloMinimize(env,obj));
	obj.end();
	std::cout << "Created " << size*size << " variables in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " wall-time" << std::endl;
	timer.reset();
	qol::CpuTimer tmr;
	for(int cnt=0;cnt < 10;++cnt){
		for(IloInt i=0;i<size;++i){
			IloExpr sum(env);
			for(IloInt j=0;j<size;++j)
				sum += x[i][j];
			model.add( sum == 1.0).setName(
				boost::str(boost::format("LHS_%d_%d") % i%cnt ).c_str());
		}
	}
	std::cout << size << "\tsimple sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt)	
		for(IloInt j=0;j<size;++j){
			IloExpr sum(env);
			for(IloInt i=0;i<size;++i)
				sum += 1.0*x[i][j];
			model.add( sum == 1.0).setName(
				boost::str(boost::format("RHS_%d_%d") % j%cnt ).c_str());
		}
		std::cout << size << "\tcoeff. sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		IloExpr sumAll(env);
		for(IloInt i=0;i<size;++i){
			IloExpr sum1(env),sum2(env);
			for(IloInt j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(IloInt j=size/3;j<size;++j)
				sum2 += x[i][j];
			sumAll += sum1 + sum2;
		}
		model.add( sumAll >= size).setName(
			boost::str(boost::format("all_%d") % cnt).c_str());
	}
	std::cout << 10 << "\texpr. sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		for(IloInt i=0;i<size;++i){
			IloExpr sum1(env),sum2(env);
			for(IloInt j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(IloInt j=size/3;j<size;++j)
				sum2 += x[i][j];
			model.add( sum1+10 >= sum2+5).setName(
				boost::str(boost::format("lhs_rhs_%d") % cnt).c_str());
		}
	}
	std::cout << size << "\texpression on both sides constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	std::cout << "Created " << 2*size*size << " constraints in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	tmr.reset();
	cplex.extract(model);
	std::cout << "Extracted model with " << cplex.getNcols() << " vars & "
			  << cplex.getNrows() << " rows into cplex\n" 
			  << "Loaded in " << tmr.elapsedSeconds() << " sec (CPU) "
			  << tmr.elapsedWallTime() << " sec wall " 
			  << timer.elapsedSeconds() << " total CPU" << std::endl;
	timer.reset();
	cplex.exportModel("CCC.lp");
 	//mip.writeLP("xxx.lp");
 	std::cout << "Wrote CCC.lp in "
 			  << timer.elapsedSeconds() << " sec (CPU) "
 			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	std::cout << "Total time: "
			  << totalTime.elapsedSeconds() << " sec (CPU) "
			  << totalTime.elapsedWallTime() << " sec wall-time" << std::endl;
} // end createLargeMIP_CPLEX

void createLargeMIP_GUROBI(GRBModel &model,int size)
// create problem with size*size variables and 2*size constraints
{
	qol::CpuTimer timer,totalTime;
	std::vector<GRBVar *> x(size);
	for(int i=0;i<size;++i){
		x[i] = model.addVars(size);
		model.update();			// is that needed here?
		for(IloInt j=0;j<size;++j){
			x[i][j].set(GRB_DoubleAttr_LB,0);
			x[i][j].set(GRB_DoubleAttr_UB,1);
			x[i][j].set(GRB_StringAttr_VarName,
						boost::str(boost::format("x_%03d_%03d")%i%j));
			x[i][j].set(GRB_CharAttr_VType,GRB_BINARY);
			x[i][j].set(GRB_DoubleAttr_Obj,(i*j/double(size)));
			
		}
	}
	model.set(GRB_IntAttr_ModelSense,1); // minimize
	model.update();
	std::cout << "Created " << size*size << " variables in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " wall-time" << std::endl;
	timer.reset();
	qol::CpuTimer tmr;
	for(int cnt=0;cnt < 10;++cnt){
		for(int i=0;i<size;++i){
			GRBLinExpr sum;
			for(int j=0;j<size;++j)
				sum += x[i][j];
			model.addConstr( sum == 1.0,
							 boost::str(boost::format("LHS_%d_%d") % i%cnt ));
		}
	}
	model.update();
	std::cout << size << "\tsimple sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt)	
		for(int j=0;j<size;++j){
			GRBLinExpr sum;
			for(int i=0;i<size;++i)
				sum += 1.0*x[i][j];
			model.addConstr( sum == 1.0,
							 boost::str(boost::format("RHS_%d_%d") % j%cnt ));
		}
	model.update();
	std::cout << size << "\tcoeff. sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		GRBLinExpr sumAll;
		for(int i=0;i<size;++i){
			GRBLinExpr sum1,sum2;
			for(int j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(int j=size/3;j<size;++j)
				sum2 += x[i][j];
			sumAll += sum1 + sum2;
		}
		model.addConstr( sumAll >= size,
						 boost::str(boost::format("all_%d") % cnt));
	}
	model.update();
	std::cout << 10 << "\texpr. sum constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	tmr.reset();
	for(int cnt=0;cnt < 10;++cnt){
		for(int i=0;i<size;++i){
			GRBLinExpr sum1,sum2;
			for(int j=0;j<2*size/3;++j)
				sum1 += x[i][j];
			for(int j=size/3;j<size;++j)
				sum2 += x[i][j];
			model.addConstr( sum1 + 10 >= sum2+5,
							 boost::str(boost::format("lhs_rhs_%d") % cnt));
		}
	}
	model.update();
	std::cout << size << "\texpression on both sides constraints in "
			  << tmr.elapsedSeconds() << " CPU sec\n";
	std::cout << "Created "
			  << model.get(GRB_IntAttr_NumConstrs) << " constraints in "
			  << timer.elapsedSeconds() << " sec (CPU) "
			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	timer.reset();
// 	std::cout << "Extracted model with "
// 			  << model.get(GRB_IntAttr_NumVars) << " vars & "
// 			  << model.get(GRB_IntAttr_NumConstrs) << " rows into cplex in " 
// 			  << timer.elapsedSeconds() << " sec (CPU) "
// 			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
 	model.write("GGG.lp");
 	std::cout << "Wrote GGG.lp in "
 			  << timer.elapsedSeconds() << " sec (CPU) "
 			  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
	std::cout << "Total time: "
			  << totalTime.elapsedSeconds() << " sec (CPU) "
			  << totalTime.elapsedWallTime() << " sec wall-time" << std::endl;
} // end createLargeMIP_GUROBI

int main(int argc, char* argv[])
{
    // QolColFormulation is not intended to solve a model. It is only used
    // to store variables and constraints
    qol::MIP *mip;
	if(argc <= 2){
		std::cout <<"USAGE: " << argv[0] << " size q|c|g|C|G\n"
				  << "Test creating large MILP using \n"
				  << "  q - QOL (data only) QolColFormulation\n"
				  << "  c - CPLEX  CplexFormulation\n"
				  << "  g - Gurobi GurobiFormulation\n"
				  << "  C - CPLEX  Concert/native C++ interface\n"
				  << "  G - Gurobi Native C++ interface\n"
			;
		return 1;
	}
	int size = atoi(argv[1]);
	if(size <=0 ) {
		std::cout << "Invalid size " << size << std::endl;
		return 2;
	}
	for(int i=2;i<argc;++i){
		switch(argv[i][0]){
			case 'c':
				std::cout << "Testing with CPLEX "
						  << "============================================\n";
				mip = new qol::CplexFormulation();
				break;
			case 'C':
				std::cout << "Testing with CPLEX Concert (native) "
						  << "============================================\n";
				mip = 0;
				{
					IloEnv envC;
					IloCplex cplex(envC);
					createLargeMIP_CPLEX(envC,cplex,size);
				}
				break;
			case 'g':
				std::cout << "Testing with Gurobi "
						  << "============================================\n";
				mip = new qol::GurobiFormulation();
				break;
			case 'G':
				std::cout << "Testing with Gurobi (native) "
						  << "============================================\n";
				mip = 0;
				{
					GRBEnv envG;
					GRBModel gurobi(envG);
					createLargeMIP_GUROBI(gurobi,size);
				}
				break;				
			case 'q':
				std::cout << "Testing with Qol Column Formulation "
						  << "============================================\n";
				mip = new qol::QolColFormulation();
				break;
			default:
				std::cout << "Unkown formulation type " << argv[i] << "\n";
				mip=0;
		}
		if(mip){
			createLargeMIP(*mip,size);
			qol::CpuTimer timer;
			if(argv[i][0] == 'q') {
				switch(argv[i][1]){
					case 'C': case 'c': {
						qol::CplexFormulation cpx(*(qol::QolColFormulation *)mip);
						std::cout << "Loaded into CPLEX from QOL in "
								  << timer.elapsedSeconds() << " sec (CPU) "
								  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
						timer.reset();
						cpx.writeLP("ccc.lp");
						std::cout << "Wrote ccc.lp in "
								  << timer.elapsedSeconds() << " sec (CPU) "
								  << timer.elapsedWallTime() << " sec wall-time\n";
					} break;
					case 'G': case 'g': {
						qol::GurobiFormulation grb(*(qol::QolColFormulation *)mip);
						std::cout << "Loaded into Gurobi from QOL in "
								  << timer.elapsedSeconds() << " sec (CPU) "
								  << timer.elapsedWallTime() << " sec wall-time" << std::endl;
						timer.reset();
						grb.writeLP("ggg.lp");
						std::cout << "Wrote ggg.lp in "
								  << timer.elapsedSeconds() << " sec (CPU) "
								  << timer.elapsedWallTime() << " sec wall-time\n";
					}break;
					default:
						mip->writeLP("xxx.lp");
						std::cout << "Wrote xxx.lp in "
								  << timer.elapsedSeconds() << " sec (CPU) "
								  << timer.elapsedWallTime() << " sec wall-time\n";
				}
			}
			delete mip;
		}
	}
    //Other Formulations provided (that provide a solve() method)
    // qol::CpxFormulation *m = new qol::CpxFormulation ();
    // qol::GurobiFormulation *m = new qol::GurobiFormulation();

	return 0;
}

