/***************************************************************
* Constraint Programming / Beam Search / Ant Colony Optimisation Hybrid
* Author: Andreas.Ernst@monash.edu  December 2015
****************************************************************/


#include "CPBeamACO.h"
#include "QolRandom.h"
#include "CpuTimer.h"
#include <boost/foreach.hpp>
#include <boost/thread.hpp>
using namespace qol;

//#define TRACE_CP(X) {std::cout << X << std::endl;}
//#define TRACE_CP(X) { if(state==OK && verify(1) ==Infeasible) { std::cout<< X << std::endl; }}
#define TRACE_CP(X)  {}// no debugging

CPstate::Consistency CPstate::updateRowBounds(Index v,BoundType bnd) // update contribution to max/min row values
{
	if(bnd&LOWER && collb[v] != 0){
		for (Index c = 0; c < mip->matrix[v].size(); ++c) {
			const double coeff = mip->matrix[v][c].first;
			const Index r = mip->matrix[v][c].second;
			if (coeff > 0){
				rowlb[r] += coeff*collb[v];
			}else{ 
				rowub[r] += coeff*collb[v];
			}	
		}
	}
	if(bnd&UPPER && colub[v] != 0) {
		for (Index c = 0; c < mip->matrix[v].size(); ++c) {
			const double coeff = mip->matrix[v][c].first;
			const Index r = mip->matrix[v][c].second;
			if (coeff > 0) {
				rowub[r] += coeff*colub[v];
			} else {
				rowlb[r] += coeff*colub[v];
			}
		}
	}
	// check all rows as a change might have occured due to reverseRowBounds()
	for (Index c = 0; c < mip->matrix[v].size(); ++c) {
		const Index r = mip->matrix[v][c].second;
		if( (mip->sense[r] != Constraint::LE && rowub[r] < mip->rhs[r] - eps)  ||
			(mip->sense[r] != Constraint::GE && rowlb[r] > mip->rhs[r] + eps) )
			return (state=Infeasible);
	}
	if( isFixed(v) )
		for (Index c = 0; c < mip->matrix[v].size(); ++c) 
			if(freeCnt[mip->matrix[v][c].second]) --freeCnt[mip->matrix[v][c].second];
	return OK;
}
void CPstate::reverseRowBounds(Index v,BoundType bnd) // reverse contribution
{
	if(bnd&LOWER && collb[v] != 0){
		for (Index c = 0; c < mip->matrix[v].size(); ++c) {
			const double coeff = mip->matrix[v][c].first;
			const Index r = mip->matrix[v][c].second;
			if (coeff > 0) {
				rowlb[r] -= coeff*collb[v];
			} else 
				rowub[r] -= coeff*collb[v];
		}
	}
	if(bnd&UPPER && colub[v] != 0) {
		for (Index c = 0; c < mip->matrix[v].size(); ++c) {
			const double coeff = mip->matrix[v][c].first;
			const Index r = mip->matrix[v][c].second;
			if (coeff > 0) {
				rowub[r] -= coeff*colub[v];
			} else 
				rowlb[r] -= coeff*colub[v];
		}
	}
}


bool CPstate::checkRowInactive(Index r) // recalculate maxChg[r],freeCnt[r] and return true if constraint can't be violated
{
	freeCnt[r] = 0; 
	maxChg[r] = 0;
	for (Index i = 0; i < mip->rowMatrix[r].size(); ++i) {
		const double coeff = mip->rowMatrix[r][i].first;
		const Variable v = mip->getVar(mip->rowMatrix[r][i].second);
		if (isFixed(v) || coeff==0) continue;
		++freeCnt[r];
		maxChg[r] = std::max(maxChg[r], fabs(coeff*(colub[v] - collb[v])));
	}
	if(freeCnt[r] <= 0) return true;
	switch(mip->getSense(r)){
	case Constraint::LE: 
		return (rowub[r] <= mip->getRHS(r)+eps);
	case Constraint::GE:
		return (rowlb[r] >= mip->getRHS(r)-eps);
	case Constraint::EQ:
		return ( mip->getRHS(r)-eps <= rowlb[r]  && rowub[r] <=  mip->getRHS(r)+eps );
	}
	return false; // never reached
}

void CPstate::init(const RowColFormulation *_mip) // initialise all data structures 
{ 
	mip = _mip;
	collb = mip->collb;
	colub = mip->colub;
	rowlb.resize(mip->nConstr(), 0);
	rowub.resize(mip->nConstr(), 0);
	maxChg.resize(mip->nConstr(),0); // max change in LHS for any row r due to a single variable moving
	freeCnt.resize(mip->nConstr(),0); // number of non-fixed variables for the constraint
	rows.reserve(mip->nConstr()); // rows that are still worth checking
	cutoff = qol::inf;
	state = OK;
	for (Index v = 0; v < mip->nVar(); ++v){
		if(fabs(mip->getVarLB(v)) > 1e10 || fabs(mip->getVarUB(v)) > 1e10){
			std::cout << "WARNING: numerical issues " << mip->getVarLB(v) 
				<< " <= " << mip->getVarName(v) << " <= " << mip->getVarUB(v) << std::endl;
		}
		updateRowBounds(v);
	}
	TRACE_CP("Initialised row bounds");
	checkRowBnds();
	state=OK; // make sure we don't have "updateRowBounds" setting state incorrectly
	for (Index v = 0; v < mip->nVar(); ++v)
		if( collb[v] > colub[v] + eps)
			state=Infeasible; // don't expect this to happen
	//changed.clear();
	for (Index r = 0; r < rowlb.size(); ++r){
		if ((mip->sense[r] != Constraint::GE && rowlb[r] > mip->rhs[r] + eps) || 
			(mip->sense[r] != Constraint::LE && rowub[r] < mip->rhs[r] - eps) ) {
			state = Infeasible;
			return;
		}
		if( ! checkRowInactive(r) ) { // re-calculates maxChg & freeCnt;
			rows.push_back(r); // this row needs to be checked still
		}
	}
	constrCnt.resize(mip->nVar(),0);
	BOOST_FOREACH(Index r,rows){
		for(size_t i=0;i<mip->rowMatrix[r].size();++i){
			Index v = mip->rowMatrix[r][i].second;
			if( ! isFixed(v) ) ++constrCnt[v];
		}
	}
	TRACE_CP("Doing initial propagate");
	//verify();
	propagate();
	TRACE_CP("CP state initialised");
}
void CPstate::setCutoff(double newUB)
{
	if (newUB < cutoff){
		cutoff = newUB;
		if (getObjectiveBound() > cutoff)
			state = Infeasible;
		else
			propagate();
	}
}
void CPstate::deactivateRow(Index r,size_t pos)// update counters - this row no longer active
	// optional input: position of r in rows[] or pos=InvalidIndex
{
	//TRACE_CP(r <<":" << mip->getConstrName(r) << " -- deactivating");
	for(size_t i=0;i<mip->rowMatrix[r].size();++i)
		if(constrCnt[mip->rowMatrix[r][i].second]) --constrCnt[mip->rowMatrix[r][i].second];
	while( rows[pos] != r && pos < rows.size()) ++pos;
	if( rows[pos]==r) {
		rows[pos] = rows.back();
		rows.pop_back();
	}
	freeCnt[r] = 0; // don't care any more
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
** WARNING: THIS CODE DOESN'T WORK CORRECTLY WHEN SOME OF THE VARIABLE 
   BOUNDS OR COEFFICIENTS ARE VERY DIFFERENT ORDERS OF MAGNITUDE. IN
   PARTICULAR IT FAILS WHEN SOME OF THE BOUNDS ARE "INFINITY" (eg 
   variables that are just constrained to positive with no upper bound)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
CPstate::Consistency CPstate::propagate()  // returns true if is feasible 
{   boost::this_thread::interruption_point();  // allows frequent checking for interruptions
	if (state != OK) return state; // no point if Infeasible or Fixed already
	//checkRowBnds();
	changed = false;//changed.clear();
	for (Index ridx = 0; ridx < rows.size(); ++ridx) {
		const Index r = rows[ridx];
		if(freeCnt[r] <= 0) { 
			//TRACE_CP(r <<":" << mip->getConstrName(r) << " all variables fixed");
			deactivateRow(r,ridx--); continue; 
		} // don't look at constraint again
		//TRACE_CP(r <<":" << mip->getConstrName(r) << " checking " << freeCnt[r] <<" free, range: "
		//	<< rowlb[r] << " to " << rowub[r] << " " << mip->getSense(r) << " " << mip->getRHS(r));
		if(freeCnt[r] > 1) { // otherwise want to update bound for remaining var & remove constraint
			if(mip->getSense(r) == Constraint::LE){
				if(rowlb[r] > mip->getRHS(r)+eps)
					return (state = Infeasible);
				if(rowub[r] <= mip->getRHS(r)+eps){ // cannot violate this constraint ever
					//TRACE_CP(r <<":" << mip->getConstrName(r) << " <= inactive");
					deactivateRow(r,ridx--); continue; 
				}
				if(rowlb[r]+maxChg[r] <= mip->rhs[r]) 
					continue; // no effect from propagation possible
			}else if(mip->getSense(r) == Constraint::GE) {
				if(rowub[r] < mip->getRHS(r)-eps)
					return (state = Infeasible);
				if(rowlb[r] >= mip->getRHS(r)-eps){ // cannot violate this constraint ever
					//TRACE_CP(r <<":" << mip->getConstrName(r) << " >= inactive");
					deactivateRow(r,ridx--); continue; 
				}
				if(rowub[r] - maxChg[r] >= mip->rhs[r])
					continue; // cannot propagate
			}else if(mip->getSense(r) == Constraint::EQ) {
				if( rowlb[r] > mip->getRHS(r)+eps || rowub[r] < mip->getRHS(r)-eps)
					return (state = Infeasible);
				if( rowub[r] - maxChg[r] >= mip->rhs[r] && rowlb[r]+maxChg[r] <= mip->getRHS(r))
					continue;// no effect from propagation possible
			}
		}
		// can't short-cut so check each coefficient & recompute maxChg
		// don't really need to recompute freeCnt but no harm done
		freeCnt[r] = 0; maxChg[r] = 0;
		for (Index i = 0; i < mip->rowMatrix[r].size() && state == OK; ++i) {
			const double coeff = mip->rowMatrix[r][i].first;
			const Variable v = mip->getVar(mip->rowMatrix[r][i].second);
			if (isFixed(v) || coeff==0) continue;
			++freeCnt[r];
			maxChg[r] = std::max(maxChg[r], fabs(coeff*(colub[v] - collb[v])));
			if (coeff > 0) {
				if (mip->sense[r] != Constraint::LE &&
					rowub[r] - coeff*(colub[v] - collb[v]) < mip->rhs[r] - eps) {
					setVarLB(v, (mip->rhs[r] - rowub[r] + coeff*colub[v]) / coeff );
					TRACE_CP(r <<":" << mip->getConstrName(r) << " tighter LB for " << 
							 v <<":" << mip->getVarName(v) << " = " << collb[v]);
				}
				if (mip->sense[r] != Constraint::GE &&
					rowlb[r] + coeff*(colub[v] - collb[v]) > mip->rhs[r] + eps) {
					setVarUB(v, (mip->rhs[r] - rowlb[r] + coeff*collb[v]) / coeff);
					TRACE_CP(r <<":" << mip->getConstrName(r) << " tighter UB for " << 
							 v <<":" << mip->getVarName(v) << " = " << colub[v]);
				}
			} else if(coeff < 0){ // coeff < 0
				if (mip->sense[r] != Constraint::LE &&
					rowub[r] + coeff*(colub[v] - collb[v]) < mip->rhs[r] - eps) {
					setVarUB(v, (mip->rhs[r] - rowub[r] + coeff*collb[v]) / coeff);
					TRACE_CP(r <<":" << mip->getConstrName(r) << " tighter UB for " << 
							 v <<":" << mip->getVarName(v) << " = " << colub[v]);
				}
				if (mip->sense[r] != Constraint::GE &&
					rowlb[r] - coeff*(colub[v] - collb[v]) > mip->rhs[r] + eps) {
					setVarLB(v, (mip->rhs[r] - rowlb[r] + coeff*colub[v]) / coeff);
					TRACE_CP(r <<":" << mip->getConstrName(r) << " tighter LB for " << 
							 v <<":" << mip->getVarName(v) << " = " << collb[v]);
				}
			}
		} // end loop over non-zeros for row
		if(freeCnt[r] <= 1){
			//TRACE_CP(r <<":" << mip->getConstrName(r) << " all bar " <<freeCnt[r] <<" variables fixed, range:"
			//	<< rowlb[r] << " to " << rowub[r] << " " << mip->getSense(r) << " " << mip->getRHS(r));
			deactivateRow(r,ridx--);
		}else {
			//TRACE_CP(r <<":" << mip->getConstrName(r) << " still " << freeCnt[r] <<" free, range: "
			//<< rowlb[r] << " to " << rowub[r] << " " << mip->getSense(r) << " " << mip->getRHS(r));
		}
	} // end loop over constraints
	//verify();
	if(state != OK) return state;
	double lb = getObjectiveBound();
	if (lb > cutoff + eps) return (state = Infeasible);
	int remainingVars=0;
	for (Index v = 0; v < mip->nVar(); ++v) {
		if (isFixed(v)) continue;
		if(constrCnt[v] <= 0) { // all constraints have only this variable free
			if (mip->getObjCoeff(v) >= 0) {
				setVarUB(v,collb[v]); // fix at lower bound
				TRACE_CP(v <<":" << mip->getVarName(v) << " only bound constrained fixing lb " << collb[v]);
			} else {
				setVarLB(v,colub[v]); // fix at upper bound
				TRACE_CP(v <<":" << mip->getVarName(v) << " only bound constrained fixing to ub " << colub[v]);
			}
		} else {
			const double c = mip->getObjCoeff(v);
			//TRACE_CP(v <<":" << mip->getVarName(v) << " cost=" << c << " cutoff=" << cutoff 
			//	<< " constrCnt=" << constrCnt[v]);
			if (c > 0 && lb + c*(colub[v] - collb[v]) > cutoff + eps) {
				setVarUB(v, (cutoff - (lb - c*collb[v])) / c);
				TRACE_CP(v <<":" << mip->getVarName(v) << " tigtening for cutoff " << cutoff << 
					     " UB = " << colub[v]);
				if (collb[v] > colub[v] + eps) return (state = Infeasible);
			} else if (c < 0 && lb - c*(colub[v] - collb[v]) > cutoff + eps) {
				setVarLB(v, (cutoff - lb + c*colub[v]) / c);
				TRACE_CP(v <<":" << mip->getVarName(v) << " tigtening for cutoff " << cutoff << 
					     " LB = " << collb[v]);
				if (collb[v] > colub[v] + eps) return (state = Infeasible);
			}
		}
		if(state != OK) return state;
		if(!isFixed(v)) ++remainingVars;
	}
	//verify();
	if (remainingVars == 0) return (state = Fixed);
	if ( ! changed ) return state;
	return propagate(); // tail recursion - check for further changes
} // end CPstate::propagate()

CPstate::Consistency CPstate::verify(int verbose) const
	// chedck current state is OK (verbose defaults to zero meaning don't print)
{ 
	if( state==Infeasible) return state;
	for (Index r = 0; r < rowlb.size(); ++r) {
		double minlhs = 0, maxlhs=0;
		for (Index i = 0; i < mip->rowMatrix[r].size(); ++i) {
			const double coeff = mip->rowMatrix[r][i].first;
			const Variable v = mip->getVar(mip->rowMatrix[r][i].second);
			if( coeff > 0){
				minlhs += coeff*collb[v];
				maxlhs += coeff*colub[v];
			}else{
				minlhs += coeff*colub[v];
				maxlhs += coeff*collb[v];
			}
		}
		if (maxlhs < mip->rhs[r] - eps && mip->sense[r] != Constraint::LE){
			if(verbose > 0)
				std::cerr << "ERROR: violating " << mip->getConstrName(r) << " by " << mip->rhs[r]-maxlhs 
						<< " min/maxlhs=" << minlhs << " / " << maxlhs 
						<< " row lb/ub = " << rowlb[r] << " / " << rowub[r] << std::endl;
			return Infeasible;
		}
		if (minlhs > mip->rhs[r] + eps && mip->sense[r] != Constraint::GE){
			if(verbose > 0)
				std::cerr << "ERROR: violating " << mip->getConstrName(r) << " by " << minlhs-mip->rhs[r] 
						<< " min/maxlhs=" << minlhs << " / " << maxlhs 
						<< " row lb/ub = " << rowlb[r] << " / " << rowub[r] << std::endl;
			return Infeasible;
		}
	}
	bool fixed=true;
	for(Index v=0;v<collb.size();++v) if( ! isFixed(v) ) fixed = false;
	return (fixed ? Fixed : OK);
}

CPstate::Consistency CPstate::checkSolution(const std::vector<double> &solution,
	std::vector<double> *LHS) const
{
	Consistency s = OK;
	if (LHS != 0) LHS->resize(rowlb.size());
	for (Index r = 0; r < rowlb.size(); ++r) {
		double lhs = 0;
		for (Index i = 0; i < mip->rowMatrix[r].size(); ++i) {
			const double coeff = mip->rowMatrix[r][i].first;
			const Variable v = mip->getVar(mip->rowMatrix[r][i].second);
			lhs += coeff*solution[v];
		}
		if (LHS) (*LHS)[r] = lhs;
		if (s == OK) { // check if we are still OK
			if (lhs < mip->rhs[r] - eps && mip->sense[r] != Constraint::LE){
				std::cerr << "ERROR: solution infeasible " << mip->getConstrName(r) << " by " << mip->rhs[r]-lhs << std::endl;
				s = Infeasible;
			}
			if (lhs > mip->rhs[r] + eps && mip->sense[r] != Constraint::GE){
				std::cerr << "ERROR: solution infeasible " << mip->getConstrName(r) << " by " << lhs-mip->rhs[r] << std::endl;
				s = Infeasible;
			}
		}
	}
	return s;
}

/// getObjectiveBound() returns a (weak) lower/upper bound for
/// a minimisation problem
double CPstate::getObjectiveBound() const
{
	double obj = 0;
	for (Index v = 0; v < mip->nVar(); ++v) {
		const double c = mip->getObjCoeff(v);
		if (c > 0)
			obj += c*collb[v];
		else if (c < 0)
			obj += c*colub[v];
	}
	return obj;
}

//************************************************************************************
//** CPBeamACO - Contraint Programming / Beam Search / Ant Colony Optimisation Hybrid
//************************************************************************************

void CPBeamACO::initVarGroups()
{
	std::vector<bool> included(model.nVar(), false);
	// find SOSequalOne (ie \sum x = 1) groups
	for (Index r = 0; r < model.sense.size(); ++r) {
		if (model.sense[r] != Constraint::EQ) continue;
		bool ok = true; Index nonFixed=0;
		for (Index i = 0; i < model.rowMatrix[r].size() && ok; ++i) {
			const qol::CoeffIdxPair &cv = model.rowMatrix[r][i];
			if (model.getVarType(cv.second) != Variable::BINARY)
				ok = false;
			else if (fabs(cv.first - model.rhs[r]) > eps) 
				ok = false;
			else if (included[cv.second])
				ok = false; // already included
			else if(state && !state->isFixed(cv.second))
				++nonFixed;
		}
		//std::cout << model.getConstrName(r) << " is " << ok << " with " << nonFixed << " nonFixed\n";
		if (ok && nonFixed > 1) {
			SOSequalOne *sos = new SOSequalOne();
			sos->setIndex(varGroup.size());
			varGroup.push_back(sos);
			for (Index i = 0; i < model.rowMatrix[r].size(); ++i){
				if( state && !state->isFixed(model.rowMatrix[r][i].second) )
					sos->addVar(model.rowMatrix[r][i].second);
				included[model.rowMatrix[r][i].second] = true;
			}
		}
	}
	if(state) 
		for(Index v=0;v<model.nVar();++v) if(state->isFixed(v)) included[v] = true;
	size_t SOS0cnt = varGroup.size();
	// find SOSatmostOne (\sum x <= 1) groups
	for (Index r = 0; r < model.sense.size(); ++r) {
		if (! ((model.sense[r] == Constraint::LE && model.rhs[r] > 0) ||
			(model.sense[r] == Constraint::GE && model.rhs[r] < 0)) )
			continue;
		bool ok = true; int varCnt =0;
		for (Index i = 0; i < model.rowMatrix[r].size() && ok; ++i) {
			const qol::CoeffIdxPair &cv = model.rowMatrix[r][i];
			if (model.getVarType(cv.second) != Variable::BINARY)
				ok = false;
			else if (fabs(cv.first - model.rhs[r]) > eps)
				ok = false;
			else if ( ! included[cv.second])// OK just ignore included
				++varCnt;
		}
		if (ok && varCnt > 1) {
			SOSatmostOne *sos = new SOSatmostOne();
			sos->setIndex(varGroup.size());
			varGroup.push_back(sos);
			for (Index i = 0; i < model.rowMatrix[r].size(); ++i)
				if (!included[model.rowMatrix[r][i].second]) {
					sos->addVar(model.rowMatrix[r][i].second);
					included[model.rowMatrix[r][i].second] = true;
				}
		}
	}
	size_t SOS1cnt=varGroup.size()-SOS0cnt;
	// add remaining binary variables
	for (Index v=0;v<model.nVar();++v){
		if(included[v]) continue; 
		// check if this variable will get fixed by others:
		// (a) if all constraints involve variables already fixed (included)
		// (b) As a special case we check if we have v >= sum(binaries in B_i) for some B_i
		//     and v <= sum( sum(binaries in B_i) for all i) with v binary 
		std::set<Index> atLeast,atMost; 
		size_t fixCnt=0;
		for (Index c = 0; c < model.matrix[v].size(); ++c) { // for each constraint
			double coeff=model.matrix[v][c].first;
			const Index r = model.matrix[v][c].second;		// row
			if(coeff == 0){ ++fixCnt; continue; }// shouldn't really happen
			bool allFixed=true;
			bool checkAtLeast=false,checkAtMost = false;
			if( model.getVarType(v) == Variable::BINARY && model.getRHS(r) == 0) {
				if( (coeff > 0 && model.getSense(r) == Constraint::GE) || // v - sum(binaries) >= 0
					(coeff < 0 && model.getSense(r) == Constraint::LE)){  // sum(binaries) - v <= 0
					checkAtLeast=true; 
				}else if((coeff > 0 && model.getSense(r) == Constraint::LE) || // v - sum(binaries) <= 0
						 (coeff < 0 && model.getSense(r) == Constraint::GE)){  // sum(binaries) - v >=0
					checkAtMost=true;  
				}
			}
			for(Index i=0;i<model.rowMatrix[r].size() && allFixed;++i){
				const Index vv = model.rowMatrix[r][i].second;
				if(vv == v) continue;
				if(!included[vv]){ allFixed=false; break;}
				if(checkAtLeast){
					if(fabs(coeff+ model.rowMatrix[r][i].first) < eps)
						atLeast.insert(vv);
					else
						atLeast.insert(InvalidIndex); // won't appear in atMost
				}else if(checkAtMost) {
					if(fabs(coeff+ model.rowMatrix[r][i].first) < eps)
						atMost.insert(vv);
					else
						atMost.insert(v); // won't appear in atLeast
				}
			}
			if(allFixed){
				if(model.getSense(r) == Constraint::EQ){ // this forces the constraint to be fixed
					fixCnt = model.matrix[v].size();
					break;
				}
				++fixCnt;
			}
		}
		if( fixCnt == model.matrix[v].size()){ // all constraints are just variable bounds
			//TRACE_CP( model.getVarName(v) << " just bounds left " );
			included[v] = true;
		}else if( ! atLeast.empty() && ! atMost.empty() && atLeast == atMost){
			included[v] = true; // binary fixed variable = max of other binaries
			//TRACE_CP(model.getVarName(v) << " dependent binary ");
		}else{
			if( model.getVarType(v) == Variable::CONTINUOUS ) { // shouldn't really happen
				const double lb=state->getVarLB(v), ub=state->getVarUB(v);
				//TRACE_CP(model.getVarName(v) << " single var in range " << lb << " to " << ub );
				varGroup.push_back( new SingleVarGroup(v,state->getVarLB(v),state->getVarUB(v),10) );
			}else{
				const int lb=state->getVarLB(v), ub=state->getVarUB(v);
				//TRACE_CP(model.getVarName(v) << " single var in range " << lb << " to " << ub);
				varGroup.push_back(new SingleVarGroup(v, double(lb), ub, ub-lb) );
			}
		}
	} // end loop over variables
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "Variable Groups " << varGroup.size() << " = " << SOS0cnt << " + " << SOS1cnt << " + " 
				<< (varGroup.size() - SOS0cnt - SOS1cnt) << " with " << included.size() << " variables"<< std::endl;

} // end initVarGroups()

void CPBeamACO::init() 
{
	if (model.rowMatrix.size() != model.rhs.size())
		model.load(); // ensure rowMatrix is set
	if (state) state->propagate(); else state = new CPstate(model);
	initVarGroups(); 
}

// main method for implementing solver
qol::Status CPBeamACO::solveBeam()
{
	bestObj = qol::inf;
	// run for given wall time limit
	const int threads = 1;// params.getParamValue(THREADS);
	const int beamWidth=4;
	const int beamMult=4;
	const int maxBeam=beamWidth*beamMult;
	Status status = qol::FAILED;
	const double tilim = params.getParamValue(TIMELIMIT);
	CpuTimer timer(std::max((2+threads)*tilim,1.0),std::max(tilim,1.0));
	int trial=-1; // not yet started
	try{
	//const int threads = std::max(omp_get_max_threads() - 2, 1);
	//const int threads = omp_get_max_threads();
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "CPBeamACO " << model.nConstr() << " rows and " << model.nVar() << " columns. Running for " << timer.timeRemaining() << " seconds" << std::endl;
	status = solveRelaxed();
	if (state->getState() == CPstate::Infeasible) return status;
	Random rand;
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "Lower bound " << state->getObjectiveBound() << std::endl;
	if(state->getState()==CPstate::Fixed) {
		updateBest(state);
		if (params.getParamValue(VERBOSITY) > 0){
			std::cout << "Found optimal at root node - all variables fixed."  << std::endl;
			std::cout << "\t Completed " << 0 << " iterations in " << timer.elapsedSeconds()
				  << " sec, found best objective " << bestObj << std::endl;
		}
		return qol::OPTIMAL;
	}
	//for(size_t vg=0;vg<varGroup.size();++vg) varGroup[vg].print(*state);
	const int itlimit = 1000000; //params.getParamValue(ITLIMIT)
	double lastPrintSec = timer.elapsedWallTime();
	//const int nThreads= std::max(omp_get_max_threads()-2, 1);
	//const int nThreads = omp_get_max_threads();
	const int nThreads = 1;// omp_get_max_threads();
	if(params.getParamValue(VERBOSITY) > 0)
		std::cout << "Running with " << nThreads << " threads\n";
	std::vector<std::vector<Index> > choices(beamWidth);
	for (trial = 0; trial < itlimit && !timer.timeLimitReached(); ++trial) {
		std::vector<CPstate *>ants;
		ants.reserve(beamWidth);
		ants.push_back(new CPstate(*state));
		double iterBest=qol::inf; int completed=0;
		for(size_t vg=0;vg<varGroup.size() && ! ants.empty();++vg){
			// children have a value and a pointer to a state
			std::vector<std::pair<double,CPstate *> > child(maxBeam,std::make_pair(qol::inf,(CPstate*)0));
			int nChild=0,nFeas=0;
			for(size_t parent=0;parent<ants.size();++parent){
				const Index nChoices=varGroup[vg].getFeasible(*ants[parent],choices[parent]);
				if(nChoices == 0){ // variable group already fixed - just pass parent to next generation
					child[parent*(size_t)beamMult].second = ants[parent];
					child[parent*(size_t)beamMult].first = ants[parent]->getObjectiveBound();
					ants[parent] = 0;
				}else{
					for(Index c=0;c<std::min(nChoices,(Index)beamMult);++c){
						CPstate *tmp=new CPstate(*ants[parent]);
						child[parent*(size_t)beamMult+c].second = tmp;
						int choice = rand.uniform_int(c,nChoices); // ignores previous choices
						if(choice != c) std::swap(choices[parent][c],choices[parent][choice]); // put to front
						//std::cout << "\t\t" << "parent " << parent << " choosing " <<choices[parent][c] << " of " << nChoices << std::endl;
						varGroup[vg].setVariables(*tmp,choices[parent][c]);
						tmp->propagate(); ++nChild;
						if(tmp->getState() != CPstate::Infeasible){
							child[parent*(size_t)beamMult+c].first = tmp->getObjectiveBound();
							++nFeas;
						}
					}
					delete ants[parent]; // replaced by children
				}
			}
			// now select best children
			std::sort(child.begin(),child.end());
			//std::cout << "\t" << trial << ":" << vg << ": " 
			//	<< ants.size() << " parents, " << nChild << " children, " << nFeas << " feas." << std::endl;
			ants.clear();
			for(size_t c=0;c<child.size();++c){
				if(child[c].first < qol::inf && child[c].second->getState() == CPstate::Fixed) {
					updateBest(child[c].second); // hurray, found a solution (don't continue with this)
					iterBest = std::min(iterBest,child[c].first); 
					++completed;
					delete child[c].second;
					child[c].second=0;
				}else if(child[c].first < qol::inf && ants.size() < beamWidth){
					ants.push_back(child[c].second);
					child[c].second = 0;
				}else if(child[c].second != 0)
					delete child[c].second;
			}

		}		
		if(params.getParamValue(VERBOSITY) > 0 && timer.elapsedWallTime() > lastPrintSec+10.0){
			std::cout << "\t Iteration " << trial << " objective = " << iterBest 
				<< " tested " << completed << " best " << bestObj << std::flush << std::endl;
			lastPrintSec = timer.elapsedWallTime();
		}
		if (bestObj <= state->getObjectiveBound() + eps) { // found optimal solution
			if (params.getParamValue(VERBOSITY) > 0)
				std::cout << "\t Optimal solution " << bestObj << " at iteration " << trial << std::endl;
			status = qol::OPTIMAL; break;
		}else if (bestObj < state->getCutoff()) {
			if (params.getParamValue(VERBOSITY) > 0)
				std::cout << "\t New best solution " << bestObj << " at iteration " << trial << std::flush << std::endl;
			state->setCutoff(bestObj);
			if (state->propagate() != CPstate::OK) {
				status = qol::OPTIMAL; // all other solutions cut
				break;
			}
		}
	} // end loop over trials
	}catch(boost::thread_interrupted){
		status = qol::ABORTED;
	}catch(...){
		status = qol::FAILED;
	}
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "\t Completed " << trial << " iterations in " << timer.elapsedSeconds() << "/" << timer.elapsedWallTime()
				  << " sec, found best objective " << bestObj << std::endl;
	return status;
} // end solveBeam()


// main method for implementing solver
qol::Status CPBeamACO::solveExact()
{
	return solveBeam(); 
	bestObj = qol::inf;
	// run for given wall time limit
	const int threads = 1;// params.getParamValue(THREADS);
	//const int threads = std::max(omp_get_max_threads() - 2, 1);
	//const int threads = omp_get_max_threads();
	const double tilim = params.getParamValue(TIMELIMIT);
	CpuTimer timer(std::max((2+threads)*tilim,1.0),std::max(tilim,1.0));
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "CPBeamACO " << model.nConstr() << " rows and " << model.nVar() << " columns. Running for " << timer.timeRemaining() << " seconds" << std::endl;
	Status status = solveRelaxed();
	if (state->getState() == CPstate::Infeasible) return status;
	Random rand;
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "Lower bound " << state->getObjectiveBound() << std::endl;
	if(state->getState()==CPstate::Fixed) {
		updateBest(state);
		if (params.getParamValue(VERBOSITY) > 0){
			std::cout << "Found optimal at root node - all variables fixed."  << std::endl;
			std::cout << "\t Completed " << 0 << " iterations in " << timer.elapsedSeconds()
				  << " sec, found best objective " << bestObj << std::endl;
		}
		return qol::OPTIMAL;
	}
	//for(size_t vg=0;vg<varGroup.size();++vg) varGroup[vg].print(*state);
	const int itlimit = 1000000; //params.getParamValue(ITLIMIT)
	double lastPrintSec = timer.elapsedWallTime();
	int trial;
	//const int nThreads= std::max(omp_get_max_threads()-2, 1);
	//const int nThreads = omp_get_max_threads();
	const int nThreads = 1;// omp_get_max_threads();
	if(params.getParamValue(VERBOSITY) > 0)
		std::cout << "Running with " << nThreads << " threads\n";
	for (trial = 0; trial < itlimit && !timer.timeLimitReached(); ++trial) {
#       if 0 // parallel version
		std::vector<CPstate *>ants(nThreads);
		for(int a=0;a<nThreads;++a) ants[a] = new CPstate(*state);
#       pragma omp parallel for
		for(int a=0;a<nThreads;++a) {		
			for(size_t vg=0;vg<varGroup.size();++vg){
				const Index nChoices=varGroup[vg].nFeasible(*ants[a]);
				if(nChoices==0) continue; // all fixed, no choice avail
				const Index choice=(trial==0 && a== 0 ? 0:rand.uniform_int(0,nChoices));
				varGroup[vg].setFeasible(*ants[a],choice);
				if (ants[a]->propagate() == CPstate::Infeasible){
					//if (params.getParamValue(VERBOSITY) > 0 && timer.elapsedWallTime() > lastPrintSec+5){
					//	std::cout << "\t Iteration " << trial << " fixed " << vg << " of " 
					//			<< varGroup.size() << std::flush << std::endl; // variables branched/fixed by CP/undecided
					//	lastPrintSec = timer.elapsedWallTime();
					//}
					break; // no point continuing
				}else if(ants[a]->getState() == CPstate::Fixed)
					break;
			}
		}
		if(params.getParamValue(VERBOSITY) > 0 && timer.elapsedWallTime() > lastPrintSec+10.0){
			std::cout << "\t Iteration " << trial << " objective = " << ants[0]->getObjectiveBound() 
				<< (ants[0]->getState() == CPstate::Fixed ? " heuristic" : (ants[0]->getState() == CPstate::Infeasible ? " infeasible" : " OK"))
				<< " best " << bestObj << std::flush << std::endl;
			lastPrintSec = timer.elapsedWallTime();
		}

		for(int a=0;a<nThreads;++a) {	
			if( ants[a]->getState() == CPstate::Fixed && ants[a]->verify() != CPstate::Infeasible) updateBest(ants[a]);
			delete ants[a];
		}
#       else // serial version
		CPstate test(*state);
		if(trial%3==2){ // see if we can find a much better solution
			double lb = test.getObjectiveBound();
			double phantomBound = lb + rand.uniform(0.1,1.0)*(bestObj-lb);
			test.setCutoff(phantomBound);
			if (test.propagate() == CPstate::Infeasible) continue; // try again
		}
		for(size_t vg=0;vg<varGroup.size();++vg){
			const Index nChoices=varGroup[vg].nFeasible(test);
			//const Index choice=rand.uniform_int(0,nChoices);
			//if (trial == 0) choice = 0;
			const Index choice= (trial == 0) ? 0 : rand.uniform_int(0,nChoices);
			varGroup[vg].setFeasible(test,choice);
			//std::cout << "\t " << trial << ": var " << vg << " set to " << choice << "/" << nChoices << std::endl;
			if (test.propagate() == CPstate::Infeasible){
				if (params.getParamValue(VERBOSITY) > 0 && timer.elapsedWallTime() > lastPrintSec+5){
					std::cout << "\t Iteration " << trial << " fixed " << vg << " of " 
							<< varGroup.size() << std::flush << std::endl; // variables branched/fixed by CP/undecided
					lastPrintSec = timer.elapsedWallTime();
				}
				break; // no point continuing
			}else if(test.verify(params.getParamValue(VERBOSITY)) == CPstate::Infeasible){
				//if (params.getParamValue(VERBOSITY) > 0){
					std::cout << "Propagate failed " << trial << ": var " << vg << " set to " << choice << "/" << nChoices 
						<< " state=" << test.getState() << std::endl;
				//}
				break;
			}else if (test.getState() == CPstate::Fixed) {
				updateBest(&test); // hurray, found a solution
				break;
			}
		}
		if(params.getParamValue(VERBOSITY) > 0 && timer.elapsedWallTime() > lastPrintSec+10.0){
			std::cout << "\t Iteration " << trial << " objective = " << test.getObjectiveBound() 
				<< (test.getState() == CPstate::Fixed ? " heuristic" : (test.getState() == CPstate::Infeasible ? " infeasible" : " OK"))
				<< " best " << bestObj << std::flush << std::endl;
			lastPrintSec = timer.elapsedWallTime();
		}
		/*
		std::vector<Index> vars;
		vars.reserve(nVar());
		for (Index v = 0; v < nVar(); ++v){
			if (model.coltype[v] != qol::Variable::CONTINUOUS)
				vars.push_back(v);
		}
		std::random_shuffle(vars.begin(), vars.end());
		const size_t nInt = vars.size();
		for (Index v = 0; v < nVar(); ++v){
			if (model.coltype[v] == qol::Variable::CONTINUOUS)
				vars.push_back(v);
		}
		if (vars.size() > nInt)
			std::random_shuffle(vars.begin() + nInt, vars.end());
		int cnt=0,fixCnt=0;
		BOOST_FOREACH(Index v,vars) {
			++cnt;
			if (test.isFixed(v)) continue;
			++fixCnt;
			if (model.coltype[v] == qol::Variable::CONTINUOUS)
				test.fixVar(v, rand.uniform(test.getVarLB(v), test.getVarUB(v)));
			else
				test.fixVar(v, rand.uniform_int(test.getVarLB(v), test.getVarUB(v)));
			if (test.propagate() == CPstate::Infeasible){
				if (params.getParamValue(VERBOSITY) > 0)
					std::cout << "\t Iteration " << trial << " infeasible " << fixCnt << " + " 
							<< (cnt-fixCnt) << " + " << (vars.size() -fixCnt) 
							<< " vars branched/fixed/remaining" << std::endl; // variables branched/fixed by CP/undecided
				break; // no point continuing
			}
			if (test.getState() == CPstate::Fixed) {
				updateBest(&test); // hurray, found a solution
				break;
			}
		}
		*/
		if(test.getState() == CPstate::OK && test.verify() != CPstate::Infeasible){
			for(Index v=0;v<model.nVar();++v)
				if( ! test.isFixed(v) )
					std::cout << "Failed to fix " << model.getVarName(v) << " in range " << test.getVarLB(v) << "," << test.getVarUB(v) << std::endl;
		}
#       endif
		if (bestObj <= state->getObjectiveBound() + eps) { // found optimal solution
			if (params.getParamValue(VERBOSITY) > 0)
				std::cout << "\t Optimal solution " << bestObj << " at iteration " << trial << std::endl;
			status = qol::OPTIMAL; break;
		}else if (bestObj < state->getCutoff()) {
			if (params.getParamValue(VERBOSITY) > 0)
				std::cout << "\t New best solution " << bestObj << " at iteration " << trial << std::flush << std::endl;
			state->setCutoff(bestObj);
			if (state->propagate() != CPstate::OK) {
				status = qol::OPTIMAL; // all other solutions cut
				break;
			}
		}
		
	}
	if (params.getParamValue(VERBOSITY) > 0)
		std::cout << "\t Completed " << trial << " iterations in " << timer.elapsedSeconds() << "/" << timer.elapsedWallTime()
				  << " sec, found best objective " << bestObj << std::endl;
	return status;
}


/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
