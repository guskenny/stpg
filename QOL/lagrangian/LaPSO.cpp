/// Generic solver for combinatorial optimisation problems using a combination
/// of Wedelin's lagrangian relaxation heuristic and particle swarm
/// optimisation. This header file follows the general pattern of the Volume
/// Algorithm implementation in the COIN-OR library to make it easy to swap
/// algorithms and make comparisons between the two methods.
///
/// The current implementation only allows for binary variables.
///
/// Autor: Andreas Ernst   $ Date: January 2010 $  $ Revision: 0.0 $
#include "LaPSO.hpp"
#include "Random.h"
#include <algorithm>
#include <stdlib.h>
#include "anyoption.h"
#include <cstdlib>
#include <time.h>
#ifdef _OPENMP
#  include <omp.h>
#else
#  include <time.h>
#  define omp_get_wtime() (double)time(NULL)
#  define omp_set_num_threads(n) // do nothing
#endif

namespace LaPSO {
    lagged_fibonacci1279 _randGen;
    void Param::parse(int argc,const char **argv)
    {
		AnyOption parser(100);
		parser.setFlag("help");
		parser.setFlag("randomSeed");
		parser.setOption("subgradFactor"); 
		parser.setOption("subgradFmult"); 
		parser.setOption("subgradFmin"); 
		parser.setOption("perturbFactor");
		parser.setOption("velocityFactor");
		parser.setOption("globalFactor");
		parser.setOption("maxIter");
		parser.setOption("maxCPU");
		parser.setOption("maxWallTime");
		parser.setOption("nParticles");
		parser.setOption("printLevel");
		parser.setOption("printFreq");
		parser.setOption("heurFreq");
		parser.setOption("eps");
		parser.setOption("absGap");
		parser.setOption("relGap");
		parser.setOption("nCPU");
		if(argc == -1)		// abuse of function
			parser.processFile(*argv);
		else
			parser.processCommandArgs(argc,argv);
		if(parser.getFlag("help"))
			parser.printAutoUsage();
		if(parser.getFlag("randomSeed"))
			randomSeed = 0;	// 
		if(parser.getValue("subgradFactor"))
			subgradFactor = atof(parser.getValue("subgradFactor"));
		if(parser.getValue("subgradFmult"))
			subgradFmult = atof(parser.getValue("subgradFmult"));
		if(parser.getValue("subgradFmin"))
			subgradFmin = atof(parser.getValue("subgradFmin"));
		if(parser.getValue("perturbFactor"))
			perturbFactor = atof(parser.getValue("perturbFactor"));
		if(parser.getValue("velocityFactor"))
			velocityFactor = atof(parser.getValue("velocityFactor"));
		if(parser.getValue("globalFactor"))
			globalFactor = atof(parser.getValue("globalFactor"));
		if(parser.getValue("maxIter")) 	    
			maxIter = atoi(parser.getValue("maxIter"));
		if(parser.getValue("maxCPU")) 	    
			maxCPU = atof(parser.getValue("maxCPU"));
		if(parser.getValue("maxWallTime")) 	    
			maxWallTime = atof(parser.getValue("maxWallTime"));
		if(parser.getValue("nParticles")) 	    
			nParticles = atoi(parser.getValue("nParticles"));
		if(parser.getValue("printLevel")) 	    
			printLevel = atoi(parser.getValue("printLevel"));
		if(parser.getValue("printFreq")) 	    
			printFreq = atoi(parser.getValue("printFreq"));
		if(parser.getValue("heurFreq")) 	    
			heurFreq = atoi(parser.getValue("heurFreq"));
		if(parser.getValue("eps"))
			eps = atof(parser.getValue("eps"));
		if(parser.getValue("absGap"))
			absGap = atof(parser.getValue("absGap"));
		if(parser.getValue("relGap"))
			relGap = atof(parser.getValue("relGap"));
		if(parser.getValue("nCPU"))
			nCPU = std::max(1,atoi(parser.getValue("nCPU")));
    }
    void Param::parse(const char *filename){ parse(-1,&filename); }
    void Particle::resize(size_t numVar,size_t numConstr){
		dual.resize(numConstr,0.0);
		perturb.resize(numVar,0.0);
		x.resize(numVar,0);
		viol.resize(numConstr,0.0);
		rc.resize(numVar,0.0);
		isFeasible = false;
		dVel.resize(numConstr,0.0);
		ub = INF;
		lb = -INF;
		pVel.resize(numVar,0.0);
    }
    Problem::Problem(int nVar,int nConstr) :
		psize(nVar), dsize(nConstr),
		dualLB(dsize,-INF), dualUB(dsize,INF),
		nIter(-1)
    { best.ub = INF;best.lb=-INF;best.isFeasible=false;
      best.perturb.resize(nVar,0.0); best.dual.resize(nConstr,0.0);}


    // Main method
    void Problem::solve(UserHooks &hooks)
    {
		const int lbCheckFreq = 10;
		initialise(hooks);		// set up the swarm
		//------------- initial solution for swarm --------------
		if( param.printLevel > 1) printf("Initialisation:\n");
		DblVec bestLB(param.nParticles);
		IntVec bestIter(param.nParticles,0);
		std::vector<Uniform>  rand(param.nParticles);  // 
		if(param.randomSeed != 0)  rand[0].seed(param.randomSeed);
		else rand[0].seedTime();
		Status status = OK;
#	pragma omp parallel for schedule(static,std::max(1,param.nParticles/param.nCPU))
		for(int idx =0;idx < param.nParticles;++idx){
			ParticleIter p(swarm,idx);
			if(p.idx > 0){ // create random see from previous generator
				rand[p.idx].seed((uint32_t) rand[p.idx-1](
									 0,std::numeric_limits<uint32_t>::max()));
			}
			for(int i=0;i<dsize;++i){ // check initial point is valid
				if( p->dual[i] < dualLB[i] ){
					if(param.printLevel)
						printf("ERROR: initial dual[%d] %.2f below LB %.2f\n",
							   i,p->dual[i],dualLB[i]);
					p->dual[i] = dualLB[i];
				}
				if(p->dual[i] > dualUB[i] ){
					if(param.printLevel)
						printf("ERROR: initial dual[%d] %.2f below UB %.2f\n",
							   i,p->dual[i],dualUB[i]);
					p->dual[i] = dualUB[i];
				}
			}
			       
			if( hooks.reducedCost(*p,p->rc) == ABORT){ status = ABORT;}
			p->rc += p->perturb; 
			if( hooks.solveSubproblem(*p) == ABORT ){status = ABORT; }
			bestLB[p.idx] = p->lb;
			if( param.printLevel > 1){
				printf("\tp%02d: LB=%g UB=%g feas=%d viol=%g -- %g\n",
					   p.idx,p->lb,p->ub,p->isFeasible,
					   p->viol.min(),p->viol.max());
				if(param.printLevel > 2){
					printf("\t\tRedCst %g - %g\n",p->rc.min(),p->rc.max());
					printf("\t\tdVel %g - %g\n",p->dVel.min(),p->dVel.max());
					printf("\t\tpVel %g - %g\n",p->pVel.min(),p->pVel.max());
					printf("\t\tdual %g - %g\n",p->dual.min(),p->dual.max());
					printf("\t\tpert %g - %g\n",
						   p->perturb.min(),p->perturb.max());
				}
			}
			p->pVel = 0;
			p->dVel = 0;
		}
		best.x.resize(psize,0);
		updateBest(hooks);
		if(status == ABORT) return;
		DblVec subgradFactor(param.nParticles,param.subgradFactor);
		DblVec perturbFactor(param.nParticles,param.perturbFactor);
		for(int i=0;i<param.nParticles;++i)
			perturbFactor[i] = i*param.perturbFactor/param.nParticles;
		if( param.printLevel > 1)
			printf("Initial LB=%g UB=%g\n",best.lb,best.ub);
		int noImproveIter = 0,nReset=0,maxNoImprove=1000;
		for(nIter=1; nIter < param.maxIter && cpuTime() < param.maxCPU &&
				  wallTime() < param.maxWallTime && status == OK &&
				  best.lb + param.absGap <= best.ub &&		
				  fabs(best.ub-best.lb/best.ub) > param.relGap;
			++nIter){
			if(param.printLevel > 1 && nIter % param.printFreq == 0)
				printf("Iteration %d --------------------\n",nIter);
#	    pragma omp parallel for schedule(static,std::max(1,param.nParticles/param.nCPU))
			for(int idx =0;idx < param.nParticles;++idx){
				ParticleIter p(swarm,idx);
				// --------- calculate step size/direction ----------
				double norm=0;
				for(int i=0;i<dsize;++i) norm += p->viol[i]*p->viol[i];
				DblVec perturbDir(psize,0.0);
				perturbationDirection(hooks,*p,perturbDir);
				double randAscent = rand[p.idx]() ;  // 0,1 uniform
				double randGlobal = rand[p.idx]() * param.globalFactor;
				const double stepSize = rand[p.idx](
					subgradFactor[p.idx]/2.0,subgradFactor[p.idx]) *
										(best.ub-bestLB[p.idx])  / norm;
// 		const double stepSize = randAscent * param.subgradFactor * 
// 					(best.ub-best.lb)  / norm;
				//-------------- update velocity (step) --------------
				for(int i=0;i<dsize;++i)
					p->dVel[i] = param.velocityFactor * p->dVel[i] +
								 stepSize * p->viol[i] +
								 randGlobal * (best.dual[i] - p->dual[i]);
		
				for(int i=0;i<psize;++i)
					p->pVel[i] = param.velocityFactor * p->pVel[i] +
								 randAscent * perturbFactor[p.idx]* //param.perturbFactor *
								 perturbDir[i] +
								 perturbFactor[p.idx]* randGlobal * (1-2*best.x[i])+
								 randGlobal * (best.perturb[i]-p->perturb[i])
								 ;
				//---------- make a step ------------------------------
				for(int i=0;i<dsize;++i){
					p->dual[i] += p->dVel[i];
					p->dual[i] = std::max(dualLB[i],std::min(dualUB[i],p->dual[i]));
				}
				p->perturb *= 0.5;
				p->perturb += p->pVel; // add step in perturbation velocity
				//---------- solve subproblems -------------------------
				if( hooks.reducedCost(*p,p->rc) == ABORT){
					status = ABORT; continue;
				}
				if( nIter % lbCheckFreq != 0){
					p->rc += p->perturb; // if we care more about ub than lb
					if( hooks.solveSubproblem(*p) == ABORT) status = ABORT;
				}else{			// temporarily set pertubation to zero
					DblVec zero(p->perturb.size(),0.0);
					std::swap(p->perturb,zero);
					if( hooks.solveSubproblem(*p) == ABORT) status = ABORT;
					std::swap(p->perturb,zero);	 // swap back
				}
				if(status == ABORT) continue;
				//if( p->isFeasible ){
				//    // make sure our next step reduces perturbation
				//    p->pVel = p->perturb;
				//    p->pVel *= -0.5/param.velocityFactor;
				//}
				//if( p->lb < bestLB[p.idx]) // downplay current velocity
				//    p->dVel *= 0.5;
				if( nIter % lbCheckFreq == 0){
					if( p->lb > bestLB[p.idx] ){
						bestLB[p.idx] = p->lb;
						bestIter[p.idx] = nIter;
						//subgradFactor[p.idx] *= 1.1;
// 			printf("\tsubgradFactor[%d] increased to %.4f (bestLB = %.2f)\n",
// 			       p.idx,subgradFactor[p.idx],bestLB[p.idx]);
					}else { //if( (nIter - bestIter[p.idx]) % 2 == 0 ){
						subgradFactor[p.idx] = // slow down
							std::max(param.subgradFmin,
									 param.subgradFmult*subgradFactor[p.idx]); 
// 			printf("\tsubgradFactor[%d] decreased to %.4f\n",
// 			       p.idx,subgradFactor[p.idx]);
					}
				}
				if( param.heurFreq > 0 && nIter % param.heurFreq == 0
					&& ! p->isFeasible)	// only run heuristic if not alreay feas.
					if(hooks.heuristics(*p) == ABORT){status=ABORT; continue;}
				if( param.printLevel > 1 && nIter % param.printFreq == 0){
					printf("\tp%02d: LB=%g UB=%g feas=%d minViol=%g\n",
						   p.idx,p->lb,p->ub,p->isFeasible,
						   p->viol.min());
					if(param.printLevel > 2){
						printf("\t\tstepSize=%g randGlobal=%g\n",
							   stepSize,randGlobal);
						printf("\t\tRedCst %g - %g\n",p->rc.min(),p->rc.max());
						printf("\t\tdVel %g - %g\n",p->dVel.min(),p->dVel.max());
						printf("\t\tpDir %g - %g\n",perturbDir.min(),perturbDir.max());
						printf("\t\tpVel %g - %g\n",p->pVel.min(),p->pVel.max());
						printf("\t\tdual %g - %g\n",p->dual.min(),p->dual.max());
						printf("\t\tpert %g - %g\n",
							   p->perturb.min(),p->perturb.max());
					}
				}
			    
			}
			if( updateBest(hooks))
				noImproveIter = 0;
			else{
				if(++noImproveIter > maxNoImprove ){ // reset duals
					++nReset;
					swarm[0]->dual = 0;
					double minLB=INF,maxLB=-INF;
					for(int idx=0;idx<param.nParticles;++idx){
						if(swarm[idx]->lb < minLB) minLB = swarm[idx]->lb;
						if(swarm[idx]->lb > maxLB) maxLB = swarm[idx]->lb;
					}
					hooks.reducedCost(*swarm[0],swarm[0]->rc);
					double maxdual = std::max(0.0,swarm[0]->rc.max()) -
									 std::min(0.0,swarm[0]->rc.min());
					if(maxdual < 1e-9) maxdual=1;
					#pragma omp parallel for
					for(int idx=0;idx<param.nParticles;++idx){
						ParticleIter p(swarm,idx);
						for(int j=0;j<dsize;++j){// generate pertubation about
							if(p->dual[j] == 0 && rand[idx]() < 2.0/dsize){
								p->dual[j] = rand[idx](
									std::max(dualLB[j],-maxdual),
									std::min(dualUB[j],maxdual));
							}else
								p->dual[j]= best.dual[j]*rand[idx](
									1-0.2/nReset,1 + 0.2/nReset);
							p->dual[j] = std::max(
								dualLB[j],std::min(dualUB[j],p->dual[j]));
						}
						p->perturb = 0.0;
						p->dVel = 0.0;
						p->pVel = 0.0;
						for(int i=0;i<psize;++i) // discourage current best
							if( rand[idx](0,1) <0.1)  // try to get diversity
								p->pVel[i] = rand[idx](0,1) *perturbFactor[idx]
											 *maxdual*(2*best.x[i]-1);
						subgradFactor[p.idx] = param.subgradFactor/nReset;
						perturbFactor[p.idx] *= 1.1;
						if(hooks.reducedCost(*p,p->rc)==ABORT){status=ABORT;}
						if(hooks.solveSubproblem(*p) == ABORT){status=ABORT;}
					}
					param.heurFreq=std::max(param.heurFreq/2,1);
					param.subgradFmin *= 0.5;
					updateBest(hooks);
					noImproveIter = 0;
					maxNoImprove += 1000;
					if(param.printLevel > 0){
						printf("%2d: Reset swarm LBs were %.2f - %.2f now",
							   nIter,minLB,maxLB);
						minLB=INF;maxLB=-INF;
						for(int idx=0;idx<param.nParticles;++idx){
							if(swarm[idx]->lb < minLB) minLB = swarm[idx]->lb;
							if(swarm[idx]->lb > maxLB) maxLB = swarm[idx]->lb;
						}
						printf(" %.2f - %.2f\n",minLB,maxLB);
					}
				}
			}
			if( param.printLevel && nIter % param.printFreq == 0)
				printf("%2d: LB=%.2f UB=%.2f\n",nIter,best.lb,best.ub);
		}
    }

    double Problem::swarmRadius() const
    { return 0; }		// not yet implemented

    double Problem::wallTime() const
    {
		return omp_get_wtime() - _wallTime;
    }
    // set up an initial swarm
    void Problem::initialise(UserHooks &hooks)
    {
		nIter = 0;
		omp_set_num_threads(param.nCPU);
		_wallTime = omp_get_wtime();
		timer.reset();
		if( swarm.size() >= (size_t)param.nParticles )
			return;		// all done
		// figure out what range of duals makes sense
		if(best.dual.size() < (size_t)dsize) best.dual.resize(dsize,0.0);
		if(best.perturb.size() < (size_t)psize) best.perturb.resize(psize,0.0);
		if(best.rc.size() < (size_t)psize) best.rc.resize(psize,0.0);
		hooks.reducedCost(best,best.rc);
		// find biggest absolute value reduced cost
		const double maxCost = std::max(
			*std::max_element(best.rc.begin(),best.rc.end()),
			- *std::min_element(best.rc.begin(),best.rc.end()));
		swarm.reserve(param.nParticles);
		Uniform rand;
		if(param.randomSeed == 0) rand.seedTime();
		else rand.seed(param.randomSeed);
		while(swarm.size() < (size_t)param.nParticles){
			Particle *p = new Particle(psize,dsize);
			swarm.push_back(p);
			for(int j=0;j<dsize;++j) // generate values up to maxCost
				p->dual[j] = rand(std::max(dualLB[j],-maxCost),
								  std::min(dualUB[j],maxCost));
			p->perturb = 0.0;
			p->dVel = 0.0;
			p->pVel = 0.0;
		}	
    }
    bool Problem::updateBest(UserHooks &hooks)
    {   Particle *bestP = 0;
		bool improved = false;
		for(ParticleIter p(swarm); ! p.done(); ++p){
			if(p->isFeasible && p->ub < best.ub){
				best.ub = p->ub;
				best.isFeasible = true;
				best.x = p->x;
				best.perturb = p->perturb; // perturbation gives good feasible
				bestP = &(*p);
			}
			if(p->lb > best.lb){
				best.lb = p->lb;
				best.dual = p->dual;
				best.viol = p->viol;
				improved = true;
			}
		}
		if(bestP) // only call once if current swarm improved optimum
			hooks.updateBest(*bestP);
		return improved || (bestP != 0);
    }


    /** Calculate a direction for the perturbation that is likely
		to lead to more feasible solutions (using hooks.fixConstraint) */
    void Problem::perturbationDirection(UserHooks &hooks,const Particle &p,
										DblVec &dir) const
    {   const double eps = param.eps;
		dir = 0.0;		// set all to zero
		for(int i=0;i<dsize;++i){
			if((p.viol[i] > eps && dualLB[i] < -eps) ||
			   (p.viol[i] < -eps && dualUB[i] > eps)) { // constraint violated
				SparseVec feas;
				hooks.fixConstraint(i,p,feas);
				for(SparseVec::iterator x=feas.begin();x!=feas.end();++x){
					if(x->second >= 1) dir[x->first] -= 1;
					if(x->second <= 0) dir[x->first] += 1;
				}
			}
		
		}
    }


    
};				// end namespace
	
/* Stuff for emacs/xemacs to conform to the "Visual Studio formatting standard"
 * Local Variables:
 * tab-width: 4
 * eval: (c-set-style "stroustrup")
 * End:
*/
