#all: pcpsp QOLpcpsp cumulative
all: merge_stpg

LR = LagrangianHeuristic
# INC = /usr/local/gurobi/6.5.1/include/
# LIB = /usr/local/gurobi/6.5.1/lib/

# dependency generation
DEPDIR := .d
$(shell mkdir -p $(DEPDIR)/source >/dev/null)
$(shell mkdir -p $(DEPDIR)/$(LR) >/dev/null)
$(shell mkdir -p $(DEPDIR)/QOL/lagrangian >/dev/null)
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.Td
POSTCOMPILE = mv -f $(DEPDIR)/$*.Td $(DEPDIR)/$*.d


CXXFLAGS= -fno-stack-protector -g $(DEPFLAGS)
# choose flags for debugging vs best performance
#CXXFLAGS+=-D_GLIBCXX_DEBUG
CXXFLAGS+=-O3 -DNDEBUG

include QOL/Makefile
CPPFLAGS+= -I$(CPLEX_ROOT)/concert/include
CXXFLAGS+=-fopenmp
CPPFLAGS+=-DIL_STD  -Wno-sign-compare
# CPPFLAGS+=-lncurses
#QOL_LIBS:= -L$(wildcard $(CPLEX_ROOT)/concert/lib/*/static_pic) -lconcert -lilocplex $(QOL_LIBS) #-lgurobi_c++ -lgurobi$(GUROBI_VER)
QOL_LIBS:= -L/opt/ibm/ILOG/CPLEX_Studio127/cplex/lib/x86-64_linux/static_pic -lilocplex -lconcert -lcplex -lm -lpthread -L/opt/ibm/ILOG/CPLEX_Studio127/concert/lib/x86-64_linux/static_pic -lilocplex -lconcert -lcplex -lm -lpthread $(QOL_LIBS) #-lgurobi_c++ -lgurobi$(GUROBI_VER)
CPPFLAGS+= -I./include -I. -I./$(LR)
$(QOL_LIB): $(wildcard $(QOL_DIR)/*.h $(QOL_DIR)/*.cpp)
	cd $(QOL_DIR); make

# LAPSO_SRC=source/Particle.cpp $(wildcard QOL/lagrangian/[^d]*.cpp)
CPPFLAGS+= -I QOL/lagrangian
#anyoption.cpp QOL/lagrangian/VolVolume.cpp QOL/lagrangian/CpuTime
SRC=source/u_graph.cpp source/STPGModel.cpp source/util.cpp source/STPGSolver.cpp source/LocalSearch.cpp source/SolutionMerger.cpp source/STPGMip.cpp source/STPGMergeMip.cpp

# SRC=source/daten.cpp source/network.cpp source/boostMaxFlow.cpp source/CumulativeModel.cpp source/util.cpp source/MaxClosure_PP.cpp $(wildcard source/MaxClosure_B*.cpp) source/MaxClosureFactory.cpp source/UpitSolver.cpp source/graph.cpp source/SinglePSolver.cpp source/Preprocess.cpp source/SinglePModel.cpp source/SolutionMerger.cpp source/MergeSolver.cpp source/MergeSolverSimple.cpp source/MergeSolverCompact.cpp source/RandomSearch.cpp source/LocalSearch.cpp source/ConeMiner.cpp #$(LAPSO_SRC) #$(LR)/solver.cpp $(LR)/mip.cpp $(LR)/ACO_solution.cpp $(LR)/Timer.cc $(LR)/Random.cc $(LR)/MaxClosure_NetworkFlow_LR.cpp $(LR)/pheromones.cpp source/lr_graph.cpp $(LR)/solver_functions.cpp source/solver_interface.cpp source/MaxClosure_NetworkFlow.cpp $(LAPSO_SRC)
OBJ=$(patsubst %.cpp,%.o,$(SRC))

LIBS = -L${LIB}$ -lgurobi_c++ -lpthread -lm -lgurobi65

%.o: %.cpp $(DEPDIR)/%.d
.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@; $(POSTCOMPILE);

merge_stpg: $(OBJ) $(QOL_LIB) source/merge_stpg.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -o $@ $^ $(QOL_LIBS)

clean::
	-rm -f source/*.o QOL/lagrangian/*.o $(LR)/*.o

$(DEPDIR)/%.d: ;
.PRECIOUS: $(DEPDIR)/%.d

-include $(patsubst %,$(DEPDIR)/%.d,$(basename $(SRC)))
