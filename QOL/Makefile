#!/usr/bin/make
#all:libqol.a
# compile into a single library: libqol.a

CXXFLAGS+=-fPIC -Wall -Wno-sign-compare -Wno-reorder -fpermissive

QOL_SRC:=$(wildcard *.cpp)

QOL_LIBS:=
UNAME=$(shell uname)
QOL_DIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))

include $(QOL_DIR)/QOL.mk

GCC_VER=$(shell gcc -dumpversion)

QOL_OBJ:=$(patsubst %.cpp,%.o,$(QOL_SRC))
CPPFLAGS+= -I$(QOL_DIR) #-I$(QOL_DIR)resources
QOL_LIBS:= $(QOL_DIR)/libqol.a $(QOL_LIBS) -lpthread -lm
QOL_LIB:= $(QOL_DIR)/libqol.a

debug: CXXFLAGS+=-DDEBUG -O0 -g
debug: CFLAGS+=-DDEBUG -O0 -g
debug: libqol.a


release: CXXFLAGS+=-O3
release: CFLAGS+=-O3
release: libqol.a

libqol.a: $(QOL_OBJ)
	ar rs libqol.a $(QOL_OBJ) ## 	ranlib libqol.a 

clean::
	-rm -f *.o libqol.a


depend:
	makedepend -Y *.cpp
# DO NOT DELETE

ACO.o: ACO.h QolOptimiser.h CpuTimer.h QolUtil.h QolRandom.h QolData.h
ACO.o: QolModel.h QolHeuristic.h QolSolutionBase.h QolMIPSolver.h QolMIP.h
ACO.o: QolSolutionPrimal.h QolSolutionFull.h QolParams.h QolColFormulation.h
ACO.o: QolMIPInstance.h LagrangianModel.h VolumeModel.h BendersModel.h
ACO.o: QolLog.h
BendersModel.o: BendersModel.h QolModel.h QolUtil.h QolRandom.h QolData.h
BendersModel.o: QolHeuristic.h QolSolutionBase.h QolMIPSolver.h QolMIP.h
BendersModel.o: QolSolutionPrimal.h QolSolutionFull.h QolParams.h
BendersModel.o: QolColFormulation.h QolMIPInstance.h CpuTimer.h
CoinFormulation.o: CoinFormulation.h QolMIPSolver.h QolMIP.h QolUtil.h
CoinFormulation.o: QolRandom.h QolData.h QolSolutionPrimal.h
CoinFormulation.o: QolSolutionBase.h QolSolutionFull.h QolParams.h
CoinFormulation.o: QolColFormulation.h QolMIPInstance.h
CplexFormulation.o: CplexFormulation.h QolMIPSolver.h QolMIP.h QolUtil.h
CplexFormulation.o: QolRandom.h QolData.h QolSolutionPrimal.h
CplexFormulation.o: QolSolutionBase.h QolSolutionFull.h QolParams.h
CplexFormulation.o: QolColFormulation.h QolMIPInstance.h
CpuTimer.o: CpuTimer.h
GurobiFormulation.o: GurobiFormulation.h QolMIPSolver.h QolMIP.h QolUtil.h
GurobiFormulation.o: QolRandom.h QolData.h QolSolutionPrimal.h
GurobiFormulation.o: QolSolutionBase.h QolSolutionFull.h QolParams.h
GurobiFormulation.o: QolColFormulation.h QolMIPInstance.h
LagrangianModel.o: LagrangianModel.h QolModel.h QolUtil.h QolRandom.h
LagrangianModel.o: QolData.h QolHeuristic.h QolSolutionBase.h QolMIPSolver.h
LagrangianModel.o: QolMIP.h QolSolutionPrimal.h QolSolutionFull.h QolParams.h
LagrangianModel.o: QolColFormulation.h QolMIPInstance.h CpuTimer.h
PSO.o: PSO.h QolOptimiser.h CpuTimer.h QolUtil.h QolRandom.h QolData.h
PSO.o: QolModel.h QolHeuristic.h QolSolutionBase.h QolMIPSolver.h QolMIP.h
PSO.o: QolSolutionPrimal.h QolSolutionFull.h QolParams.h QolColFormulation.h
PSO.o: QolMIPInstance.h LagrangianModel.h VolumeModel.h BendersModel.h
PSO.o: QolLog.h
QolCommData.o: QolCommData.h QolData.h QolUtil.h QolRandom.h
QolMIP.o: QolMIP.h QolUtil.h QolRandom.h QolData.h
QolMIPSolver.o: QolMIPSolver.h QolMIP.h QolUtil.h QolRandom.h QolData.h
QolMIPSolver.o: QolSolutionPrimal.h QolSolutionBase.h QolSolutionFull.h
QolMIPSolver.o: QolParams.h QolColFormulation.h QolMIPInstance.h
QolOptimiser.o: QolOptimiser.h CpuTimer.h QolUtil.h QolRandom.h QolData.h
QolOptimiser.o: QolModel.h QolHeuristic.h QolSolutionBase.h QolMIPSolver.h
QolOptimiser.o: QolMIP.h QolSolutionPrimal.h QolSolutionFull.h QolParams.h
QolOptimiser.o: QolColFormulation.h QolMIPInstance.h LagrangianModel.h
QolOptimiser.o: VolumeModel.h BendersModel.h QolLog.h
VolumeModel.o: VolumeModel.h QolModel.h QolUtil.h QolRandom.h QolData.h
VolumeModel.o: QolHeuristic.h QolSolutionBase.h QolMIPSolver.h QolMIP.h
VolumeModel.o: QolSolutionPrimal.h QolSolutionFull.h QolParams.h
VolumeModel.o: QolColFormulation.h QolMIPInstance.h CpuTimer.h
