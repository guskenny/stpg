CXXFLAGS= -O4 -Wall -fopenmp#-pg
LDFLAGS=#-pg

OBJ= LaPSO.o anyoption.o CpuTimer.o VolVolume.o

build: dcmst # $(TEST) 

dcmst: dcmst.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o dcmst dcmst.o $(OBJ)

dcmst.o: dcmst.cpp LaPSO.hpp VolVolume.hpp

LaPSO.o: LaPSO.cpp LaPSO.hpp anyoption.h

clean:
	-rm -f *.o 


depend:
	makedepend -Y $(CPPFLAGS) *.cpp >& /dev/null

# DO NOT DELETE

CpuTimer.o: CpuTimer.h
LaPSO.o: LaPSO.hpp Random.h anyoption.h
Random.o: Random.h
anyoption.o: anyoption.h
dcmst.o: LaPSO.hpp CpuTimer.h
