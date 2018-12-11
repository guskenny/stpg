// CpuTimer.cpp: implementation of the CpuTimer class.
//
//////////////////////////////////////////////////////////////////////

#include "CpuTimer.h"

#include <limits>
#include <math.h>


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CpuTimer::CpuTimer(double maxCpuTimer,double maxWallTime)
{
     // Default length is 1 hour
     reset(maxCpuTimer,maxWallTime);
}

#ifdef CPUTIMER_USE_CPU64
CpuTimer::CpuTimer(const CpuTimer &other) :
	startProcKernel (other.startProcKernel),
	startProcUser (other.startProcUser),
	endTime (other.endTime)
{
}
#else
CpuTimer::CpuTimer(const CpuTimer& other) :
    startTime (other.startTime),
    endTime (other.endTime)
{
}
#endif

CpuTimer::~CpuTimer()
{

}

double CpuTimer::elapsedSeconds() const
{
# if defined(CPUTIMER_USE_CPU64)
  return elapsedSeconds64Bit();
# elif defined(CPM_CPUTIMER_USE_CLOCK)
  const clock_t t = clock();
  if(t < startTime){
    const double clock_max = pow (2.0, 8.0*sizeof(clock_t) );    
    return ((double)t - (double) startTime  + clock_max)  / CLOCKS_PER_SEC;
  }
  return (double)(t - startTime) / CLOCKS_PER_SEC;
# elif defined(CPUTIMER_USE_PROCESSTIMES)
  FILETIME cpu,sys,create,exit;
  GetProcessTimes(GetCurrentProcess(),&create,&exit,&sys,&cpu);
  return (cpu-startTime)/1.0e7;
# else // POSIX style resource usage is more accurate
  struct rusage ru;
  getrusage(RUSAGE_SELF,&ru);
  ru.ru_utime.tv_sec -= startTime.ru_utime.tv_sec;
  if(ru.ru_utime.tv_usec < startTime.ru_utime.tv_usec){
        ru.ru_utime.tv_sec--;
        ru.ru_utime.tv_usec += 1000000 - startTime.ru_utime.tv_usec;
  }else 
        ru.ru_utime.tv_usec -=  startTime.ru_utime.tv_usec;
  return (ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec*1e-6);
# endif
}
    
double CpuTimer::timeRemaining() const
{
# if defined(CPUTIMER_USE_CPU64)
	return endTime - elapsedSeconds();
# elif defined(CPM_CPUTIMER_USE_CLOCK)
  const clock_t t = clock();
  if(t < startTime){		// laped counter
    const double clock_max = pow (2.0, 8.0*sizeof(clock_t) );
    return (endTime - (double)t	- clock_max) / CLOCKS_PER_SEC;
  }
  return (double)(endTime - t) / CLOCKS_PER_SEC;
# else
  return endTime - elapsedSeconds();
# endif
}
    
void CpuTimer::setTimeRemaining(double timeRemaining,double maxWallTime)
{
	startTime_t = time(0);
	wallTimeLimit = maxWallTime;
# if defined(CPUTIMER_USE_CPU64)
	FILETIME procCreation, procExit;
	GetProcessTimes(GetCurrentProcess(), &procCreation, &procExit, &startProcKernel, &startProcUser);
	endTime = timeRemaining;
#   elif defined(CPM_CPUTIMER_USE_CLOCK)
    startTime = clock();
    endTime  = (double)startTime + (timeRemaining * CLOCKS_PER_SEC);
#   elif defined(CPUTIMER_USE_PROCESSTIMES)
	FILETIME cpu,sys,create,exit;
	GetProcessTimes(GetCurrentProcess(),&create,&exit,&sys,&cpu);
	startTime = cpu;
	endTime = timeRemaining;
#   else
    getrusage(RUSAGE_SELF,&startTime);
    endTime = timeRemaining;
#   endif
}

void CpuTimer::setTimeLimit(double timeRemaining,double maxWallTime)
{
	if(maxWallTime < 0) wallTimeLimit = maxWallTime; // infinite
	else wallTimeLimit = elapsedWallTime() + maxWallTime;
# if defined(CPUTIMER_USE_CPU64)
	endTime = elapsedSeconds() + timeRemaining;
# elif defined(CPM_CPUTIMER_USE_CLOCK)
  const clock_t t = clock();
  endTime = (double)t + (timeRemaining * CLOCKS_PER_SEC);
  if(t < startTime){ // laped counter
    const double clock_max = pow (2.0, 8.0*sizeof(clock_t) );
    endTime  += clock_max;
  }
# else
  endTime = elapsedSeconds() + timeRemaining;
# endif
}

bool CpuTimer::timeLimitReached() const
{
	if(wallTimeLimit >=0 && elapsedWallTime() > wallTimeLimit) return true;
# if defined(CPUTIMER_USE_CPU64)
	return elapsedSeconds() >= endTime;
# elif defined(CPM_CPUTIMER_USE_CLOCK)
  const clock_t t = clock();
  if(t < startTime){
    const double clock_max = pow (2.0, 8.0*sizeof(clock_t) );    
    return ((double)t +  clock_max) >= endTime;
  }
  return double(t) >= endTime;
# else
  return elapsedSeconds() >= endTime;
# endif
}

std::ostream& CpuTimer::print(std::ostream& os) const
{
    os << "Time elapsed: " << elapsedSeconds() <<
        " remaining: " << timeRemaining();
    return os;
}

std::ostream& operator<<(std::ostream& os, const CpuTimer& a) 
{
    return a.print (os);
}

#ifdef CPUTIMER_USE_CPU64
double CpuTimer::elapsedSeconds64Bit() const
{
	// Current process time (kernel + user)
	FILETIME procCreation, procExit, procKernel, procUser;
	GetProcessTimes(GetCurrentProcess(), &procCreation, &procExit, &procKernel, &procUser);
	// Differences between previous times for kernel + user
	LARGE_INTEGER la, lb;
	la.LowPart = procKernel.dwLowDateTime;
	la.HighPart = procKernel.dwHighDateTime;
	lb.LowPart = startProcKernel.dwLowDateTime;
	lb.HighPart = startProcKernel.dwHighDateTime;
	ULONGLONG procKernelDiff = la.QuadPart - lb.QuadPart;
	la.LowPart = procUser.dwLowDateTime;
	la.HighPart = procUser.dwHighDateTime;
	lb.LowPart = startProcUser.dwLowDateTime;
	lb.HighPart = startProcUser.dwHighDateTime;
	ULONGLONG procUserDiff = la.QuadPart - lb.QuadPart;
	// Sum differences to get total CPU time
	ULONGLONG procTotal = procKernelDiff + procUserDiff;
	double eTime = procTotal * 1.0e-7;
	return eTime;
}
#endif

// Emacs settings to make indentation match with Visual C++
// Local Variables:
// c-file-style: "stroustrup"
// End:
