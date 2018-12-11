// CpuTimer.h: interface for the CpuTimer class.
//             Roughly equivalent to common::ComputationTime module
//             but does not depend on the rest of the common library
//
//////////////////////////////////////////////////////////////////////

#ifndef __CPM_CPUTIMER_H__
#define __CPM_CPUTIMER_H__


#if defined(_M_X64) || defined(_M_IA64) || defined(_M_AMD64) // Windows (x64)
#define CPUTIMER_USE_CPU64
#include <Windows.h>
#ifdef min  // don't let windows redefine min & max!
#undef min
#undef max
#endif
#elif defined(_MSC_VER)	// Windows (x86)
#define CPM_CPUTIMER_USE_CLOCK
//#define CPUTIMER_USE_PROCESSTIMES
//#include <Windows.h>
#include <time.h>
#else // POSIX systems
#include <sys/resource.h>
#endif
#include <ctime>
#include <iostream>


class CpuTimer
{
public:
    //CpuTimer();
    CpuTimer(double maxCpuTimer=3600,double maxWallTime=-1);
    CpuTimer(const CpuTimer& other);
    virtual ~CpuTimer();

    double elapsedSeconds() const;
	double elapsedWallTime() const { return difftime(time(0),startTime_t);}
    double timeRemaining() const;
    void setTimeRemaining(double timeRemaining,double maxWallTime=-1);
	void reset(double timeRemaining=0.0,double maxWallTime=-1) {
		setTimeRemaining(timeRemaining,maxWallTime); }
    void setTimeLimit(double timeRemaining,double maxWallTime=-1) ;
    /// like setTimeRemaining but doesn't reset elapsed. Just limits remaining
    /// CPU &/or wall-time
    bool timeLimitReached() const;
    std::ostream& print(std::ostream& os) const;
protected:
	time_t startTime_t;
	double wallTimeLimit;
#   if defined(CPUTIMER_USE_CPU64)
	FILETIME startProcKernel, startProcUser;
	double elapsedSeconds64Bit() const;
#   elif defined(CPM_CPUTIMER_USE_CLOCK)
    clock_t startTime;
#   elif defined(CPUTIMER_USE_PROCESSTIMES)
	FILETIME startTime;
#   else
    struct rusage startTime;
#   endif
    double endTime;// clock_t doesn't hold enough seconds
};

std::ostream& operator<<(std::ostream& os, const CpuTimer& a);
#endif

