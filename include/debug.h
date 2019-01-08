#ifndef __DEBUG_H__
#define __DEBUG_H__

#define GENERAL_MESSAGES 1
#define VERBOSE 0
#define WORKING_DEBUG_1 0
#define WORKING_DEBUG_2 0
#define DEF_PAUSE 0


// self explanatory
#if GENERAL_MESSAGES
#  define PE(x) std::cout << x << std::endl;
#  define PF(x) std::cout << x << std::flush;
#  define CHECK(x) std::cout << #x << ": " << x << std::endl;
#else
#  define PE(x)
#  define PF(x)
#endif

// for preprocessing and reconstructing
#if VERBOSE 
#  define VE(x) std::cout << x << std::endl;
#  define VF(x) std::cout << x << std::flush;
#else
#  define VE(x)
#  define VF(x)
#endif

// for preprocessing and reconstructing
#if WORKING_DEBUG_1 
#  define PF1(x) std::cout << x << std::flush;
#  define PE1(x) std::cout << x << std::endl;
#else
#  define PE1(x)
#  define PF1(x)
#endif

//
#if WORKING_DEBUG_2 
#  define PF2(x) std::cout << x << std::flush;
#  define PE2(x) std::cout << x << std::endl;
#else
#  define PE2(x)
#  define PF2(x)
#endif

#if DEF_PAUSE
#  define PAUSE std::cin.get();
#else
#  define PAUSE 
#endif


#endif