#ifndef __SETTINGSHANDLER_H__
#define __SETTINGSHANDLER_H__

#include <debug.h>
#include <fstream>
#include <iostream>

class SettingsHandler {

  public:
    std::string SETTINGS_FILE = "settings.ini";
    bool PREPROCESS=1;
    bool INIT_SOL_TYPE=0;
    int  NUM_MOVES=1000;
    int  JUMP_FREQ=0;
    int  NUM_JUMPS=1;
    int  DEAD_MOVES=1;
    bool RECORD_DATA=1;
    bool TEST_GROUPS=0;
    bool WARM_START=false;
    size_t MIP_TIME=0;
    size_t NUM_ITER=1;

    size_t TIMER_LENGTH=0;
    int RANDOM_SEARCH=0;
    int NUM_TOTAL_RUNS=1;
    bool QUIET=0;
    bool TEMP_FLAG;
    bool TEST_SOL_MERGE=0;
    bool RECORD_RUN=1;
    bool RECORD_PASS_TIME=0;
    double SCALED_RESOURCE=1.0;
    double MIP_GAP=0.0;
    size_t FIND_DEPTH_TRIALS=1000;
    size_t BOUND_DEPTH=1;
    double SA_ALPHA=0.99;
    double SA_T_MIN=0.00001;
    int SA_ITER=100;
    int FULL_RUNS=10;
    bool GREEDY_PLUS=0;
    bool LOAD_SEEDS=1;
    size_t SOL_IDX=0;
    bool AUTO_REPAIR=0;
    size_t INIT_SEEDS=1;
    size_t NUM_SEEDS=1;
    size_t MERGE_TYPE=1;
    size_t ITER_INCR=0;
    size_t FIX_BEST_GROUP=false;
    bool GROUP_MERGE=false;
    size_t MERGE_POP_SIZE=4;
    size_t SWAP_POP=1;
    size_t MAX_SUB_SWAPS=1000;
    size_t WINDOW_SIZE=2;
    size_t VNS_WINDOW_LIMIT=0;
    double REDUCED_BLOCKS=0.0;
    double MODIFIER_WEIGHT=0.0;
    double PROCESS_THRESHOLD=0.0;
    double MIN_DIFF_TO_PROCESS=0.0;
    bool ONLY_NEGATIVE=false;
    size_t WINDOW_SEARCH_TIME = 0;
    double TOP_PERCENT = 1.0;
    double LIMIT_MULTIPLIER = 1.25;
    double P_SELECTION = 1.0;
    int NUM_SUB_PERIODS = 1;
    double RES_CUTOFF = 0.3;
    double OBJ_UB = 9e20;
    size_t WALL_SEARCH_TIME = 0;
    size_t CPU_SEARCH_TIME = 0;
    int CONVERGE_STOP=0;

    SettingsHandler(){

      std::ifstream in(SETTINGS_FILE);

      std::string in_line;
      std::string delimiter = "=";

      while ((std::getline(in, in_line))) {

        switch (in_line[0]) {
          case '#':       // skip comment lines
            continue;
          case '\n':      // skip empty lines
            continue;
          case '\0':      // skip empty lines at end of file
            continue;
        }
        int pos = in_line.find(delimiter);
        std::string value = in_line.substr(pos+1,in_line.length());
        std::string setting = in_line.substr(0, pos);

        if (setting.compare("NUM_MOVES") == 0){
          NUM_MOVES = stoi(value);
          continue;
        }
        else if (setting.compare("JUMP_FREQ") == 0){
          JUMP_FREQ = stoi(value);
          continue;
        }
        else if (setting.compare("NUM_JUMPS") == 0){
          NUM_JUMPS = stoi(value);
          continue;
        }
        else if (setting.compare("DEAD_MOVES") == 0){
          DEAD_MOVES = stoi(value);
          continue;
        }
        else if (setting.compare("TEST_GROUPS") == 0){
          TEST_GROUPS = stoi(value);
          continue;
        }
        else if (setting.compare("PREPROCESS") == 0){
          PREPROCESS = stoi(value);
          continue;
        }
        else if (setting.compare("MIP_TIME") == 0){
          MIP_TIME = stoi(value);
          continue;
        }
        else if (setting.compare("INIT_SOL_TYPE") == 0){
          INIT_SOL_TYPE= stoi(value);
          continue;
        }
        else if (setting.compare("NUM_TOTAL_RUNS") == 0){
          NUM_TOTAL_RUNS = stoi(value);
          continue;
        }
        else if (setting.compare("RECORD_DATA") == 0){
          RECORD_DATA = stoi(value);
          continue;
        }
        else if (setting.compare("QUIET") == 0){
          QUIET = stoi(value);
          continue;
        }
        else if (setting.compare("TEMP_FLAG") == 0){
          TEMP_FLAG = stoi(value);
          continue;
        }
        else if (setting.compare("TEST_SOL_MERGE") == 0){
          TEST_SOL_MERGE = stoi(value);
          continue;
        }
				else if (setting.compare("RECORD_RUN") == 0){
          RECORD_RUN = stoi(value);
          continue;
        }
				else if (setting.compare("RECORD_PASS_TIME") == 0){
          RECORD_PASS_TIME = stoi(value);
          continue;
        }
        else if (setting.compare("RANDOM_SEARCH") == 0){
          RANDOM_SEARCH = stoi(value);
          continue;
        }
        else if (setting.compare("SCALED_RESOURCE") == 0){
          SCALED_RESOURCE = stof(value);
          continue;
        }
        else if (setting.compare("MIP_GAP") == 0){
          MIP_GAP = stof(value);
          continue;
        }
        else if (setting.compare("FIND_DEPTH_TRIALS") == 0){
          FIND_DEPTH_TRIALS = stoi(value);
          continue;
        }
        else if (setting.compare("BOUND_DEPTH") == 0){
          BOUND_DEPTH = stoi(value);
          continue;
        }
        else if (setting.compare("SA_ALPHA") == 0){
          SA_ALPHA = stof(value);
          continue;
        }
        else if (setting.compare("SA_T_MIN") == 0){
          SA_T_MIN = stof(value);
          continue;
        }
        else if (setting.compare("SA_ITER") == 0){
          SA_ITER = stoi(value);
          continue;
        }
        else if (setting.compare("FULL_RUNS") == 0){
          FULL_RUNS = stoi(value);
          continue;
        }
        else if (setting.compare("GREEDY_PLUS") == 0){
          GREEDY_PLUS = stoi(value);
          continue;
        }
        else if (setting.compare("FIX_BEST_GROUP") == 0){
          FIX_BEST_GROUP = stoi(value);
          continue;
        }
        else if (setting.compare("GROUP_MERGE") == 0){
          GROUP_MERGE = stoi(value);
          continue;
        }
        else if (setting.compare("MERGE_POP_SIZE") == 0){
          MERGE_POP_SIZE = stoi(value);
          continue;
        }
        else if (setting.compare("SWAP_POP") == 0){
          SWAP_POP = stoi(value);
          continue;
        }
        else if (setting.compare("MAX_SUB_SWAPS") == 0){
          MAX_SUB_SWAPS = stoi(value);
          continue;
        }
        else if (setting.compare("LOAD_SEEDS") == 0){
          LOAD_SEEDS = stoi(value);
          continue;
        }
        else if (setting.compare("SOL_IDX") == 0){
          SOL_IDX = stoi(value);
          continue;
        }
        else if (setting.compare("AUTO_REPAIR") == 0){
          AUTO_REPAIR = stoi(value);
          if (AUTO_REPAIR){
            MAX_SUB_SWAPS = 0;
          }
          continue;
        }
        else if (setting.compare("INIT_SEEDS") == 0){
          INIT_SEEDS = stoi(value);
          continue;
        }
        else if (setting.compare("NUM_SEEDS") == 0){
          NUM_SEEDS = stoi(value);
          continue;
        }
        else if (setting.compare("NUM_ITER") == 0){
          NUM_ITER = stoi(value);
          continue;
        }
        else if (setting.compare("MERGE_TYPE") == 0){
          MERGE_TYPE = stoi(value);
          continue;
        }
        else if (setting.compare("ITER_INCR") == 0){
          ITER_INCR = stoi(value);
          continue;
        }
        else if (setting.compare("WINDOW_SIZE") == 0){
          WINDOW_SIZE = stoi(value);
          continue;
        }
        else if (setting.compare("VNS_WINDOW_LIMIT") == 0){
          VNS_WINDOW_LIMIT = stoi(value);
          continue;
        }
        else if (setting.compare("REDUCED_BLOCKS") == 0){
          REDUCED_BLOCKS = stof(value);
          continue;
        }
        else if (setting.compare("MODIFIER_WEIGHT") == 0){
          MODIFIER_WEIGHT = stof(value);
          continue;
        }
        else if (setting.compare("PROCESS_THRESHOLD") == 0){
          PROCESS_THRESHOLD = stof(value);
          continue;
        }
        else if (setting.compare("MIN_DIFF_TO_PROCESS") == 0){
          MIN_DIFF_TO_PROCESS = stof(value);
          continue;
        }
        else if (setting.compare("WARM_START") == 0){
          WARM_START = stoi(value);
          continue;
        }
        else if (setting.compare("ONLY_NEGATIVE") == 0){
          ONLY_NEGATIVE = stoi(value);
          continue;
        }
        else if (setting.compare("TOP_PERCENT") == 0){
          TOP_PERCENT = stof(value);
          continue;
        }
        else if (setting.compare("LIMIT_MULTIPLIER") == 0){
          LIMIT_MULTIPLIER = stof(value);
          continue;
        }
        else if (setting.compare("P_SELECTION") == 0){
          P_SELECTION = stof(value);
          continue;
        }
        else if (setting.compare("NUM_SUB_PERIODS") == 0){
          NUM_SUB_PERIODS = stoi(value);
          continue;
        }
        else if (setting.compare("RES_CUTOFF") == 0){
          RES_CUTOFF = stof(value);
          continue;
        }
        else if (setting.compare("OBJ_UB") == 0){
          OBJ_UB = stof(value);
          continue;
        }
        else if (setting.compare("WALL_SEARCH_TIME") == 0){
          WALL_SEARCH_TIME = stoi(value);
          continue;
        }
        else if (setting.compare("CPU_SEARCH_TIME") == 0){
          CPU_SEARCH_TIME = stoi(value);
          continue;
        }
        else if (setting.compare("CONVERGE_STOP") == 0){
          CONVERGE_STOP = stoi(value);
          continue;
        }
        else{
          std::cout << "Unrecognised setting: '" << setting << "'!" << std::endl;
        }
      }
    }

    void printSettings(){
      std::cout << std::endl
                << "Algorithm settings:" << std::endl
                << "-------------------" << std::endl
                << "RECORD_DATA: ";
      if (RECORD_DATA)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "RECORD_RUN: ";
      if (RECORD_RUN)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "TEMP_FLAG: ";
      if (QUIET)
        std::cout << "true";
      else
        std::cout << "false";
      if (TEMP_FLAG)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "TEST_SOL_MERGE: ";
      if (TEST_SOL_MERGE)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "RECORD_PASS_TIME: ";
      if (RECORD_PASS_TIME)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "RANDOM_SEARCH: ";
      if (RANDOM_SEARCH)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "NUM_TOTAL_RUNS: " << NUM_TOTAL_RUNS
                << std::endl << "SCALED_RESOURCE: " << SCALED_RESOURCE
                << std::endl << "MIP_GAP: " << MIP_GAP
                << std::endl << "FIND_DEPTH_TRIALS: " << FIND_DEPTH_TRIALS
                << std::endl << "BOUND_DEPTH: " << BOUND_DEPTH
                << std::endl << "SA_ALPHA: " << SA_ALPHA
                << std::endl << "SA_T_MIN: " << SA_T_MIN
                << std::endl << "SA_ITER: " << SA_ITER
                << std::endl << "FULL_RUNS: " << FULL_RUNS
                << std::endl << "MERGE_POP_SIZE: " << MERGE_POP_SIZE
                << std::endl << "SWAP_POP: " << SWAP_POP
                << std::endl << "MAX_SUB_SWAPS: " << MAX_SUB_SWAPS
                << std::endl << "GREEDY_PLUS: ";
      if (GREEDY_PLUS)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "PROCESS_THRESHOLD: " << PROCESS_THRESHOLD
                << std::endl << "MIN_DIFF_TO_PROCESS: " << MIN_DIFF_TO_PROCESS
                << std::endl << "TIMER_LENGTH: ";
      if (TIMER_LENGTH > 0)
        std::cout << TIMER_LENGTH;
      else
        std::cout << "NONE";
      std::cout << std::endl << "GROUP_MERGE: ";
        if (GROUP_MERGE)
          std::cout << "true";
        else
          std::cout << "false";
      std::cout << std::endl << "NUM_TOTAL_RUNS: " << NUM_TOTAL_RUNS;
      std::cout << std::endl << "LOAD_SEEDS: ";
        if (LOAD_SEEDS)
          std::cout << "true";
        else
          std::cout << "false";
      std::cout << std::endl << "AUTO_REPAIR: ";
        if (AUTO_REPAIR)
          std::cout << "true";
        else
          std::cout << "false";
      std::cout << std::endl << "MERGE_POP_SIZE: " << MERGE_POP_SIZE << std::endl
                << "FIX_BEST_GROUP: " << FIX_BEST_GROUP << std::endl
                << "INIT_SEEDS: " << INIT_SEEDS << std::endl
                << "NUM_SEEDS: " << NUM_SEEDS << std::endl
                << "NUM_ITER: " << NUM_ITER << std::endl
                << "MERGE_TYPE: " << MERGE_TYPE << std::endl
                << "ITER_INCR: " << ITER_INCR << std::endl
                << "WINDOW_SIZE: " << WINDOW_SIZE << std::endl
                << "VNS_WINDOW_LIMIT: " << VNS_WINDOW_LIMIT << std::endl
                << "REDUCED_BLOCKS: " << REDUCED_BLOCKS << std::endl
                << "MODIFIER_WEIGHT: " << MODIFIER_WEIGHT << std::endl
                << "WARM_START: ";
      if (WARM_START)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "ONLY_NEGATIVE: ";
      if (ONLY_NEGATIVE)
        std::cout << "true";
      else
        std::cout << "false";
      std::cout << std::endl << "WINDOW_SEARCH_TIME: " << WINDOW_SEARCH_TIME
                << std::endl << "TOP_PERCENT: " << TOP_PERCENT
                << std::endl << "LIMIT_MULTIPLIER: " << LIMIT_MULTIPLIER
                << std::endl << "P_SELECTION: " << P_SELECTION
                << std::endl << "NUM_SUB_PERIODS: " << NUM_SUB_PERIODS
                << std::endl << "RES_CUTOFF: " << RES_CUTOFF
                << std::endl << "OBJ_UB: " << OBJ_UB
                << std::endl << "WALL_SEARCH_TIME: " << WALL_SEARCH_TIME
                << std::endl << "CPU_SEARCH_TIME: " << CPU_SEARCH_TIME
                << std::endl << "CONVERGE_STOP: " << CONVERGE_STOP;

      std::cout << std::endl << std::endl;

    }
};

#endif
