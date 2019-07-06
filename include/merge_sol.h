#ifndef __merge_sol_H__
#define __merge_sol_H__

#include <debug.h>
#include <vector>
#include "set_obj.h"

class merge_sol{
  public:
    int n_edges;
    int n_nodes;
    // vector of length E + N to store mapping from variables to groups
    std::vector<int> group_map;
    // vector of set_obj to store groups
    std::vector<set_obj> groups;

    merge_sol(int n_edges, int n_nodes) : n_edges(n_edges), n_nodes(n_nodes){
      group_map = std::vector<int>(n_edges+n_nodes, -1);
    };

    void clear(){
      group_map = std::vector<int>(n_edges+n_nodes, -1);
      groups.clear();
    }

    bool map_valid(){
      if (group_map.empty()){
        return false;
      }
      for (int i = 0; i < group_map.size(); ++i){
        if (group_map[i] < 0){
          return false;
        }
      }
      return true;
    }

    bool groups_valid(){
      if (groups.empty()){
        return false;
      }
      int sum = 0;
      for (int i = 0; i < groups.size(); ++i){
        sum += groups[i].size();
      }
      return (sum == (n_edges * n_nodes));
    }

    int find_group(int var){
      return group_map[var];
    }

    int find_node_group(int var){
      return group_map[var + n_edges];
    }
    
    int size(){
      return groups.size();
    }

    const set_obj get_group(int group){
      return groups[group];
    }
    
    const set_obj& operator[] (int group) {
        return groups[group];
    }

    ~merge_sol(){};
};
#endif
