#ifndef __set_obj_H__
#define __set_obj_H__

#include <debug.h>
#include <set>
#include <vector>
#include <random>
#include <numeric>

class set_obj{
  public:
    int num_elements;
    std::vector<int> idx_data;
    std::vector<int> set_data;

    set_obj(){};

    set_obj(int num_elements) : num_elements(num_elements){
      idx_data = std::vector<int>(num_elements, -1);
      set_data = std::vector<int>();
    };

    void clear(){
      idx_data = std::vector<int>(num_elements, -1);
      set_data.clear();
    };

    void fill(){
      idx_data = std::vector<int>(num_elements, 1);
      set_data = std::vector<int>(num_elements);
      std::iota (std::begin(set_data), std::end(set_data), 0);
    };

    void fill(int _num_elements){
      num_elements = _num_elements;
      idx_data = std::vector<int>(num_elements, 1);
      set_data = std::vector<int>(num_elements);
      std::iota(std::begin(set_data), std::end(set_data), 0);
    };

    void clear(int num_elements_){
      num_elements = num_elements_;
      idx_data = std::vector<int>(num_elements, -1);
      set_data.clear();
    };

    void setAll(){
      set_data.clear();
      for (int i =0; i < num_elements; ++i){
        idx_data[i] = i;
        set_data.push_back(i);
      }
    };

    std::vector<int> getSet(){
      return set_data;
    };

    bool empty(){
      return set_data.empty();
    };

    bool is_element(int idx){
      return (idx_data[idx] > -1);
    };

    int size(){
      return set_data.size();
    };

    int idx_size(){
      return idx_data.size();
    };

    int get(int idx){
      if (idx < set_data.size()){
        return set_data[idx];
      }
      else{
        return -1;
      }
    };

    int getRandomElement(std::mt19937 &rng){
      if (set_data.empty()){
        return -1;
      }
      std::uniform_int_distribution<int> uni(0,set_data.size()-1);
      int idx = uni(rng);
      return set_data[idx];
    };

    void addElement(const int idx){
      if (idx_data[idx] < 0){
        idx_data[idx] = set_data.size();
        set_data.push_back(idx);
      }
    };

    void addElements(const std::vector<int> &idxs){
      for (int i = 0; i < idxs.size(); ++i){
        if (idx_data[idxs[i]] < 0){
          idx_data[idxs[i]] = set_data.size();
          set_data.push_back(idxs[i]);   
        }
      }
    };

    void removeElement(const int idx){
      if (idx_data[idx] > -1){
        // replace current element with back element
        set_data[idx_data[idx]] = set_data.back();
        // change idx_data for back element
        idx_data[set_data.back()] = idx_data[idx];
        // remove from idx data
        idx_data[idx] = -1;
        // remove back element
        set_data.pop_back();
      }
    };

    void removeElements(const std::vector<int> &idxs){
      for (int i = 0; i < idxs.size(); ++i){
        if (idx_data[idxs[i]] > -1){
          // replace current element with back element
          set_data[idx_data[idxs[i]]] = set_data.back();
          // change idx_data for back element
          idx_data[set_data.back()] = idx_data[idxs[i]];
          // remove from idx data
          idx_data[idxs[i]] = -1;
          // remove back element
          set_data.pop_back();
        }
      }
    };

    void makeUnion(const set_obj &src){
      for (int i = 0; i < src.set_data.size(); ++i){
        addElement(src.set_data[i]);
      }
    };

    // void setMinus(const SetObj &src){
    //   for (int i = 0; i < src.set_data.size(); ++i){
    //     removeElement(src.set_data[i]);
    //   }
    // };

    // SetObj getIntersection(const SetObj &src){
    //   SetObj temp(num_elements);

    //   if (src.size() > set_data.size()){
    //     for (int i = 0; i < src.set_data.size(); ++i){
    //       if (is_element(src.set_data[i])){
    //         temp.addElement(src.set_data[i]);
    //       }
    //     }
    //   }
    //   else{
    //     for (int i = 0; i < set_data.size(); ++i){
    //       if (src.is_element(set_data[i])){
    //         temp.addElement(set_data[i]);
    //       }
    //     }
    //   }
    //   return temp;
    // };

    ~set_obj(){};
};
#endif
