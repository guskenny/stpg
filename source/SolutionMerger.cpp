#include <SolutionMerger.h>

void SolutionMerger::merge(const std::vector<set_obj> &sols, std::vector<std::vector<int> > &groups, std::vector<int> &group_map){

	PE("merging " << sols.size() << " STPG solutions using merge()")

	qol::CpuTimer merge_timer;

	int nE = sols[0].num_elements;

	groups = std::vector<std::vector<int> >();
    group_map = std::vector<int>(nE,0);

    // create first group with all edges in it
    groups.push_back(std::vector<int>());
    groups[0].reserve(nE);
    for(int j = 0; j < nE; ++j){
      groups[0].push_back(j);
    }

    PF("getting groups.. ")

    // iterate over all solutions
    for (int s = 0; s < sols.size(); ++s){
    	// in case groups are added, have fixed iteration stopping point
    	int groups_size = groups.size();
    	// iterate over all groups
    	for (int g = 0; g < groups_size; ++g){
    		std::vector<int> temp_group;
    		// value for group
    		int group_val = sols[s].is_element(groups[g][0]);
    		// iterate over every element in the group to check if same still
    		for (int i = 1; i < groups[g].size(); ++i){
    			// if current value not same as group value, move to new group
    			if (group_val != sols[s].is_element(groups[g][i])){
    				// index group_map to new group
    				group_map[groups[g][i]] = groups.size();
    				// push current variable index to new group
    				temp_group.push_back(groups[g][i]);
    				// remove current variable from current group
    				// by overwriting variable with last variable in group
    				groups[g][i] = groups[g].back();
    				// remove copy of last variable in group
    				groups[g].pop_back();
    				i--; // decrement i to test newly inserted index
    			}
    		}
    		// if new group created, add to list of groups
    		if (!temp_group.empty()){
    			groups.push_back(temp_group);
    		}
    	}
    }

    PE("done!")

    if (sh.TEST_GROUPS){
	    int groups_count = 0;
	    // int total_count = 0;

	    for (int group = 0; group < groups.size(); ++group){
	      if (groups[group].size() > 1){
	        groups_count++;
	      }
	    }

	    int zero_count = 0;

	    for (int i = 0; i < nE; ++i){
	      if (group_map[i] == 0){
	        zero_count++;
	      }
	    }


	    std::cout << "groups[0].size(): " << groups[0].size() << ", zero_count: " << zero_count << std::endl;


	    int total_vars = 0;


	    std::priority_queue<std::pair<double, int>> q;
	    for (int i = 0; i < groups.size(); ++i) {
	      total_vars += groups[i].size();
	      q.push(std::pair<double, int>(groups[i].size(), i));
	    }
	    int k = 5; // number of indices we need
	    std::cout << "\nTop " << k << " group sizes:\n------------------" << std::endl;
	    for (int i = 0; i < k; ++i) {
	      int ki = q.top().second;
	      double ks = q.top().first;
	      std::cout << i+1 << ": group[" << ki << "] = " << ks << std::endl;
	      q.pop();
	    }

	    std::cout << "\ntotal_vars: " << total_vars << ", group_map.size(): " << group_map.size() << std::endl;

	    std::cout << "\nMerge time: " << merge_timer.elapsedSeconds() << std::endl << std::endl;

	    std::cout << "merge finished!" << std::endl << nE << " variables found in " << groups.size() << " groups and " << groups_count << " groups with size greater than 1"<<std::endl;
	}
}