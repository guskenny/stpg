// the ant class

#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <new>
#include <fstream>
#include <vector>

#define Infinity 1000000

using namespace std;

class DATA{
public:
	int machines;
	int max_power;
	int prec_no;
	int total_jobs;
	int max_time; // = compute_time_line();
	vector<int> dur;
	vector<int> rel;
	vector<int> due;
	vector<int> id;
	vector<int> power;
	vector<double> weight;
	vector<vector<int> > precedences;
	vector<int> jobs_in_machines;
	vector<int> first_job_in_machine;
	vector<int> job_belongs_to_machine;
	double lower_bound;
	
	DATA();
	DATA* copy();
	void ReadFile(ifstream &in);
	void display();
	int compute_time_line();
	void compute_lower_bound();
	~DATA();
};

#endif
