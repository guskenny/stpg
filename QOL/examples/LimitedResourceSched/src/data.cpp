
#include <iostream>
#include <new>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

#include "data.h"
//#include "parameters.h"

DATA::DATA(){}

void DATA::ReadFile(ifstream &in){ 
	string t;

	getline(in,t);
	//cout << t << endl;
	t="";
	in >> this->machines;
	//cout << this->machines << endl;
	getline(in,t);
	//cout << t << endl;
	t="";
	getline(in,t);
	//cout << t << endl;
	t="";
	in >> this->max_power;
	//cout << this->max_power << endl;
	//cout << "\nJid\tRelease\tProc\tDue\tPower\tWeight: " << "\n";
	int job_counter=0;
	for(int i=0;i<this->machines;i++){
		t="";
		getline(in,t);
		//cout << t << endl;
		t="";
		getline(in,t);
		//cout << t << endl;
		int jobs;
		in >> jobs;
		jobs_in_machines.push_back(jobs);
		//cout << jobs << endl;
		for(int j=0;j<jobs;j++){
			job_belongs_to_machine.push_back(i);
			string jid;
			in >> jid;
			id.push_back(job_counter);
			if(j==0) first_job_in_machine.push_back(job_counter);
			int release;
			in >> release;
			rel.push_back(release);
			int proc;
			in >> proc;
			dur.push_back(proc);
			int duet;
			in >> duet;		
			due.push_back(duet);
			int pow;
			in >> pow;
			power.push_back(pow);
			double wei;
			in >> wei;
			weight.push_back(wei);
			//cout << this->id[job_counter] << ", \t" << this->rel[job_counter] << ", \t" << this->due[job_counter] << ", \t" << this->dur[job_counter] << "\n";
			job_counter+=1;
		}
	}
	this->total_jobs = int(this->id.size());
	//cout << "\nDependencies:\n";
	t="";
	getline(in,t);
	//cout << t << endl;
	t="";
	getline(in,t);
	//cout << t << endl;
	in >> this->prec_no;
	//cout << this->prec_no << endl;
	for(int i=0;i<prec_no;i++) this->precedences.push_back(vector<int>(2,0)); // for a job that depends on another
	for(int i=0;i<prec_no;i++) {
		int job_pre;
		int job_post;
		in >> job_pre;
		this->precedences[i][0] = job_pre;
		in >> job_post;
		this->precedences[i][1] = job_post;
		//cout << this->precedences[i][0] << "\t" << this->precedences[i][1] << endl;
	}
	compute_lower_bound();
	max_time = compute_time_line();
}

DATA *DATA::copy(){
	DATA *new_sol = new DATA();
	new_sol->machines=machines;
	new_sol->max_power=max_power;
	new_sol->prec_no=prec_no;
	new_sol->total_jobs=total_jobs;
	new_sol->dur=dur;
	new_sol->rel=rel;
	new_sol->due=due;
	new_sol->id=id;
	new_sol->power=power;
	new_sol->weight=weight;
	new_sol->precedences=precedences;
	new_sol->jobs_in_machines=jobs_in_machines;
	new_sol->first_job_in_machine=first_job_in_machine;
	new_sol->job_belongs_to_machine=job_belongs_to_machine;
	new_sol->max_time = max_time;
	return new_sol;
}

/*
	returns the time line computed by considering max cumulative duration by machine pairs
*/

int DATA::compute_time_line(){
	int time_line=0;
	int max_release=0;
	for(int m = 0; m < machines; m += 2){
		int first_estimate=0;
		int second_estimate=0;		
		for(int i = first_job_in_machine[m];i<first_job_in_machine[m]+jobs_in_machines[m];i++) {
			first_estimate+=dur[i];
			if(rel[i]>max_release) max_release=rel[i];
		}
		if(m+1<machines){
			for(int i = first_job_in_machine[m+1];i<first_job_in_machine[m+1]+jobs_in_machines[m+1];i++) 
				second_estimate+=dur[i];
		}
		if(first_estimate>second_estimate) time_line += first_estimate;
		else time_line +=second_estimate;
	}
	time_line += max_release;//+time_line;
	//cout << "\nTime line: " << time_line << endl;
	return time_line;
}

void DATA::compute_lower_bound(){
	
}

void DATA::display(){
	/*cout << "\nJid\tRelease\tDue\tDur\t: " << n << "\n";
	for(int i=0;i<n;i++){
		cout << this->id[i] << ", \t" << this->rel[i] << ", \t" << this->due[i] << ", \t" << this->dur[i] << "\n";
	}

	for(int i=0;i<n;i++){
		cout << "\n";
		for(int j=0;j<n;j++) cout << this->setup[i][j] << " ";
	}*/
}

DATA::~DATA(){;}

