#ifndef INSTANCE_H 
	#define INSTANCE_H
	
	#include <iostream> 
	#include <iomanip> 
	#include <fstream>
	#include <sstream>
	#include <vector>
	#include <string>
	#include <set>
	#include <math.h> 
	#include <cmath>  
	#include "gurobi_c++.h"
	#include "time.h"

	using namespace std;

	class Job{
	public:
		int dur;
		vector<double> det;
		void print();
	};

	class Info{
	public:
		bool opt;
		vector<double> time;
		double LB;
		double UB;
		int nbVar;
		int nbCons;	
		int nbNZ;
		double cont;
	};

	class Instance{
	public:
		string name;
		int nbJobs;
		int nbMachines;
		vector<Job> jobs;
		vector<int> maintenances;
		vector<vector<int> > order;
		Info infos;

		void load(const string& path, const string& filein);
		void print();
		void createBlocks();
		void exact();
		void printSolution();
		void printInfo(const string& pathAndFileout);

	};

#endif 
