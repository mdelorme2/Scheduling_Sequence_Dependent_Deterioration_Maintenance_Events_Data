#ifndef INSTANCE_H 
	#define INSTANCE_H
	
	#include <iostream> 
	#include <iomanip> 
	#include <fstream>
	#include <sstream>
	#include <vector>
	#include <string>
	#include <map>
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

	class Arc{
	public:
		double dur;
		int job;
		int init;
		int end;
		void print();
	};

	class Block{
	public:
		double dur;
		vector<int> jobs;
		double det;
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
		vector<vector<Arc> > arcs;
		vector<vector<double> > nodes;
		Info infos;

		void load(const string& path, const string& filein);
		void print();
		void createIndexDet();
		void createArcs();
		void exact();
		void printSolution();
		void printInfo(const string& pathAndFileout);

	};

#endif 
