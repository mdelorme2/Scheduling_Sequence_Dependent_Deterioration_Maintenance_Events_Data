#include "Instance.h"

float EPSILON = 0.00001;

void Job::print(){
	cout << dur << " ";
	for(int i=0;i<det.size();i++) cout << det[i] << " ";
	cout << endl;
}

void Arc::print(){
	cout << setprecision(15) << dur << " " << init << " " << end << " " << job << endl;
}

void Instance::load(const string& path, const string& filein){
	
	// Local variables 
	istringstream iss;
	string parser;
	string garbage;
	string nameFile = path + filein;

	// File opening
	ifstream file(nameFile.c_str(), ios::in);

	// File lecture
	if (file){
		name = filein;
		
		// Read the number of machines 
		getline(file, parser); iss.str(parser); iss >> nbMachines; iss.clear(); 
		maintenances.resize(nbMachines);

		// Read the number of jobs 
		getline(file, parser); iss.str(parser); iss >> nbJobs; iss.clear(); 

		// Read the job durations 
		getline(file, parser); iss.str(parser); 
		for(int i=0; i<nbJobs;i++){
			Job j; 
			iss >> j.dur;
			jobs.push_back(j);
		}
		iss.clear(); 

		// Read the maintenance durations 
		getline(file, parser); iss.str(parser); 
		for(int i=0; i<nbMachines;i++) iss >> maintenances[i];
		iss.clear(); 

		// Read the job deteriorations
		for(int i=0;i<nbJobs;i++){
			jobs[i].det.resize(nbMachines);
			getline(file, parser); iss.str(parser); 
			for(int j=0; j<nbMachines;j++) iss >> jobs[i].det[j];
			iss.clear();
		}
		file.close();

		// Order
		order.resize(nbMachines);
		for(int j=0; j<nbMachines;j++){
			vector<bool> waP(nbJobs,false);
			bool cont = true;
			while (cont){
				cont = false;
				int maxIdx = -1;
				double maxVal = 0.0;
				for(int i=0;i<nbJobs;i++){
					if(waP[i] == false && (double(double(jobs[i].dur) / double(double(jobs[i].det[j]) - 1.0)) > maxVal)){
						maxVal = double(double(jobs[i].dur) / double(double(jobs[i].det[j]) - 1.0));
						maxIdx = i;
						cont = true;
					}
				}
				if(cont){
					order[j].push_back(maxIdx);
					waP[maxIdx] = true;
				}
			}
		}
	}
	else cout << "Could not open the file " << nameFile << endl;
}

void Instance::print(){
	cout << name << " " << nbJobs << " " << nbMachines << endl;
	for(int i=0; i<nbMachines;i++) cout << maintenances[i] << " ";
	cout << endl;
	for(int i=0; i<order.size();i++){
		for(int j=0;j<order[i].size();j++){
			cout << order[i][j] << "\t";
		}
		cout << endl;
	}
	for(int i=0; i<nbJobs;i++){
		cout << i << " ";
		jobs[i].print();
	}
	/*for(int i=0; i<nbMachines;i++){
		for(int j=0;j<arcs[i].size();j++){
			cout << i << " " << j << " ";
			arcs[i][j].print();
		}
	}*/
}

/*
void Instance::createArcs(){
	// Variables
	vector<vector<Block> > blocks;
	arcs.resize(nbMachines);
	nodes.resize(nbMachines);
	blocks.resize(nbMachines);
	for(int i = 0; i<nbMachines;i++){
		map<double,int> mymap;
		int idx = 0;
		// Empty block
		Block b; b.dur = 0.0; b.det=1.0;
		mymap[1.0]=idx; idx++;
		blocks[i].push_back(b);		
		// For each job
		for(int j=0;j<nbJobs;j++){
			set<vector<int> > myset;
			vector<Block> addB;
			// For each existing block
			for(int k=0;k<blocks[i].size();k++){
				Block newB = blocks[i][k];
				bool isGood = true;
				newB.dur = newB.dur + newB.det * double(jobs[order[i][j]].dur);
				newB.det = newB.det * jobs[order[i][j]].det[i];
				newB.jobs.push_back(order[i][j]);
				// Test if the block is iteresting
				for(int l = 0; l < newB.jobs.size();l++){	
					double durL = 0.0; double durR = 0.0;
					double detL = 1.0; double detR = 1.0;  
					// Get time for the left
					for(int m = 0; m <= l; m++){	
						durL += jobs[newB.jobs[m]].dur * detL;				
						detL = detL * jobs[newB.jobs[m]].det[i];
					}
					// Get time for the right
					for(int m = l+1; m < newB.jobs.size(); m++){	
						durR += jobs[newB.jobs[m]].dur * detR;				
						detR = detR * jobs[newB.jobs[m]].det[i];
					}
					if(newB.dur >= durL + durR + double(maintenances[i])) isGood = false;
				}
				// Block is valid, so add the arc
				if(isGood){
					addB.push_back(newB);
					Arc newA;
					newA.dur = blocks[i][k].det * double(jobs[order[i][j]].dur);
					if(blocks[i][k].det == 1.0) newA.dur += maintenances[i];
					newA.job = order[i][j];
					newA.init = mymap.find(blocks[i][k].det)->second;
					if(mymap.count(newB.det)<=0){
						mymap[newB.det]=idx;
						newA.end = idx;
						idx++;
					}
					else{
						newA.end = mymap.find(newB.det)->second;
					}
					vector<int> test(2); test[0] = newA.init; test[1] = newA.end;
					if(myset.count(test) == 0) {
						arcs[i].push_back(newA);
						myset.insert(test);
					}
				}	
			}
			for(int k=0;k<addB.size();k++) blocks[i].push_back(addB[k]);
		}
		for(map<double,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
			nodes[i].push_back(it->first);
		}
	}
}*/

void Instance::createArcs(){
	// Variables
	arcs.resize(nbMachines);
	nodes.resize(nbMachines);
	for(int i = 0; i<nbMachines;i++){
		vector<vector<double> > lossT(nbMachines);
		map<double,int> mymap;
		map<double,int> copyMymap;
		vector<double> dead;
		int idx = 0;
		// Empty block
		mymap[1.0]=idx; idx++; lossT[i].push_back(0);
		dead.push_back(99999.0);
		// For each job
		for(int j=0;j<nbJobs;j++){
			set<vector<int> > myset;
			copyMymap = mymap;
			// For each existing node
			for(map<double,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
				Arc newA;
				newA.dur = it->first  * double(jobs[order[i][j]].dur);
				if(it->first == 1.0) newA.dur += maintenances[i];
				newA.job = order[i][j];
				newA.init = it->second;
				// Update deadline
				double newDead = min(dead[it->second] - double(jobs[order[i][j]].dur), maintenances[i] / (it->first * jobs[order[i][j]].det[i] - 1.0));
				if(newDead > 0.0){
					if(copyMymap.count(it->first  * jobs[order[i][j]].det[i])<=0){
						copyMymap[it->first  * jobs[order[i][j]].det[i]]=idx;
						newA.end = idx;
						dead.push_back(newDead);
						idx++;
					}
					else{
						newA.end = copyMymap.find(it->first  * jobs[order[i][j]].det[i])->second;
						dead[newA.end] = max(newDead, dead[newA.end]);
					}
					vector<int> test(2); test[0] = newA.init; test[1] = newA.end;
					if(myset.count(test) == 0) {
						arcs[i].push_back(newA);
						myset.insert(test);
					}
				}
			}	
			mymap = copyMymap;
		}
		for(map<double,int>::iterator it=mymap.begin(); it!=mymap.end(); ++it){
			nodes[i].push_back(it->first);
		}
	}
}

void Instance::exact(){

	GRBEnv env = GRBEnv();

	// Model
	try{
		GRBModel model = GRBModel(env);

		// Declaration of the variables for the model
		vector<vector<GRBVar> > isArcActive(nbMachines);
		vector<GRBLinExpr> isJobTaken(nbJobs);
		vector<vector<GRBLinExpr> > cIn(nbMachines);
		vector<vector<GRBLinExpr> > cOut(nbMachines);
		vector<GRBLinExpr> timeOnMachineI(nbMachines);
		GRBVar makespan = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBLinExpr obj = 0.0;

		// Initialization
		for(int i = 0; i<nbMachines;i++){	
			isArcActive[i].resize(arcs[i].size());
			for(int j = 0; j<arcs[i].size();j++){
				isArcActive[i][j] = model.addVar(0, 1, 0, GRB_INTEGER);
			}
			cIn[i].resize(nodes[i].size());
			cOut[i].resize(nodes[i].size());
			for(int k = 0; k<nodes[i].size();k++){
				cIn[i][k] = 0.0;
				cOut[i][k] = 0.0;
			}
			timeOnMachineI[i] = - maintenances[i];
		}

		for(int i = 0; i<nbJobs;i++){	
			isJobTaken[i] = 0.0;
		}

		model.update();

		// Objective function
		obj += makespan;
		model.setObjective(obj, GRB_MINIMIZE);

		// Perform values
		for(int i = 0; i<nbMachines;i++){	
			for(int j = 0; j<arcs[i].size();j++){
				isJobTaken[arcs[i][j].job] += isArcActive[i][j];
				timeOnMachineI[i] += isArcActive[i][j] * arcs[i][j].dur;
				cIn[i][arcs[i][j].end] += isArcActive[i][j];
				cOut[i][arcs[i][j].init] += isArcActive[i][j];
			}
		}

		// Constraints -- Makespan grater than the durations of each machine
		for(int i = 0; i<nbMachines;i++)
		model.addConstr(timeOnMachineI[i] <= makespan);

		// Constraints -- All items need to be packed
		for(int i = 0; i<nbJobs;i++) model.addConstr(isJobTaken[i] == 1);

		// Constraints -- Flow conservation
		for(int i = 0; i<nbMachines;i++){
			for(int k=0;k<nodes[i].size();k++){
				if(k >= 1){
					model.addConstr(cIn[i][k] >= cOut[i][k]);
				}
			}
		}
	
		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
	//	model.getEnv().set(GRB_IntParam_Method, 2);
		model.getEnv().set(GRB_IntParam_Threads, 1);
		model.getEnv().set(GRB_DoubleParam_MIPGap, 0);
		model.getEnv().set(GRB_DoubleParam_IntFeasTol, 0.0000001);
		
		model.optimize();

		// Fill info
		infos.opt = false;
		infos.LB  =  model.get(GRB_DoubleAttr_ObjBound);
		infos.nbVar =  model.get(GRB_IntAttr_NumVars);
		infos.nbCons = model.get(GRB_IntAttr_NumConstrs);
		infos.nbNZ = model.get(GRB_IntAttr_NumNZs);

		// If no solution found
		if (model.get(GRB_IntAttr_SolCount) < 1)
			cout << "Failed to optimize ILP. " << endl;
		else{
			infos.UB = model.get(GRB_DoubleAttr_ObjVal);
			if(infos.UB - infos.LB < 0.001) infos.opt = true;
			for(int i = 0; i<nbMachines;i++){	
				cout << "Machine " << i << ":" << endl;
				for(int j = 0; j<arcs[i].size();j++){
					if(ceil(isArcActive[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						cout << j << " ";
						arcs[i][j].print();
					}
				}
			}
		}
	}
	
	// Exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}
}
	
void Instance::printInfo(const string& pathAndFileout){
	string nameFile = pathAndFileout;
	std::ofstream file(nameFile.c_str(), std::ios::out | std::ios::app);
	file << name << "\t" << infos.opt << "\t" << infos.time[0] << "\t" << infos.time[1] << "\t" << setprecision(9) << infos.LB << "\t" << infos.UB  << "\t" << infos.nbVar << "\t" << infos.nbCons << "\t" << infos.nbNZ  << endl;
	file.close();
}
