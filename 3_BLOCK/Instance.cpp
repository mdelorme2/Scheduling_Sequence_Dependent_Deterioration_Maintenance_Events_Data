#include "Instance.h"

float EPSILON = 0.00001;

void Job::print(){
	cout << dur << " ";
	for(int i=0;i<det.size();i++) cout << det[i] << " ";
	cout << endl;
}

void Block::print(){
	cout << setprecision(15) << dur << " " << det << " ";
	for(int i=0;i<jobs.size();i++) cout << jobs[i] << " ";
	cout << endl;
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
/*	for(int i=0; i<nbMachines;i++){
		for(int j=0;j<blocks[i].size();j++){
			cout << i << " " << j << " ";
			blocks[i][j].print();
		}
	}*/
}

void Instance::createBlocks(){
	blocks.resize(nbMachines);
	for(int i = 0; i<nbMachines;i++){
		// Empty block
		Block b; b.dur = 0.0; b.det=1;
		blocks[i].push_back(b);
		// For each job
		for(int j=0;j<nbJobs;j++){
			vector<Block> addB;
			// For each existing block
			for(int k=0;k<blocks[i].size();k++){
				if((blocks[i][k].det - 1.0) * double(jobs[order[i][j]].dur) < double(maintenances[i])){
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
					if(isGood) addB.push_back(newB);
				}
			}
			for(int k=0;k<addB.size();k++) blocks[i].push_back(addB[k]);
		}
	}
}


void Instance::exact(){

	GRBEnv env = GRBEnv();


	// Model
	try{
		GRBModel model = GRBModel(env);

		// Declaration of the variables for the model
		vector<vector<GRBVar> > isBlockActive(nbMachines);
		vector<GRBLinExpr> isJobInBlock(nbJobs);
		GRBVar makespan = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBLinExpr obj = 0;

		// Initialization
		for(int i = 0; i<nbMachines;i++){	
			isBlockActive[i].resize(blocks[i].size());
			for(int j = 0; j<blocks[i].size();j++){
				isBlockActive[i][j] = model.addVar(0, 1, 0, GRB_INTEGER);
			}
		}

		for(int i = 0; i<nbJobs;i++){	
			isJobInBlock[i] = 0;
		}

		model.update();

		// Objective function
		obj += makespan;
		model.setObjective(obj, GRB_MINIMIZE);

		// Perform values
		for(int i = 0; i<nbMachines;i++){	
			GRBLinExpr duration = - maintenances[i];
			for(int j = 0; j<blocks[i].size();j++){
				duration += maintenances[i] * isBlockActive[i][j];
				duration += isBlockActive[i][j] * blocks[i][j].dur;
				for(int k=0;k<blocks[i][j].jobs.size();k++) isJobInBlock[blocks[i][j].jobs[k]] += isBlockActive[i][j];
			}
			// Constraints -- Makespan grater than the durations of each machine
			model.addConstr(duration <= makespan);
		}

		// Constraints -- All items need to be packed
		for(int i = 0; i<nbJobs;i++) model.addConstr(isJobInBlock[i] == 1);

		// Setting of Gurobi
		model.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);
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
				for(int j = 0; j<blocks[i].size();j++){
					if(ceil(isBlockActive[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						cout << j << " ";
						blocks[i][j].print();
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
