#include "Instance.h"

float EPSILON = 0.00001;

void Job::print(){
	cout << dur << " ";
	for(int i=0;i<det.size();i++) cout << det[i] << " ";
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
}

void Instance::exact(){

	GRBEnv env = GRBEnv();

	// Big-M 2
	vector<double> bigM2(nbMachines,0.0);
	for(int i = 0; i<nbMachines;i++){
		double products = 1.0;
		int minD = 999;
		for(int j = 0; j<nbJobs;j++){
			minD = min(minD, jobs[j].dur);
			products *= jobs[j].det[i];
		}
		bigM2[i] = min(products, 1.0 + double(double(maintenances[i])/double(minD)));
	}
	cout << endl;

	// Model
	try{
		GRBModel model = GRBModel(env);
		int nbPositions = 2 * (nbJobs - nbMachines) + 1;
		// Declaration of the variables for the model
		vector<vector<vector<GRBVar> > > xijh(nbMachines);
		vector<vector<vector<GRBVar> > > aijh(nbMachines);
		vector<vector<GRBVar> > sih(nbMachines);
		vector<vector<GRBVar> > dih(nbMachines);
		GRBVar makespan = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBLinExpr obj = 0;

		// Initialization
		for(int i = 0; i<nbMachines;i++){	
			xijh[i].resize(nbJobs);
			aijh[i].resize(nbJobs);
			for(int j = 0; j<nbJobs;j++){
				xijh[i][j].resize(nbPositions);
				aijh[i][j].resize(nbPositions);
				for(int k = 0; k<nbPositions;k++){
					xijh[i][j][k] = model.addVar(0, 1, 0, GRB_INTEGER);
					aijh[i][j][k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
				}
			}
			sih[i].resize(nbPositions);
			dih[i].resize(nbPositions);
			for(int k = 0; k<nbPositions;k++){
				sih[i][k] = model.addVar(0, 1, 0, GRB_INTEGER);
				dih[i][k] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
			}			
		}

		model.update();

		// Objective function
		obj += makespan;
		model.setObjective(obj, GRB_MINIMIZE);

		// Constraints -- At most one activity per time slot
		for(int i = 0; i<nbMachines;i++){	
			for(int k = 0; k<nbPositions;k++){
				GRBLinExpr nbActivity = sih[i][k];
				for(int j = 0; j<nbJobs;j++){
					nbActivity += xijh[i][j][k];
				}
				model.addConstr(nbActivity <= 1);
			}
		}

		// Constraints -- Each job is scheduled
		for(int j = 0; j<nbJobs;j++){
			GRBLinExpr isScheduled = 0;
			for(int i = 0; i<nbMachines;i++){	
				for(int k = 0; k<nbPositions;k++){
					isScheduled+= xijh[i][j][k];
				}
			}
			model.addConstr(isScheduled == 1);
		}

		// Constraints -- M1
		for(int i = 0; i<nbMachines;i++){	
			for(int j = 0; j<nbJobs;j++){
				for(int k = 0; k<nbPositions;k++){
					model.addConstr(aijh[i][j][k] >= jobs[j].dur * dih[i][k] - bigM2[i] * jobs[j].dur * (1 - xijh[i][j][k]));
				}
			}
		}

		// Constraints -- Makespan grater than the durations of each machine
		for(int i = 0; i<nbMachines;i++){
			GRBLinExpr duration = 0;
			for(int j = 0; j<nbJobs;j++){
				for(int k = 0; k<nbPositions;k++){
					duration += aijh[i][j][k];
				}
			}
			for(int k = 0; k<nbPositions;k++){
				duration += sih[i][k] * maintenances[i];
			}
			model.addConstr(duration <= makespan);
		}

		// Constraints -- No free maintenance
		for(int i = 0; i<nbMachines;i++){
			for(int k = 1; k<nbPositions;k++){
				GRBLinExpr left = sih[i][k];
				GRBLinExpr right = sih[i][k-1];
				for(int j = 0; j<nbJobs;j++){
					 left += xijh[i][j][k];
					 right += xijh[i][j][k-1];
				}
				model.addConstr(left <= right);
			}
		}

		// Constraints -- Get the delta
		for(int i = 0; i<nbMachines;i++){
			for(int j = 0; j<nbJobs;j++){
				for(int k = 1; k<nbPositions;k++){
					model.addConstr(jobs[j].det[i] * dih[i][k-1] - bigM2[i] * jobs[j].det[i]  * (1 + sih[i][k-1] - xijh[i][j][k-1]) <= dih[i][k]);
				}
			}
		}

		// Constraints -- LB on delta
		for(int i = 0; i<nbMachines;i++){
			for(int k = 0; k<nbPositions;k++){
				model.addConstr(dih[i][k] >= 1.0);
			}
		}

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
			if(infos.UB == infos.LB) infos.opt = true;
			for(int i = 0; i<nbMachines;i++){	
				cout << "Machine " << i << ":" << endl;
				for(int k = 0; k<nbPositions;k++){
					if(ceil(sih[i][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						cout << " | ";
					}
					else{
						int c = 0;
						for(int j = 0; j<nbJobs;j++){
							if(ceil(xijh[i][j][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
								cout << j << " (" << aijh[i][j][k].get(GRB_DoubleAttr_X) << ") ";
								c++;
							}
						}
						if(c == 0) cout << " x ";
					}
					cout << "( " << dih[i][k].get(GRB_DoubleAttr_X) << ") " ;
				}
				cout << endl;
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
