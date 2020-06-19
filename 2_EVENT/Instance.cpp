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
		double maxD = 1.0;
		for(int j = 0; j<nbJobs;j++){
			minD = min(minD, jobs[j].dur);
			maxD = max(maxD, jobs[j].det[i]);
			products *= jobs[j].det[i];
		}
		bigM2[i] = min(products,  maxD * (1.0 + double(double(maintenances[i])/double(minD))));
	}

	// Model
	try{
		GRBModel model = GRBModel(env);
		// Declaration of the variables for the model
		vector<vector<vector<GRBVar> > > xijj(nbMachines);
		vector<vector<GRBVar> > sij(nbMachines);
		vector<vector<GRBVar> > aij(nbMachines);
		vector<GRBVar> dj(nbJobs);
		GRBVar makespan = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		GRBLinExpr obj = 0;

		// Initialization
		for(int i = 0; i<nbMachines;i++){	
			xijj[i].resize(nbJobs);
			aij[i].resize(nbJobs);
			sij[i].resize(nbJobs);
			for(int j = 0; j<nbJobs;j++){
				aij[i][j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
				sij[i][j] = model.addVar(0, 1, 0, GRB_INTEGER);
				xijj[i][j].resize(nbJobs);
				for(int k = 0; k<nbJobs;k++){
					if(double(double(jobs[j].dur) / double(double(jobs[j].det[i]) - 1.0)) >= double(double(jobs[k].dur) / double(double(jobs[k].det[i]) - 1.0)))
						xijj[i][j][k] = model.addVar(0, 1, 0, GRB_INTEGER);
					else 
						xijj[i][j][k] = model.addVar(0, 0, 0, GRB_INTEGER);
				}
			}
		}

		for(int j = 0; j<nbJobs;j++){
			dj[j] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS);
		}

		model.update();

		// Objective function
		obj += makespan;
		model.setObjective(obj, GRB_MINIMIZE);

		// Constraints -- Each job is scheduled
		for(int j = 0; j<nbJobs;j++){
			GRBLinExpr isScheduled = 0;
			for(int i = 0; i<nbMachines;i++){	
				isScheduled += sij[i][j];
				for(int k = 0; k<nbJobs;k++){
					if(k != j) isScheduled+= xijj[i][k][j];
				}
			}
			model.addConstr(isScheduled == 1);
		}

		// Constraints -- Flow conservation
		for(int j = 0; j<nbJobs;j++){
			for(int i = 0; i<nbMachines;i++){
				GRBLinExpr in = 0; GRBLinExpr out = 0; 
				for(int k = 0; k<nbJobs;k++){
					if(k != j){
						out += xijj[i][j][k];
						in += xijj[i][k][j];
					}
				}
				in += sij[i][j];
				model.addConstr(in >= out);
			}
		}

		// Constraints -- M1 and M2
		for(int i = 0; i<nbMachines;i++){	
			for(int j = 0; j<nbJobs;j++){
				for(int k = 0; k<nbJobs;k++){
					model.addConstr(aij[i][k] >= jobs[k].dur * dj[k] - bigM2[i] * jobs[k].dur * (1 - xijj[i][j][k]));
					model.addConstr(dj[k] >= dj[j] * jobs[j].det[i] - bigM2[i] * jobs[j].det[i] * (1 - xijj[i][j][k]));
				}
				model.addConstr(aij[i][j] >= jobs[j].dur * sij[i][j]);
				model.addConstr(dj[j] >= 1);
			}
		}

		// Constraints -- Makespan grater than the durations of each machine
		for(int i = 0; i<nbMachines;i++){
			GRBLinExpr duration = -maintenances[i];
			for(int j = 0; j<nbJobs;j++){
				duration += aij[i][j];
				duration += sij[i][j] * maintenances[i];
			}
			model.addConstr(duration <= makespan);
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
			if(infos.UB - infos.LB < 0.001) infos.opt = true;
			for(int i = 0; i<nbMachines;i++){	
				cout << "Machine " << i << ":" << endl;
				for(int j = 0; j<nbJobs;j++){
					if(ceil(sij[i][j].get(GRB_DoubleAttr_X) - EPSILON) == 1){
						cout << j << " " << aij[i][j].get(GRB_DoubleAttr_X) << "(" << dj[j].get(GRB_DoubleAttr_X) << ") ";
						int curr = j;
						bool hbf = true;
						while(hbf){
							hbf = false;
							for(int k = 0; k<nbJobs;k++){
								if(ceil(xijj[i][curr][k].get(GRB_DoubleAttr_X) - EPSILON) == 1){
									cout << k << " -- " << aij[i][k].get(GRB_DoubleAttr_X) << "(" << dj[k].get(GRB_DoubleAttr_X) << ") --";
									curr = k;
									hbf = true;
								}
							}
						}
						cout << endl;
					}
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
