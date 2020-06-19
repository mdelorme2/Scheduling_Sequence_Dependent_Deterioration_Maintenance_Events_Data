#include "main.h"

/*	*************************************************************************************
	*************************************  MAIN *****************************************
	************************************************************************************* */

int main(int argc, char **argv){
	
	double initTimeModelCPU = getCPUTime();
	// local variables
	Instance instance;
	string filein = argv[2];
	string path = argv[1];
	string pathAndFileout = argv[3];

	// functions
	instance.load(path, filein);
	instance.createArcs();
	instance.print();
	instance.infos.time.resize(2);
	instance.infos.time[1] = getCPUTime() - initTimeModelCPU;
	instance.exact();
	instance.infos.time[0] = getCPUTime() - initTimeModelCPU;
	instance.printInfo(pathAndFileout);

	return 0;
}


/*	*************************************************************************************
	**********************************	FUNCTIONS  **************************************
	************************************************************************************* */


