#include <stdlib.h>
#include <TStopwatch.h>
#include "UnigenQA.h"

using namespace std;

void RunUnigenQA (TString filePath = "/home/ogolosov/Desktop/analysis/mc/model_root/dcmqgsm_12.root",
								 TString qaPath = "/home/ogolosov/Desktop/analysis/mc/model_qa/botvina_12agev.root") 
{
	TStopwatch timer;
	timer.Reset();
	timer.Start();
		
  cout << "Input:" << filePath << endl;
  cout << "Output:" << qaPath << endl;
	
	qa::UnigenQA qa;
	qa.Init (filePath, "events");
	qa.Init_Histograms();
	qa.Run();
	qa.Write_Histograms(qaPath);    
	
	timer.Stop();
	printf("Real time: %f\n",timer.RealTime());
	printf("CPU time: %f\n",timer.CpuTime());
}

# ifndef __CINT__
int main (int argc, char **argv) {
		if (argc == 1) RunUnigenQA ();
		else if (argc == 2) RunUnigenQA (argv [1]);
		else RunUnigenQA (argv [1], argv [2]);
    return 0;
}
# endif