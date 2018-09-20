#include <stdlib.h>
#include <TStopwatch.h>
#include "UnigenQA.h"

using namespace std;

void RunUnigenQA(TString filePath = "/home/ogolosov/Desktop/analysis/mc/root/dcmqgsm_12.root",
                 TString qaPath = "/home/ogolosov/Desktop/analysis/mc/qa/qa_botvina_T0_12agev.root") {
  TStopwatch timer;
  timer.Reset();
  timer.Start();

  qa::UnigenQA qa;
  qa.Init(filePath, "events");
  qa.Init_Histograms();
  qa.Run();
  qa.Write_Histograms(qaPath);

  timer.Stop();
  printf("Real time: %f\n", timer.RealTime());
  printf("CPU time: %f\n", timer.CpuTime());
}

# ifndef __CINT__

int main(int argc, char **argv) {
  TString in;
  TString out;

  cout << argc << endl;
  if (argc < 2) in = "*.root";
  else in = argv[1];
  if (argc < 3) out = "UnigenQA.root";
  else out = argv[2];

  cout << "Input:" << in << endl;
  cout << "Output:" << out << endl;

  RunUnigenQA(in, out);
  return 0;

}
# endif