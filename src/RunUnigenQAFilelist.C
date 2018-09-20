#include <stdlib.h>
#include <TStopwatch.h>
#include <TFileCollection.h>
#include <THashList.h>
#include "UnigenQA.h"

using namespace std;

int main(int argc, char **argv) {

  const char *tree_name = "events";

  const char *file_list = nullptr;
  const char *ref_list = nullptr;
  const char *output_file = nullptr;

  if (argc == 3) {
    file_list = argv[1];
    output_file = argv[2];
  } else if (argc == 4) {
    file_list = argv[1];
    ref_list = argv[2];
    output_file = argv[3];
  } else {
    cerr << "Bad input" << endl;
    return 1;
  }

  TFileCollection fileCollection("inputFiles", "", file_list);
  fileCollection.SetDefaultTreeName(tree_name);
  fileCollection.Print();

  TChain chain(tree_name);
  chain.AddFileInfoList(fileCollection.GetList());

  qa::UnigenQA unigenQA;
  unigenQA.Init(&chain);

  TChain ref_chain(tree_name);
  if (ref_list) {
    TFileCollection refFileCollection("refFiles", "", ref_list);
    refFileCollection.SetDefaultTreeName(tree_name);
    ref_chain.AddFileInfoList(refFileCollection.GetList());
    unigenQA.SetReferenceChain(&ref_chain);
    ref_chain.ls();
  }

  unigenQA.Init_Histograms();
  unigenQA.Run();
  unigenQA.Write_Histograms(output_file);

  return 0;

}
