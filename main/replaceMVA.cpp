#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/TreeUtils.h"
#include "interface/AnalysisUtils.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include "TMVA/Reader.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> addMVA.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  
  //-------------------------
  // open files and get trees
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  TChain* t = new TChain("chain","the chain");
  for(unsigned int ii = 0; ii < input.size()/2; ++ii)
  {
    std::string inFileName = input.at(0+ii*2);
    std::string treeName = input.at(1+ii*2);
    t -> Add((inFileName+"/"+treeName).c_str());
  }
  long int nEntries = t->GetEntries();
  std::cout << "Added " << nEntries << " entries to the chain" << std::endl;
  
  int nJobs = opts.GetOpt<int>("Input.nJobs");
  int jobId = opts.GetOpt<int>("Input.jobId");
  long int nEntriesPerJob = int(nEntries/nJobs);
  long int firstEntry = (jobId-1)*nEntriesPerJob;
  long int lastEntry = firstEntry + nEntriesPerJob;
  if( jobId == nJobs ) lastEntry = nEntries;

  
  TreeVars treeVars;
  InitTreeVars(t,treeVars);
  
  std::map<std::string,float*> varMap;
  varMap["dipho_leadEta"] = &treeVars.dipho_leadEta;
  varMap["dipho_subleadEta"] = &treeVars.dipho_subleadEta;
  varMap["dipho_lead_ptoM"] = &treeVars.dipho_lead_ptoM;
  varMap["dipho_sublead_ptoM"] = &treeVars.dipho_sublead_ptoM;
  varMap["dipho_leadIDMVA"] = &treeVars.dipho_leadIDMVA;
  varMap["dipho_subleadIDMVA"] = &treeVars.dipho_subleadIDMVA;
  varMap["dipho_lead_sigmaEoE"] = &treeVars.dipho_lead_sigmaEoE;
  varMap["dipho_sublead_sigmaEoE"] = &treeVars.dipho_sublead_sigmaEoE;
  varMap["dipho_sigmaWV"] = &treeVars.dipho_sigmaWV;
  varMap["dipho_vtxProb"] = &treeVars.dipho_vtxProb;
  varMap["dipho_sigmaRV"] = &treeVars.dipho_sigmaRV;
  varMap["dipho_cosDeltaphi"] = &treeVars.dipho_cosDeltaphi;
  
  
  //---------------
  // clone the tree
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  TFile* outputFile = TFile::Open(outputFileName.c_str(),"RECREATE");
  outputFile -> cd();
  TTree* newTree = t -> CloneTree(0);
  
  float dipho_mva;
  t -> SetBranchAddress("dipho_mva",&dipho_mva);
  
  
  //-----
  // TMVA
  std::string diphoMVA_method = opts.GetOpt<std::string>(Form("Input.method"));
  std::string weightsFile = opts.GetOpt<std::string>(Form("Input.weightsFile"));
  std::vector<std::string> inputVariables = opts.GetOpt<std::vector<std::string> >(Form("Input.inputVariables"));
  
  TMVA::Reader* diphoMVAReader = new TMVA::Reader( "!Color:!Silent" );
  
  for(unsigned int jj = 0; jj < inputVariables.size(); ++jj)
  {
    std::string inputVariable = inputVariables.at(jj);
    
    diphoMVAReader -> AddVariable(inputVariable.c_str(),varMap[inputVariable.c_str()]);
  }
  diphoMVAReader -> BookMVA( diphoMVA_method,weightsFile.c_str() );
  
  
  //-----------------
  // loop over events
  std::cout << "jobId: " << jobId << " / " << nJobs << std::endl;
  std::cout << "tot entries: " << nEntries << std::endl;
  std::cout << "first entry: " << firstEntry << std::endl;
  std::cout << " last entry: " << lastEntry << std::endl;
  for(long int ii = firstEntry; ii < lastEntry; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << std::endl;
    t -> GetEntry(ii);
    
    // evaluate dipho MVA
    dipho_mva = diphoMVAReader -> EvaluateMVA(diphoMVA_method.c_str());
    
    newTree -> Fill();
  }
  
  
  newTree -> AutoSave();
  outputFile -> Close();
  
  
  return 0;
}
