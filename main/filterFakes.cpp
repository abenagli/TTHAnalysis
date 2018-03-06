#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h" 
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TLorentzVector.h"
#include "TGraph.h"
#include "TLegend.h"



int main(int argc, char** argv)
{
  if( argc < 2 )
  {
    std::cerr << ">>>>> filterFakes.cpp::usage:   " << argv[0] << " configFileName" << std::endl;
    return -1;
  }

  
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  //--- open files and get trees
  std::vector<std::string> input = opts.GetOpt<std::vector<std::string> >("Input.input");
  
  TChain* t = new TChain("chain","the chain");
  for(unsigned int ii = 0; ii < input.size()/2; ++ii)
  {
    std::string inFileName = input.at(0+ii*2);
    std::string treeName = input.at(1+ii*2);
    std::cout << ">>> adding " << inFileName+"/"+treeName << std::endl;
    t -> Add((inFileName+"/"+treeName).c_str());
  }
  long int nEntries = t->GetEntries();
  std::cout << "Added " << nEntries << " entries to the chain" << std::endl;
  
  float dipho_leadGenMatch;
  float dipho_subleadGenMatch;
  t -> SetBranchAddress("dipho_leadGenMatch",    &dipho_leadGenMatch);
  t -> SetBranchAddress("dipho_subleadGenMatch", &dipho_subleadGenMatch);
  
  
  //-------------------------------------------------------
  //--- create a new file + a clone of old tree in new file
  std::string outputFileName = opts.GetOpt<std::string>("Output.outputFileName");
  std::string outputTreeName = opts.GetOpt<std::string>("Output.outputTreeName");
  
  int doPromptFake = opts.GetOpt<int>("Output.doPromptFake");
  int doFakeFake   = opts.GetOpt<int>("Output.doFakeFake");
  
  if( doPromptFake && doFakeFake )
  {
    std::cerr << ">>>>> filterFakes.cpp::ERROR: set either doPromptFake or doFakeFake, but not both!" << std::endl;
    return -1;
  }
  
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile -> cd();
  
  TTree* outputTree = t->CloneTree(0);
  outputTree -> SetName(outputTreeName.c_str());
  
  for(long int ii = 0; ii < nEntries; ++ii)
  {
    if( ii%1000 == 0 ) std::cout << ">>> Reading entry " << ii << " / " << nEntries << "\r" << std::flush;
    t -> GetEntry(ii);
    
    bool accepted = true;
    if( doFakeFake )   accepted = (dipho_leadGenMatch != 1 && dipho_subleadGenMatch != 1);
    if( doPromptFake ) accepted = ( !(dipho_leadGenMatch == 1 && dipho_subleadGenMatch == 1) && 
                                     (dipho_leadGenMatch == 1 || dipho_subleadGenMatch == 1) );
    
    if( !accepted ) continue;
    
    outputTree -> Fill();
  }
  std::cout << std::endl;
  
  
  outputTree -> Print();
  outputTree -> AutoSave();
  outputFile -> Close();
  
  
  return 0;
}
